#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP
#define PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP


#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact {


extern "C"
void OP_convectionDiffusionSOR(
    const int& dimens,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const nL,
    const int* const nU,
    const int* const SS,
    const int* const NN,
    const short int* const dir,
    const short int* const loopOrder,
    const double* const c1D,
    const double* const c2D,
    const double* const c3D,
    const double* const c1U,
    const double* const c2U,
    const double* const c3U,
    const double* const c11,
    const double* const c22,
    const double* const c33,
    const double* const phiU,
    const double* const phiV,
    const double* const phiW,
    const double* const b,
    double* const phi,
    const double& mulI,
    const double& mulC,
    const double& mulL,
    const double& om );



/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
template<class OperatorT>
class ConvectionDiffusionSORSmoother {

public:

  using SpaceT = typename OperatorT::SpaceT;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using FluxFieldT = Teuchos::Tuple< Teuchos::RCP< ScalarField<SpaceT> >, 3 >;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

  Scalar omega_;
  int nIter_;

  int ordering_;

  Teuchos::Tuple<short int,3> dirs_;

  Teuchos::Tuple<short int,3> loopOrder_;


  const Teuchos::RCP<const OperatorT> op_;

public:

  ConvectionDiffusionSORSmoother(
      const Teuchos::RCP<const OperatorT>& op,
      Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
        omega_( pl->get("omega", 1. ) ),
        nIter_( pl->get("numIters", 1 ) ),
        ordering_( pl->get("Ordering",1 ) ),
        loopOrder_( Teuchos::tuple<short int>(1,2,3) ),
        //space_(op->space_),
        op_(op) {

    if( 4==SpaceT::dimNC )
      if( 0==op->space()->rankST() )
        std::cout << "Warning!!! ConvectionDiffusionSORSmoother strange behavior for dimNC=4, problems at outflow\n";

    if( 0==ordering_ ) {
      dirs_[0] = pl->get<short int>( "dir X", 1 );
      dirs_[1] = pl->get<short int>( "dir Y", 1 );
      dirs_[2] = pl->get<short int>( "dir Z", 1 );
    }
  }



  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mulI, Scalar mulC, Scalar mulL, const Add& add=Add::No ) const { std::cout << "not implmented\n"; }

  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, const Add& add=Add::No ) const {

    // testing field consistency
    assert( z.getType() == y.getType() );

    for( int i=0; i<SpaceT::sdim; ++i )
      assert( x[i]->getType() == y.getType() );

    // exchange wind and "rhs"
    for( int vel_dir=0; vel_dir<SpaceT::sdim; ++vel_dir )
      x[vel_dir]->exchange();

    for( int i=0; i<nIter_; ++i ) {

      if( ordering_==0 )
        apply(x,y,z,dirs_,loopOrder_ );
      else
        applyNPoint( x, y, z );

    }
  }

protected:

	void applyNPoint( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z ) const {

		Teuchos::Tuple<int,SpaceT::dimension> ib = space()->getProcGrid()->getNP();
	
		Teuchos::Tuple<int,3>	dirS;
		Teuchos::Tuple<int,3>	inc;

		for( int i=0; i<3; ++i ) {
			if( (ib[i]-1)%2 == 0 ) {
			 dirS[i] = -1; 
			 inc[i] = 2;
			}
			else {
			 dirS[i] = 1; 
			 inc[i] = -2;
			}
		}

		dirs_[0] = -1;
		if( 3==SpaceT::sdim )
			for( dirs_[2]=dirS[2]; std::abs(dirs_[2])<=1; dirs_[2]+=inc[2] )
				for( dirs_[1]=dirS[1]; std::abs(dirs_[1])<=1; dirs_[1]+=inc[1] )
//					for( dirs_[0]=dirS[0]; std::abs(dirs_[0])<=1; dirs_[0]+=inc[0] )
						apply( x, y, z, dirs_, loopOrder_ );
		else {
			dirs_[2] = 1 ;
			for( dirs_[1]=-1; dirs_[1]<2; dirs_[1]+=2 )
				for( dirs_[0]=-1; dirs_[0]<2; dirs_[0]+=2 )
					apply( x, y, z, dirs_, loopOrder_ );
		}
	}


  /// \brief little helper
  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z,
      const Teuchos::Tuple<short int,3>& dirs,
      const Teuchos::Tuple<short int,3>& loopOrder ) const {

		const int sdim = SpaceT::sdim;

    z.exchange();

		applyBC( y, z);

    OP_convectionDiffusionSOR(
        sdim,
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->nl(),
        space()->nu(),
        space()->sInd(z.getType()),
        space()->eInd(z.getType()),
        dirs.getRawPtr(),
        loopOrder.getRawPtr(),
        op_->getConvSOp()->getCD( X, z.getType() ),
        op_->getConvSOp()->getCD( Y, z.getType() ),
        op_->getConvSOp()->getCD( Z, z.getType() ),
        op_->getConvSOp()->getCU( X, z.getType() ),
        op_->getConvSOp()->getCU( Y, z.getType() ),
        op_->getConvSOp()->getCU( Z, z.getType() ),
        op_->getHelmOp()->getC( X, z.getType() ),
        op_->getHelmOp()->getC( Y, z.getType() ),
        op_->getHelmOp()->getC( Z, z.getType() ),
        x[X]->getConstRawPtr(),
        x[Y]->getConstRawPtr(),
        x[Z]->getConstRawPtr(),
        y.getConstRawPtr(),
        z.getRawPtr(),
        op_->getMulI(),
        op_->getMulC(),
        op_->getMulL(),
        omega_ );

    z.changed();
  }


	/// \brief implements smoothing for Dirichlet boundary conditions as identity
	/// in tangential / velocity direction or interpolation in wand normal
	/// direction
	void applyBC( const DomainFieldT& b, RangeFieldT& y	) const {

		assert( b.getType()==y.getType() );

		const Scalar& omegaBC = omega_;
		// U-field
		if( U==y.getType() ) {
			if( DirichletBC==space()->bcl(Y) ) {
				Ordinal j = space()->begin(U,Y,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i )
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
			}
			if( DirichletBC==space()->bcu(Y) ) {
				Ordinal j = space()->end(U,Y,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i )
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
			}

			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->begin(U,Z,With::B);
				for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i )
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
			}
			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->end(U,Z,With::B);
				for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i )
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
			}

			if( DirichletBC==space()->bcl(X) ) {
				Ordinal i = space()->begin(U,X,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j ) {
						y(i,j,k) = 0.;
						for( Ordinal ii=space()->dl(X); ii<=space()->du(X); ++ii )
							y(i,j,k) += space()->getInterpolateV2S()->getC( X, i+1, ii )*y(1+i+ii,j,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( X, i+1, 0 );
					}
			}
			if( DirichletBC==space()->bcu(X) ) {
				Ordinal i = space()->end(U,X,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j ) {
						y(i,j,k) = 0.;
						for( Ordinal ii=space()->dl(X); ii<=space()->du(X); ++ii )
							y(i,j,k) += space()->getInterpolateV2S()->getC( X, i, ii )*y(i+ii,j,k);
						y(i,j,k) = y(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( X, i, 0 );
					}
			}
		}

		// V-field
		if( V==y.getType() ) {
			if( DirichletBC==space()->bcl(X) ) {
				Ordinal i = space()->begin(V,X,With::B);
				for( Ordinal k=space()->begin(V,Z,With::noB); k<=space()->end(V,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(V,Y,With::B); j<=space()->end(V,Y,With::B); ++j )
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
			}
			if( DirichletBC==space()->bcu(X) ) {
				Ordinal i = space()->end(V,X,With::B);
				for( Ordinal k=space()->begin(V,Z,With::noB); k<=space()->end(V,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(V,Y,With::B); j<=space()->end(V,Y,With::B); ++j )
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
			}

			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->begin(V,Z,With::B);
				for( Ordinal j=space()->begin(V,Y,With::noB); j<=space()->end(V,Y,With::noB); ++j )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
					}
			}
			if( DirichletBC==space()->bcu(Z) ) {
				Ordinal k = space()->end(V,Z,With::B);
				for( Ordinal j=space()->begin(V,Y,With::noB); j<=space()->end(V,Y,With::noB); ++j )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
					}
			}

			if( DirichletBC==space()->bcl(Y) ) {
				Ordinal j = space()->begin(V,Y,With::B);
				for( Ordinal k=space()->begin(V,Z,With::B); k<=space()->end(V,Z,With::B); ++k )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal jj=space()->dl(Y); jj<=space()->du(Y); ++jj )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Y, j+1, jj )*y(i,1+j+jj,k);
						y(i,j,k) = y(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( Y, j+1, 0 );
					}
			}
			if( DirichletBC==space()->bcu(Y) ) {
				Ordinal j = space()->end(V,Y,With::B);
				for( Ordinal k=space()->begin(V,Z,With::B); k<=space()->end(V,Z,With::B); ++k )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal jj=space()->dl(Y); jj<=space()->du(Y); ++jj )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Y, j, jj )*y(i,j+jj,k);
						y(i,j,k) = y(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( Y, j, 0 );
					}
			}
		}

		// W-field
		if( W==y.getType() ) {
			if( DirichletBC==space()->bcl(X) ) {
				Ordinal i = space()->begin(W,X,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
			}
			if( DirichletBC==space()->bcu(X) ) {
				Ordinal i = space()->end(W,X,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
			}

			if( DirichletBC==space()->bcl(Y) ) {
				Ordinal j = space()->begin(W,Y,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
					}
			}
			if( DirichletBC==space()->bcu(Y) ) {
				Ordinal j = space()->end(W,Y,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) += omegaBC*( b(i,j,k) - y(i,j,k) );
					}
			}

			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->begin(W,Z,With::B);
				for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal kk=space()->dl(Z); kk<=space()->du(Z); ++kk )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Z, k+1, kk )*y(i,j,1+k+kk);
						y(i,j,k) = y(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( Z, k+1, 0 );
					}
			}
			if( DirichletBC==space()->bcu(Z) ) {
				Ordinal k = space()->end(W,Z,With::B);
				for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal kk=space()->dl(Z); kk<=space()->du(Z); ++kk )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Z, k, kk )*y(i,j,k+kk);
						y(i,j,k) = y(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( Z, k, 0 );
					}
			}
		}
	}

public:

  void print( std::ostream& out=std::cout ) const {

    out << "--- " << getLabel() << "---\n";
    op_->print();
  }


  bool hasApplyTranspose() const { return( false ); }

  constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); }

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}
   

	constexpr const std::string getLabel() const { return( "ConvectionDiffusionSORSmoother" ); };


}; // end of class ConvectionDiffusionSORSmoother



} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::ConvectionDiffusionSORSmoother< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::ConvectionDiffusionSORSmoother< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::ConvectionDiffusionSORSmoother< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::ConvectionDiffusionSORSmoother< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP
