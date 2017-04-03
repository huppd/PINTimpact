#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP
#define PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP


#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"




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

  using FluxFieldT = ScalarField<SpaceT>[3];

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

	using SW = typename SpaceT::SW;

  ST omega_;
  int nIter_;

  int ordering_;

  Teuchos::Tuple<short int,3> dirs_;

  Teuchos::Tuple<short int,3> loopOrder_;


  const Teuchos::RCP<const OperatorT> op_;

	constexpr const ST& getHC( const ECoord& dir, const F& ftype, OT i, OT ii ) {
		return( op_->getHelmOp()->getC(dir,ftype,i,ii) );
	}

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



  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, ST mulI, ST mulC, ST mulL, const Add& add=Add::N ) const { std::cout << "not implmented\n"; }

  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, const Add& add=Add::N ) const {

    // testing field consistency
    assert( z.getType() == y.getType() );

    for( int i=0; i<SpaceT::sdim; ++i )
      assert( x[i].getType() == y.getType() );

    // exchange wind and "rhs"
    for( int vel_dir=0; vel_dir<SpaceT::sdim; ++vel_dir )
      x[vel_dir].exchange();

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
					for( dirs_[0]=dirS[0]; std::abs(dirs_[0])<=1; dirs_[0]+=inc[0] )
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

		applyBC( y, z );
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
        x[X].getConstRawPtr(),
        x[Y].getConstRawPtr(),
        x[Z].getConstRawPtr(),
        y.getConstRawPtr(),
        z.getRawPtr(),
        op_->getMulI(),
        op_->getMulC(),
        op_->getMulL(),
        omega_ );

		applyBC( y, z );
    z.changed();
  }


	/// \brief implements smoothing for Dirichlet boundary conditions as identity
	/// in tangential / velocity direction or interpolation in wand normal
	/// direction
	void applyBC( const DomainFieldT& b, RangeFieldT& y	) const {

		assert( b.getType()==y.getType() );

		const F& f = y.getType();

		const ST& omegaBC = omega_;
		//const ST& omegaBC = 0.9;
	
		// U-field
		if( F::U==f ) {

			// tangential direction: Y
			if( 0<space()->bcl(Y) ) {
				OT j = space()->si(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
						ST temp = 0.;
						for( OT jj=0; jj<=SW::BU(Y); ++jj )
							temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
					}
			}
			if( 0<space()->bcu(Y) ) {
				OT j = space()->ei(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
						ST temp = 0.;
						for( OT jj=SW::BL(Y); jj<=0; ++jj )
							temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
					}
			}

			// tangential direction: Z
			if( 0<space()->bcl(Z) ) {
				OT k = space()->si(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
					for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
						ST temp = 0.;
						for( OT kk=0; kk<=SW::BU(Z); ++kk )
							temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
					}
			}
			if( 0<space()->bcu(Z) ) {
				OT k = space()->ei(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
					for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
						ST temp = 0.;
						for( OT kk=SW::BL(Z); kk<=0; ++kk )
							temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
					}
			}

			// normal direction: X
			if( 0<space()->bcl(X) ) {
				OT i = space()->si(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
						ST temp = 0.;
						for( OT ii=0; ii<=SW::BU(X); ++ii )
							temp += getHC(X,f,i,ii)*y(i+ii,j,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
					}
			}
			if( 0<space()->bcu(X) ) {
				OT i = space()->ei(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
						ST temp = 0.;
						for( OT ii=SW::BL(X); ii<=0; ++ii )
							temp += getHC(X,f,i,ii)*y(i+ii,j,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
					}
			}
		}

		// V-field
		if( F::V==f ) {

			// tangential direction: X
			if( 0<space()->bcl(X) ) {
				OT i = space()->si(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j ) {
						ST temp = 0.;
						for( OT ii=0; ii<=SW::BU(X); ++ii )
							temp += getHC(X,f,i,ii)*y(i+ii,j,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
					}
			}
			if( 0<space()->bcu(X) ) {
				OT i = space()->ei(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j ) {
						ST temp = 0.;
						for( OT ii=SW::BL(X); ii<=0; ++ii )
							temp += getHC(X,f,i,ii)*y(i+ii,j,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
					}
			}

			// tangential direction: Z
			if( 0<space()->bcl(Z) ) {
				OT k = space()->si(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						ST temp = 0.;
						for( OT kk=0; kk<=SW::BU(Z); ++kk )
							temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
					}
			}
			if( 0<space()->bcu(Z) ) {
				OT k = space()->ei(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						ST temp = 0.;
						for( OT kk=SW::BL(Z); kk<=0; ++kk )
							temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
					}
			}

			// normal direction: Y
			if( 0<space()->bcl(Y) ) {
				OT j = space()->si(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						ST temp = 0.;
						for( OT jj=0; jj<=SW::BU(Y); ++jj )
							temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
					}
			}
			if( 0<space()->bcu(Y) ) {
				OT j = space()->ei(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						ST temp = 0.;
						for( OT jj=SW::DL(Y); jj<=SW::DU(Y); ++jj )
							temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
					}
			}
		}

		// W-field
		if( F::W==f ) {

			// tangential direction: X
			if( 0<space()->bcl(X) ) {
				OT i = space()->si(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
					for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
						ST temp = 0.;
						for( OT ii=0; ii<=SW::BU(X); ++ii )
							temp += getHC(X,f,i,ii)*y(i+ii,j,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
					}
			}
			if( 0<space()->bcu(X) ) {
				OT i = space()->ei(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
					for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
						ST temp = 0.;
						for( OT ii=SW::BL(X); ii<=0; ++ii )
							temp += getHC(X,f,i,ii)*y(i+ii,j,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
					}
			}

			// tangential direction: Y
			if( 0<space()->bcl(Y) ) {
				OT j = space()->si(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						ST temp = 0.;
						for( OT jj=0; jj<=SW::BU(Y); ++jj )
							temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
					}
			}
			if( 0<space()->bcu(Y) ) {
				OT j = space()->ei(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						ST temp = 0.;
						for( OT jj=SW::BL(Y); jj<=0; ++jj )
							temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
					}
			}

			// normal direction: Z
			if( 0<space()->bcl(Z) ) {
				OT k = space()->si(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						ST temp = 0.;
						for( OT kk=0; kk<=SW::BU(Z); ++kk )
							temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
					}
			}
			if( 0<space()->bcu(Z) ) {
				OT k = space()->ei(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						ST temp = 0.;
						for( OT kk=SW::BL(Z); kk<=0; ++kk )
							temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
						y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
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
