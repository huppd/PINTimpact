#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP
#define PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP


#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact {


extern "C" 
void OP_convectionDiffusionJSmoother(
    const int& dimens,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const nL,
    const int* const nU,
    const int* const SS,
    const int* const NN,
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
    const double* const phi,
          double* const phio,
    const double& mulI,
    const double& mulC,
    const double& mulL,
    const double& om );



/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
/// \note todo merge with SORSmoother or make interface
template<class OperatorT>
class ConvectionDiffusionJSmoother {

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

  const Teuchos::RCP<const OperatorT> op_;

public:

	/// \brief constructor
	///
  /// These options include the following:
	/// - "omega" - damping parameter
  /// - "numIters" - an \c int specifying the maximum number of iterations the 
	ConvectionDiffusionJSmoother(
			const Teuchos::RCP<const OperatorT>& op,
			Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
		omega_( pl->get<Scalar>("omega", 0.5 ) ),
		nIter_( pl->get("numIters", 1 ) ),
		op_(op) {}



	void apply( const FluxFieldT& wind, const DomainFieldT& x, RangeFieldT& y, Scalar mul, Scalar mulI, Scalar mulC, Scalar mulL ) const { std::cout << "not implmented\n"; }


protected:

  void applyStep( const FluxFieldT& wind, const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y ) const {

		const EField& m = y.getType();
		const With& wB = With::noB;

		x.exchange();

		applyBC( b, x, y );

		if( 3==SpaceT::sdim ) {
			for( Ordinal k=space()->begin(m,Z,wB); k<=space()->end(m,Z,wB); ++k )
				for( Ordinal j=space()->begin(m,Y,wB); j<=space()->end(m,Y,wB); ++j )
					for( Ordinal i=space()->begin(m,X,wB); i<=space()->end(m,X,wB); ++i ) {
						Scalar diag =
							op_->getMulI() 
								+ op_->getMulC() * op_->getConvSOp()->innerDiag3D(
									(*wind[0])(i,j,k),
									(*wind[1])(i,j,k),
									(*wind[2])(i,j,k), m, i, j, k )
								- op_->getMulL() * op_->getHelmOp()->innerDiag3D( m, i, j, k) ;
						assert( diag!=0 );
						y(i,j,k) = x(i,j,k) + omega_*( b(i,j,k)
							- op_->getMulI() * x(i,j,k)
							- op_->getMulC() * op_->getConvSOp()->innerStenc3D(
									(*wind[0])(i,j,k),
									(*wind[1])(i,j,k),
									(*wind[2])(i,j,k), x, i, j, k )
							+ op_->getMulL() * op_->getHelmOp()->innerStenc3D( x, m, i, j, k) ) / diag;
					}
		}
		else {
		}
		y.changed();
	}

	/// \brief implements smoothing for Dirichlet boundary conditions as identity
	/// in tangential / velocity direction or interpolation in wand normal
	/// direction
	void applyBC( const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y	) const {

		assert( b.getType()==y.getType() );
		assert( x.getType()==y.getType() );

		const Scalar& omegaBC = omega_;
		// U-field
		if( U==x.getType() ) {
			if( DirichletBC==space()->bcl(Y) ) {
				Ordinal j = space()->begin(U,Y,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i )
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
			}
			if( DirichletBC==space()->bcu(Y) ) {
				Ordinal j = space()->end(U,Y,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i )
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
			}

			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->begin(U,Z,With::B);
				for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i )
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
			}
			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->end(U,Z,With::B);
				for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i )
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
			}

			if( DirichletBC==space()->bcl(X) ) {
				Ordinal i = space()->begin(U,X,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j ) {
						y(i,j,k) = 0.;
						for( Ordinal ii=space()->dl(X); ii<=space()->du(X); ++ii )
							y(i,j,k) += space()->getInterpolateV2S()->getC( X, i+1, ii )*x(1+i+ii,j,k);
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( X, i+1, 0 );
					}
			}
			if( DirichletBC==space()->bcu(X) ) {
				Ordinal i = space()->end(U,X,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j ) {
						y(i,j,k) = 0.;
						for( Ordinal ii=space()->dl(X); ii<=space()->du(X); ++ii )
							y(i,j,k) += space()->getInterpolateV2S()->getC( X, i, ii )*x(i+ii,j,k);
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( X, i, 0 );
					}
			}
		}

		// V-field
		if( V==x.getType() ) {
			if( DirichletBC==space()->bcl(X) ) {
				Ordinal i = space()->begin(V,X,With::B);
				for( Ordinal k=space()->begin(V,Z,With::noB); k<=space()->end(V,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(V,Y,With::B); j<=space()->end(V,Y,With::B); ++j )
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
			}
			if( DirichletBC==space()->bcu(X) ) {
				Ordinal i = space()->end(V,X,With::B);
				for( Ordinal k=space()->begin(V,Z,With::noB); k<=space()->end(V,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(V,Y,With::B); j<=space()->end(V,Y,With::B); ++j )
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
			}

			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->begin(V,Z,With::B);
				for( Ordinal j=space()->begin(V,Y,With::noB); j<=space()->end(V,Y,With::noB); ++j )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
					}
			}
			if( DirichletBC==space()->bcu(Z) ) {
				Ordinal k = space()->end(V,Z,With::B);
				for( Ordinal j=space()->begin(V,Y,With::noB); j<=space()->end(V,Y,With::noB); ++j )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
					}
			}

			if( DirichletBC==space()->bcl(Y) ) {
				Ordinal j = space()->begin(V,Y,With::B);
				for( Ordinal k=space()->begin(V,Z,With::B); k<=space()->end(V,Z,With::B); ++k )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal jj=space()->dl(Y); jj<=space()->du(Y); ++jj )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Y, j+1, jj )*x(i,1+j+jj,k);
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( Y, j+1, 0 );
					}
			}
			if( DirichletBC==space()->bcu(Y) ) {
				Ordinal j = space()->end(V,Y,With::B);
				for( Ordinal k=space()->begin(V,Z,With::B); k<=space()->end(V,Z,With::B); ++k )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal jj=space()->dl(Y); jj<=space()->du(Y); ++jj )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Y, j, jj )*x(i,j+jj,k);
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( Y, j, 0 );
					}
			}
		}

		// W-field
		if( W==x.getType() ) {
			if( DirichletBC==space()->bcl(X) ) {
				Ordinal i = space()->begin(W,X,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
			}
			if( DirichletBC==space()->bcu(X) ) {
				Ordinal i = space()->end(W,X,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
			}

			if( DirichletBC==space()->bcl(Y) ) {
				Ordinal j = space()->begin(W,Y,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
					}
			}
			if( DirichletBC==space()->bcu(Y) ) {
				Ordinal j = space()->end(W,Y,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - x(i,j,k) );
					}
			}

			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->begin(W,Z,With::B);
				for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal kk=space()->dl(Z); kk<=space()->du(Z); ++kk )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Z, k+1, kk )*x(i,j,1+k+kk);
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( Z, k+1, 0 );
					}
			}
			if( DirichletBC==space()->bcu(Z) ) {
				Ordinal k = space()->end(W,Z,With::B);
				for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal kk=space()->dl(Z); kk<=space()->du(Z); ++kk )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Z, k, kk )*x(i,j,k+kk);
						y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - y(i,j,k) )/space()->getInterpolateV2S()->getC( Z, k, 0 );
					}
			}
		}
	}

public:

  void apply( const FluxFieldT& wind, const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::No ) const {

		const EField& m = y.getType();

		DomainFieldT temp( space(), true, m );

    assert( y.getType() == x.getType() );

    for( int i =0; i<SpaceT::sdim; ++i )
      assert( wind[i]->getType() == x.getType() );

    for( int vel_dir=0; vel_dir<SpaceT::sdim; ++vel_dir )
      wind[vel_dir]->exchange();

    for( int i=0; i<nIter_; ++i ) {

			applyStep( wind, x, y, temp );
			applyStep( wind, x, temp, y );

    }
  }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    op_->print();
  }


  bool hasApplyTranspose() const { return( false ); }

	constexpr const std::string getLabel() const { return( "ConvectionDiffusionJSmoother " ); };

}; // end of class ConvectionDiffusionJSmoother





} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::ConvectionDiffusionJSmoother< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::ConvectionDiffusionJSmoother< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::ConvectionDiffusionJSmoother< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::ConvectionDiffusionJSmoother< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP
