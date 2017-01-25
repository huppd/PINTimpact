#pragma once
#ifndef PIMPACT_DIVGRADO2OP_HPP
#define PIMPACT_DIVGRADO2OP_HPP


// for EV
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact{


extern "C" void Op_getCDG_dir(
    const int& M,
    const int& N,
    const int& BL,
    const int& BU,
    const int& BCL,
    const int& BCU,
    const double* const y1u,
    const double* const x1p,
    const double* const x1u,
    double* const cdg1 );



/// \brief "Laplace" for pressure 2nd Order.
///
/// independent of \c StencilWidths
/// \ingroup BaseOperator
template<class ST>
class DivGradO2Op {

public:

  using SpaceT = ST;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

	using Stenc = Stencil< Scalar, Ordinal, 1, -1, 1 >;
  using TS = const Teuchos::Tuple< Stenc, ST::sdim >;

  const Teuchos::RCP<const SpaceT> space_;

  TS c_;

public:

	DivGradO2Op( const Teuchos::RCP<const SpaceT>& space ): space_(space) {

		for( int dir=0; dir<ST::sdim; ++dir ) {
			// allocate stencil
			c_[dir] = Stenc( space_->nLoc(dir) );

			Op_getCDG_dir(
					space_->nGlo( dir ),
					space_->nLoc( dir ),
					space_->bl( dir ),
					space_->bu( dir ),
					space_->getBCLocal()->getBCL( dir ),
					space_->getBCLocal()->getBCU( dir ),
					space_->getCoordinatesGlobal()->getX( static_cast<F>(dir), static_cast<ECoord>(dir) ),
					space_->getCoordinatesLocal()->getX( F::S, static_cast<ECoord>(dir) ),
					space_->getCoordinatesLocal()->getX( static_cast<F>(dir), static_cast<ECoord>(dir) ),
					c_[dir].get() );
		}
	}


	void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N ) const {

		x.exchange();

		if( 3==SpaceT::sdim ) {
			for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
				for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
					for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i ) {
						if( Add::N==add ) y(i,j,k) = 0.;
						y(i,j,k) += innerStenc3D(x, i,j,k);
					}
		}
		else {
			for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
				for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
					for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i ) {
						if( Add::N==add ) y(i,j,k) = 0.;
						y(i,j,k) += innerStenc2D(x, i,j,k);
					}
		}
	}


	void computeResidual( const RangeFieldT& b, const DomainFieldT& x,
			RangeFieldT& res ) const {

		x.exchange();
		// inner stencil
		if( 3==SpaceT::sdim ) {
			for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
				for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
					for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i )
						res(i,j,k) = b(i,j,k) - innerStenc3D(x, i,j,k);
		}
		else {
			for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
				for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
					for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i )
						res(i,j,k) = b(i,j,k) - innerStenc2D(x, i,j,k);
		}

		res.changed();
	}


	void applyInvDiag( const DomainFieldT& x, RangeFieldT& y ) const {

		const Scalar& eps = 0.1;

		if( 3==SpaceT::sdim ) {
			for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
				for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
					for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i ) {

						const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(F::S,X) ) ||
							(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(F::S,X) ) ;
						const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(F::S,Y) ) ||
							(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(F::S,Y) ) ;
						const bool bcZ = (space()->getBCLocal()->getBCL(Z) > 0 && k==space()->begin(F::S,Z) ) ||
							(               space()->getBCLocal()->getBCU(Z) > 0 && k==space()->end(F::S,Z) ) ;

						const Scalar epsX = ( (bcY||bcZ)?eps:1. );
						const Scalar epsY = ( (bcX||bcZ)?eps:1. );
						const Scalar epsZ = ( (bcX||bcY)?eps:1. );

						Scalar diag = std::fabs( epsX*getC(X,i,0) + epsY*getC(Y,j,0) + epsZ*getC(Z,k,0) );
						y(i,j,k) = x(i,j,k)/diag;
					}
		}
		else {
			for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
				for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
					for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i ) {
						const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(F::S,X) ) ||
							(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(F::S,X) ) ;
						const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(F::S,Y) ) ||
							(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(F::S,Y) ) ;

						const Scalar epsX = ( bcY?eps:1. );
						const Scalar epsY = ( bcX?eps:1. );

						y(i,j,k) = x(i,j,k)/std::fabs( epsX*getC(X,i,0) + epsY*getC(Y,j,0) );
					}
		}

		y.changed();
	}

	/// \name setter
	/// @{ 

  void assignField ( const DomainFieldT& mv ) const {};

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	///  @} 

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << " --- stencil: ---";
    for( int dir=0; dir<ST::sdim; ++dir ) {
      out << "\ndir: " << dir << "\n";
			c_[dir].print( out );
    }
  }

  void print2Mat(  ) const {

    for( int dir=0; dir<ST::sdim; ++dir ) {
			std::string fn = "A_" + toString( static_cast<ECoord>(dir) ) + "_" + std::to_string(space_->nLoc(dir)) + ".txt";

			Teuchos::RCP<std::ostream> out = Pimpact::createOstream( fn );
      for( int i=1; i<=space_->nLoc(dir); ++i ) {
        for( int k=-1; k<=1; ++k ) {
					*out << getC(dir,i,k) << "\t" ;
        }
        *out << "\n";
      }
    }
  }


	constexpr Scalar innerStenc3D( const DomainFieldT& x, const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(F::S,X) ) ||
			(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(F::S,Y) ) ||
			(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(F::S,Y) ) ;
		const bool bcZ = (space()->getBCLocal()->getBCL(Z) > 0 && k==space()->begin(F::S,Z) ) ||
			(               space()->getBCLocal()->getBCU(Z) > 0 && k==space()->end(F::S,Z) ) ;

		const Scalar& eps = 1.e-1;

		const Scalar epsX = ( (bcY||bcZ)?eps:1. );
		const Scalar epsY = ( (bcX||bcZ)?eps:1. );
		const Scalar epsZ = ( (bcX||bcY)?eps:1. );

		return( 
				epsX*getC(X,i,-1)*x(i-1,j  ,k  ) + epsX*getC(X,i,1)*x(i+1,j  ,k  ) +
				epsY*getC(Y,j,-1)*x(i  ,j-1,k  ) + epsY*getC(Y,j,1)*x(i  ,j+1,k  ) +
				epsZ*getC(Z,k,-1)*x(i  ,j  ,k-1) + epsZ*getC(Z,k,1)*x(i  ,j  ,k+1) +
				( epsX*getC(X,i,0) + epsY*getC(Y,j,0) + epsZ*getC(Z,k,0) )*x(i,j,k)
				);
	}

	constexpr Scalar innerStenc2D( const DomainFieldT& x, const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(F::S,X) ) ||
		           (space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(F::S,Y) ) ||
		           (space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(F::S,Y) ) ;

		const Scalar& eps = 1.e-1;

		const Scalar epsX = (bcY)?eps:1.;
		const Scalar epsY = (bcX)?eps:1.;

		return( 
				epsX*getC(X,i,-1)*x(i-1,j  ,k  ) + epsX*getC(X,i,1)*x(i+1,j  ,k  ) +
				epsY*getC(Y,j,-1)*x(i  ,j-1,k  ) + epsY*getC(Y,j,1)*x(i  ,j+1,k  ) +
				( epsX*getC(X,i,0) + epsY*getC(Y,j,0) )*x(i,j,k)
				);
	}

	constexpr Scalar innerDiag3D( const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(F::S,X) ) ||
			(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(F::S,Y) ) ||
			(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(F::S,Y) ) ;
		const bool bcZ = (space()->getBCLocal()->getBCL(Z) > 0 && k==space()->begin(F::S,Z) ) ||
			(               space()->getBCLocal()->getBCU(Z) > 0 && k==space()->end(F::S,Z) ) ;

		const Scalar& eps = 1.e-1;

		const Scalar epsX = ( (bcY||bcZ)?eps:1. );
		const Scalar epsY = ( (bcX||bcZ)?eps:1. );
		const Scalar epsZ = ( (bcX||bcY)?eps:1. );

		return( epsX*getC(X,i,0) + epsY*getC(Y,j,0) + epsZ*getC(Z,k,0) );
	}

	constexpr Scalar innerDiag2D( const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(F::S,X) ) ||
			(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(F::S,Y) ) ||
			(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(F::S,Y) ) ;

		const Scalar& eps = 1.e-1;

		const Scalar epsX = ( bcY?eps:1. );
		const Scalar epsY = ( bcX?eps:1. );

		return( epsX*getC(X,i,0) + epsY*getC(Y,j,0) );
	}

	/// \name getters
	/// @{ 

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	constexpr const Scalar* getC( const ECoord& dir) const  {
		return( getC( static_cast<const int&>(dir) ) );
  }

  constexpr const Scalar* getC( const int& dir) const  {
		return( c_[dir].get() );
  }

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const  {
		return( getC( static_cast<const int&>(dir), i, off ) );
  }

	constexpr const Scalar& getC( const int& dir, Ordinal i, Ordinal off ) const  {
		return( c_[dir](i,off) );
  }

	const std::string getLabel() const { return( "DivGradO2" ); };

	///  @} 


}; // end of class DivGradO2Op





} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_DIVGRADO2OP_HPP
