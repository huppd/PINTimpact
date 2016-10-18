#pragma once
#ifndef PIMPACT_DIVGRADO2OP_HPP
#define PIMPACT_DIVGRADO2OP_HPP


// for EV
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"




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

  using TS = const Teuchos::Tuple<Scalar*,3>;

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

  const Teuchos::RCP<const SpaceT> space_;

  TS c_;

public:

	DivGradO2Op( const Teuchos::RCP<const SpaceT>& space ): space_(space) {

		for( int dir=0; dir<3; ++dir ) {
			// allocate stencil
			Ordinal nTemp = 3*( space_->nLoc(dir) - 1 + 1 );
			c_[dir] = new Scalar[ nTemp ];

			if( dir<space->dim() )
				Op_getCDG_dir(
						space_->nGlo( dir ),
						space_->nLoc( dir ),
						space_->bl( dir ),
						space_->bu( dir ),
						space_->getBCLocal()->getBCL( dir ),
						space_->getBCLocal()->getBCU( dir ),
						space_->getCoordinatesGlobal()->getX( static_cast<ECoord>(dir), static_cast<EField>(dir) ),
						space_->getCoordinatesLocal()->getX( static_cast<ECoord>(dir), EField::S ),
						space_->getCoordinatesLocal()->getX( static_cast<ECoord>(dir), static_cast<EField>(dir) ),
						c_[dir] );
		}
	}


	void apply( const DomainFieldT& x, RangeFieldT& y, Belos::ETrans
			trans=Belos::NOTRANS ) const {

		x.exchange();

		if( 3==space()->dim() ) {
			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i )
						y.at(i,j,k) = innerStenc3D(x, i,j,k);
		}
		else {
			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i )
						y.at(i,j,k) = innerStenc2D(x, i,j,k);
		}
	}


	void computeResidual( const RangeFieldT& b, const DomainFieldT& x,
			RangeFieldT& res ) const {

		x.exchange();
		// inner stencil
		if( 3==space()->dim() ) {
			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i )
						res.at(i,j,k) = b.at(i,j,k) - innerStenc3D(x, i,j,k);
		}
		else {
			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i )
						res.at(i,j,k) = b.at(i,j,k) - innerStenc2D(x, i,j,k);
		}

		res.changed();
	}


	void applyInvDiag( const DomainFieldT& x, RangeFieldT& y ) const {

		if( 3==space()->dim() ) {
			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i ) {
						Scalar diag = std::abs( getC(X,i,0) + getC(Y,j,0) + getC(Z,k,0) );
						std::cout << i << "\t" << j << "\t" << k << "\t" << diag << "\n";
						y.at(i,j,k) = x.at(i,j,k)/diag;
					}
		}
		else {
			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i )
						y.at(i,j,k) = x.at(i,j,k)/std::abs( getC(X,i,0) + getC(Y,j,0) );
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
    for( int dir=0; dir<3; ++dir ) {
      out << "\ndir: " << dir << "\n";
      for( int i=1; i<=space_->nLoc(dir); ++i ) {
        out << "\ni: " << i << "\t(";
        for( int k=-1; k<=1; ++k ) {
					out << getC(dir,i,k) << "\t" ;
        }
        out << ")\n";
      }
      out << "\n";
    }
  }

  void print2Mat(  ) const {

    for( int dir=0; dir<3; ++dir ) {
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

protected:

	inline constexpr Scalar innerStenc3D( const DomainFieldT& x, const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(S,X) ) ||
		           (space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(S,Y) ) ||
		           (space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(S,Y) ) ;
		const bool bcZ = (space()->getBCLocal()->getBCL(Z) > 0 && k==space()->begin(S,Z) ) ||
		           (space()->getBCLocal()->getBCU(Z) > 0 && k==space()->end(S,Z) ) ;

		const Scalar& eps = 1.e-1;

		const Scalar epsX = (bcY||bcZ)?eps:1.;
		const Scalar epsY = (bcX||bcZ)?eps:1.;
		const Scalar epsZ = (bcX||bcY)?eps:1.;

		return( 
				epsX*getC(X,i,-1)*x.at(i-1,j  ,k  ) + epsX*getC(X,i,1)*x.at(i+1,j  ,k  ) +
				epsY*getC(Y,j,-1)*x.at(i  ,j-1,k  ) + epsY*getC(Y,j,1)*x.at(i  ,j+1,k  ) +
				epsZ*getC(Z,k,-1)*x.at(i  ,j  ,k-1) + epsZ*getC(Z,k,1)*x.at(i  ,j  ,k+1) +
				( epsX*getC(X,i,0) + epsY*getC(Y,j,0) + epsZ*getC(Z,k,0) )*x.at(i,j,k)
				);
	}

	inline constexpr Scalar innerStenc2D( const DomainFieldT& x, const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(S,X) ) ||
		           (space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(S,Y) ) ||
		           (space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(S,Y) ) ;

		const Scalar& eps = 1.e-4;

		const Scalar epsX = (bcY)?eps:1.;
		const Scalar epsY = (bcX)?eps:1.;

		return( 
				epsX*getC(X,i,-1)*x.at(i-1,j  ,k  ) + epsX*getC(X,i,1)*x.at(i+1,j  ,k  ) +
				epsY*getC(Y,j,-1)*x.at(i  ,j-1,k  ) + epsY*getC(Y,j,1)*x.at(i  ,j+1,k  ) +
				( epsX*getC(X,i,0) + epsY*getC(Y,j,0) )*x.at(i,j,k)
				);
	}

public:

	/// \name getters
	/// @{ 

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	inline constexpr const Scalar* getC( const ECoord& dir) const  {
		return( getC( static_cast<const int&>(dir) ) );
  }

  inline constexpr const Scalar* getC( const int& dir) const  {
		return( c_[dir] );
  }

	inline constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const  {
		return( getC( static_cast<const int&>(dir), i, off ) );
  }

	inline constexpr const Scalar& getC( const int& dir, Ordinal i, Ordinal off ) const  {
		return( c_[dir][ off + 1 + (i-1)*3 ] );
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
