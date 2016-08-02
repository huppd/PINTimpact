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


extern "C" {

void Op_getCDG(
    const int& dimens,
    const int* const M,
    const int* const N,
    const int* const BL,
    const int* const BU,
    const int* const BCL,
    const int* const BCU,
    const double* const y1u,
    const double* const y2v,
    const double* const y3w,
    const double* const x1p,
    const double* const x2p,
    const double* const x3p,
    const double* const x1u,
    const double* const x2v,
    const double* const x3w,
    double* const cdg1,
    double* const cdg2,
    double* const cdg3 );

void Op_getCDG_dir(
    const int& dir,
    const int& dimens,
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

}



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
			Op_getCDG_dir(
					dir+1,
					space_->dim(),
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

		//Op_getCDG(
				//space_->dim(),
				//space_->nGlo(),
				//space_->nLoc(),
				//space_->bl(),
				//space_->bu(),
				//space_->getBCLocal()->getBCL(),
				//space_->getBCLocal()->getBCU(),
				//space_->getCoordinatesGlobal()->getX( ECoord::X, EField::U ),
				//space_->getCoordinatesGlobal()->getX( ECoord::Y, EField::V ),
				//space_->getCoordinatesGlobal()->getX( ECoord::Z, EField::W ),
				//space_->getCoordinatesLocal()->getX( ECoord::X, EField::S ),
				//space_->getCoordinatesLocal()->getX( ECoord::Y, EField::S ),
				//space_->getCoordinatesLocal()->getX( ECoord::Z, EField::S ),
				//space_->getCoordinatesLocal()->getX( ECoord::X, EField::U ),
				//space_->getCoordinatesLocal()->getX( ECoord::Y, EField::V ),
				//space_->getCoordinatesLocal()->getX( ECoord::Z, EField::W ),
				//c_[0],
				//c_[1],
				//c_[2] );
	}


	void apply( const DomainFieldT& x, RangeFieldT& y, Belos::ETrans
			trans=Belos::NOTRANS ) const {

		x.exchange();

		if( 3==space()->dim() ) {
			for( Ordinal k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( Ordinal j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( Ordinal i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i )
						y.at(i,j,k) = innerStenc3D(x, i,j,k);
		}
		else {
			for( Ordinal k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( Ordinal j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( Ordinal i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i )
						y.at(i,j,k) = innerStenc2D(x, i,j,k);
		}

		y.changed();
	}


	void computeResidual( const RangeFieldT& b, const DomainFieldT& x,
			RangeFieldT& res ) const {

		x.exchange();
		// inner stencil
		if( 3==space()->dim() ) {
			for( Ordinal k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( Ordinal j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( Ordinal i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i )
						res.at(i,j,k) = b.at(i,j,k) - innerStenc3D(x, i,j,k);
		}
		else {
			for( Ordinal k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( Ordinal j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( Ordinal i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i )
						res.at(i,j,k) = b.at(i,j,k) - innerStenc2D(x, i,j,k);
		}

		res.changed();
	}


	void applyInvDiag( const DomainFieldT& x, RangeFieldT& y ) const {

		if( 3==space()->dim() ) {
			for( Ordinal k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( Ordinal j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( Ordinal i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i ) {
						Scalar diag = std::abs( getC(X,i,0) + getC(Y,j,0) + getC(Z,k,0) );
						std::cout << i << "\t" << j << "\t" << k << "\t" << diag << "\n";
						y.at(i,j,k) = x.at(i,j,k)/diag;
					}
		}
		else {
			for( Ordinal k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( Ordinal j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( Ordinal i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i )
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

		return( 
				getC(X,i,-1)*x.at(i-1,j  ,k  ) + getC(X,i,1)*x.at(i+1,j  ,k  ) +
				getC(Y,j,-1)*x.at(i  ,j-1,k  ) + getC(Y,j,1)*x.at(i  ,j+1,k  ) +
				getC(Z,k,-1)*x.at(i  ,j  ,k-1) + getC(Z,k,1)*x.at(i  ,j  ,k+1) +
				( getC(X,i,0) + getC(Y,j,0) + getC(Z,k,0) )*x.at(i,j,k)
				);
	}

	inline constexpr Scalar innerStenc2D( const DomainFieldT& x, const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		return( 
				getC(X,i,-1)*x.at(i-1,j  ,k  ) + getC(X,i,1)*x.at(i+1,j  ,k  ) +
				getC(Y,j,-1)*x.at(i  ,j-1,k  ) + getC(Y,j,1)*x.at(i  ,j+1,k  ) +
				( getC(X,i,0) + getC(Y,j,0) )*x.at(i,j,k)
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
