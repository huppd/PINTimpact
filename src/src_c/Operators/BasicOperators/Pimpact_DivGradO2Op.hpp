#pragma once
#ifndef PIMPACT_DIVGRADO2OP_HPP
#define PIMPACT_DIVGRADO2OP_HPP


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

void OP_DivGradO2Op(
    const int& dimens,
    const int* const N,
    const int* const SR,
    const int* const ER,
    const int* const BL,
    const int* const BU,
    const int* const BCL,
    const int* const BCU,
    const double* const cdg1,
    const double* const cdg2,
    const double* const cdg3,
    const double* const phi,
          double* const Lap );

}



/// \brief "laplace" for pressure 2nd Order.
///
/// independent of \c StencilWidths
/// \ingroup BaseOperator
/// \todo instead of hardcode 2nd Order it would be pretty to use new space with \c StencilWidths<3,2>
/// \todo handle corner
template<class ST>
class DivGradO2Op {

public:

  using SpaceT = ST;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

  using TO = const Teuchos::Tuple<Ordinal,3>;

protected:

  using TS = const Teuchos::Tuple<Scalar*,3>;

  const Teuchos::RCP<const SpaceT> space_;

  TS c_;

	TO SR_;
	TO ER_;

public:


	DivGradO2Op( const Teuchos::RCP<const SpaceT>& space ):
		space_(space) {

			for( int i=0; i<3; ++i ) {
				// alocate stencilt
				Ordinal nTemp = 3*( space_->nLoc(i) - 1 + 1 );
				c_[i] = new Scalar[ nTemp ];

				// inner field bounds
				SR_[i] = 1;
				ER_[i] = space_->nLoc(i) - 1;

				if( space_->getBCLocal()->getBCL(i) >  0 ) SR_[i] = 2;
				if( space_->getBCLocal()->getBCL(i) == 0 ) SR_[i] = 1;
				if( space_->getBCLocal()->getBCU(i) >  0 ) ER_[i] = space_->nLoc(i) - 1;
				if( space_->getBCLocal()->getBCU(i) == 0 ) ER_[i] = space_->nLoc(i);
			}

			Op_getCDG(
					space_->dim(),
					space_->nGlo(),
					space_->nLoc(),
					space_->bl(),
					space_->bu(),
					space_->getBCLocal()->getBCL(),
					space_->getBCLocal()->getBCU(),
					space_->getCoordinatesGlobal()->getX( ECoord::X, EField::U ),
					space_->getCoordinatesGlobal()->getX( ECoord::Y, EField::V ),
					space_->getCoordinatesGlobal()->getX( ECoord::Z, EField::W ),
					space_->getCoordinatesLocal()->getX( ECoord::X, EField::S ),
					space_->getCoordinatesLocal()->getX( ECoord::Y, EField::S ),
					space_->getCoordinatesLocal()->getX( ECoord::Z, EField::S ),
					space_->getCoordinatesLocal()->getX( ECoord::X, EField::U ),
					space_->getCoordinatesLocal()->getX( ECoord::Y, EField::V ),
					space_->getCoordinatesLocal()->getX( ECoord::Z, EField::W ),
					c_[0],
					c_[1],
					c_[2] );
		}


	void apply(const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		x.exchange();
		OP_DivGradO2Op(
				space_->dim(),
				space_->nLoc(),
				getSR(),
				getER(),
				space_->bl(),
				space_->bu(),
				space_->getBCLocal()->getBCL(),
				space_->getBCLocal()->getBCU(),
				getC(X),
				getC(Y),
				getC(Z),
				x.getConstRawPtr(),
				y.getRawPtr() );

		//y.setCornersZero();

		y.changed();

	}

  void assignField ( const DomainFieldT& mv ) const {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(space_); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << " --- stencil: ---";
		out << " sr: " << SR_ << "\n";
		out << " er: " << ER_ << "\n";
    for( int dir=0; dir<3; ++dir ) {
      out << "\ndir: " << dir << "\n";
      for( int i=1; i<space_->nLoc(dir); ++i ) {
        out << "\ni: " << i << "\t(";
        for( int k=-1; k<=1; ++k ) {
					out << getC(dir,i,k) << "\t" ;
        }
        out << ")\n";
      }
      out << "\n";
    }
  }

	/// \name getters
	/// @{ 

	const Scalar* getC( const ECoord& dir) const  {
		return( getC( static_cast<const int&>(dir) ) );
  }

  const Scalar* getC( const int& dir) const  {
		return( c_[dir] );
  }

	const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const  {
		return( getC( static_cast<const int&>(dir), i, off ) );
  }

	const Scalar& getC( const int& dir, Ordinal i, Ordinal off ) const  {
		return( c_[dir][ off + 1 + (i-1)*3 ] );
  }

	const Ordinal* getSR() const { return( SR_.getRawPtr() ); }
	const Ordinal* getER() const { return( ER_.getRawPtr() ); }

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
