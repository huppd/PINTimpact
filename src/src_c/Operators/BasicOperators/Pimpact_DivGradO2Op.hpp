#pragma once
#ifndef PIMPACT_DIVGRADO2OP_HPP
#define PIMPACT_DIVGRADO2OP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_ScalarField.hpp"




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

  template<class SpaceTT>
  friend class DivGradO2JSmoother;

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  const Teuchos::RCP<const SpaceT> space_;

  TO c_;

public:


  DivGradO2Op( const Teuchos::RCP<const SpaceT>& space ):
    space_(space) {

    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = 3*( space_->nLoc(i) - 1 + 1 );
      c_[i] = new Scalar[ nTemp ];
    }

        Op_getCDG(
            space_->dim(),
            space_->nGlo(),
            space_->nLoc(),
            space_->bl(),
            space_->bu(),
            space_->getDomain()->getBCLocal()->getBCL(),
            space_->getDomain()->getBCLocal()->getBCU(),
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
        space_->bl(),
        space_->bu(),
        space_->getDomain()->getBCLocal()->getBCL(),
        space_->getDomain()->getBCLocal()->getBCU(),
        c_[0],
        c_[1],
        c_[2],
        x.getConstRawPtr(),
        y.getRawPtr() );

	 SF_handle_corner(
			 space_->nLoc(),
			 space_->bl(),
			 space_->bu(),
			 space_->getDomain()->getBCLocal()->getBCL(),
			 space_->getDomain()->getBCLocal()->getBCU(),
			 y.getRawPtr() );

    y.changed();

  }

  void assignField ( const DomainFieldT& mv ) const {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(space_); };
  Teuchos::RCP<const SpaceT> getSpace() const { return( space_ ); }

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << " --- stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ndir: " << i << "\n";
      Ordinal nTemp1 = space_->nLoc(i) - 1 + 1;
      Ordinal nTemp2 = 3;
      for( int j=0; j<nTemp1; ++j ) {
        out << "\ni: " << j+1 << "\t(";
        for( int k=0; k<nTemp2; ++k ) {
          out << c_[i][k+nTemp2*j] <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }
  }


}; // end of class DivGradO2Op



/// \relates DivGradO2Op
template<class SpaceT>
Teuchos::RCP<const DivGradO2Op<SpaceT> > createDivGradO2Op(
    const Teuchos::RCP<const SpaceT>& space ) {
  return(
      Teuchos::rcp( new DivGradO2Op<SpaceT>(space) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVGRADO2OP_HPP
