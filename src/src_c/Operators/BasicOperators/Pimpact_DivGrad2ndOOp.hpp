#pragma once
#ifndef PIMPACT_DIVGRAD2NDOP_HPP
#define PIMPACT_DIVGRAD2NDOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_ScalarField.hpp"




namespace Pimpact{


extern "C" {

//  void SF_level( double* const phi );
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

void OP_DivGrad2ndOOp(
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
/// \ingroup BaseOperator
/// \todo instead of hardcode 2nd Order it would be pretty to use new space with StencilWidth<3,2>
/// \todo handle corner
template<class Scalar,class Ordinal, int dimension=3>
class DivGrad2ndOOp {

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  const Teuchos::RCP<const Space<Scalar,Ordinal,dimension> > space_;

  TO c_;

public:

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  RangeFieldT;


  DivGrad2ndOOp( const Teuchos::RCP<const Space<Scalar,Ordinal,dimension> >& space ):
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
    OP_DivGrad2ndOOp(
        space_->dim(),
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->getDomain()->getBCLocal()->getBCL(),
        space_->getDomain()->getBCLocal()->getBCU(),
        c_[0],
        c_[1],
        c_[2],
        x.s_,
        y.s_ );
    //    SF_level( y.s_ );
    // handlecorners
    y.changed();

  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

  void print( std::ostream& out=std::cout ) const {
    out << " --- stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ndir: " << i << "\n";
      Ordinal nTemp1 = space_->nLoc(i) - 1 + 1;
      Ordinal nTemp2 = 3;
      for( int j=0; j<nTemp1; ++j ) {
        out << "\ni: " << j << "\t(";
        for( int k=0; k<nTemp2; ++k ) {
          out << c_[i][k+nTemp2*j] <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }
  }

}; // end of class DivGrad2ndOOp



/// \relates DivGrad2ndOOp
template<class S, class O, int d=3>
Teuchos::RCP< DivGrad2ndOOp<S,O,d> > createDivGrad2ndOOp(
    const Teuchos::RCP<const Space<S,O,d> >& space ) {
  return(
      Teuchos::rcp( new DivGrad2ndOOp<S,O,d>(space) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVGRAD2NDOP_HPP
