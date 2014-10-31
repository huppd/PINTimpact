#pragma once
#ifndef PIMPACT_INTERPOLATES2VDOP_HPP
#define PIMPACT_INTERPOLATES2VDOP_HPP

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"



namespace Pimpact{


extern "C" {

void OP_grad(
    const int& dir,
    const int* const N,
    const int* const bl,
    const int* const bu,
    const int* const gl,
    const int* const gu,
    const int* const ss,
    const int* const nn,
    const double* const c,
    const double* const phi,
    double* const grad );

}


/// \ingroup BaseOperator
template<class Scalar, class Ordinal, int dimension=3>
class InterpolateS2V {

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  Teuchos::RCP< const Space<Scalar,Ordinal,dimension> > space_;

  TO c_;

public:

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  RangeFieldT;

  InterpolateS2V( const Teuchos::RCP< const Space<Scalar,Ordinal,dimension> >& space):
    space_(space) {
    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->gu(i) - space_->gl(i) + 1);

      c_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->gl(i),
            space_->gu(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            2,
            i+1,
            0,
            0,
            true,
            space_->getStencilWidths()->getDimNcbG(i),
            space_->getStencilWidths()->getNcbG(i),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            space_->getCoordinatesLocal()->getX( i, i ),
            c_[i] );
    }

  };


  ~InterpolateS2V() {
    for( int i=0; i<3; ++i ) {
      delete[] c_[i];
    }
  }



  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    TEUCHOS_TEST_FOR_EXCEPTION(
        x.fType_ != S,
        std::logic_error,
        "Pimpact::InterpolateV2S:: can only interpolate from VectorField!!!\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
        y.fType_ == S,
        std::logic_error,
        "Pimpact::InterpolateV2S:: can only interpolate to Scalar!!!\n");

    int m = (int)y.fType_;

//    int dim = x.dim();
    x.exchange(m);

    OP_grad(
        m+1,
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->gl(),
        space_->gu(),
        space_->sInd(m),
        space_->eInd(m),
        c_[m],
        x.s_,
        y.s_ );
    y.changed();
  }

  void assignField( const RangeFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }


  void print( std::ostream& out=std::cout ) const {
    out << " --- stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ni: " << i << "\n( ";
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->gu(i) - space_->gl(i) + 1);
      for( int j=0; j<nTemp; ++j )
        out << c_[i][j] <<"\t";
      out << ")\n";
    }
  }

}; // end of class InterpolateS2V



/// \relates InterpolateS2V
template<class S=double, class O=int, int d=3>
Teuchos::RCP< InterpolateS2V<S,O,d> > createInterpolateS2V( const Teuchos::RCP< const Space<S,O,d> >& space ) {
  return( Teuchos::rcp( new InterpolateS2V<S,O,d>(space) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTERPOLATES2VDOP_HPP
