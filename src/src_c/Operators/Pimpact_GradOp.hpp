#pragma once
#ifndef PIMPACT_GRADOP_HPP
#define PIMPACT_GRADOP_HPP

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{


extern "C" {
//  void OP_grad( const int& m, double* phi, double *grad );
  void OP_grad(
      const int& dir,
//      const int& dimens,
      const int* const N,
      const int* const bl,
      const int* const bu,
      const int* const gl,
      const int* const gu,
      const int* const BC_L,
      const int* const BC_U,
      const int* const ss,
      const int* const nn,
      const int* const sb,
      const int* const nb,
      const double* const c,
      const double* const phi,
      double* const grad );

  void OP_bc_extrapolation( const int& m, double* phi );
}


/// \ingroup BaseOperator
template<class Scalar, class Ordinal, int dimension=3>
class GradOp {

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  Teuchos::RCP< const Space<Scalar,Ordinal,dimension> > space_;

  TO c_;

public:

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal,dimension>  RangeFieldT;

  GradOp( const Teuchos::RCP< const Space<Scalar,Ordinal,dimension> >& space):
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
            1,
            0,
            true,
            space_->getFieldSpace()->getDimNcbG(i),
            space_->getFieldSpace()->getNcbG(i),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            space_->getCoordinatesLocal()->getX( i, i ),
            c_[i] );
    }

  };


  ~GradOp() {
    for( int i=0; i<3; ++i ) {
      delete[] c_[i];
    }
  }



  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    int dim = x.dim();
    for( int i=0; i<dim; ++i) {
      x.exchange(i);
//      OP_grad( i+1, x.s_, y.vec_[i] );
//      OP_grad( i+1, x.s_, y.vec(i) );
//        c_[0],
//        c_[1],
//        c_[2],
//        x.vecC(0),
//        x.vecC(1),
//        x.vecC(2),
//        y.s_ );
  OP_grad(
      i+1,
//      space_->dim(),
      space_->nLoc(),
      space_->bl(),
      space_->bu(),
      space_->gl(),
      space_->gu(),
      space_->getDomain()->getBCLocal()->getBCL(),
      space_->getDomain()->getBCLocal()->getBCU(),
      y.getField(i).sInd(),
      y.getField(i).eInd(),
      y.getField(i).sIndB(),
      y.getField(i).eIndB(),
      c_[i],
      x.s_,
      y.vec(i) );
//      OP_bc_extrapolation( i+1, y.vec_[i] ); // doesnot work with Schurcomplement, not cleary what it does anyway
    }
    y.changed();
  }

  void assignField( const RangeFieldT& mv ) {};
  void assignField( const DomainFieldT& mv ) {};

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

}; // end of class GradOp



/// \relates GradOp
template<class S=double, class O=int, int d=3>
Teuchos::RCP< GradOp<S,O,d> > createGradOp( const Teuchos::RCP< const Space<S,O,d> >& space ) {
  return( Teuchos::rcp( new GradOp<S,O,d>(space) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRADOP_HPP
