#pragma once
#ifndef PIMPACT_GRADOP_HPP
#define PIMPACT_GRADOP_HPP

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"




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

void OP_SetBCZero(
    const int* const N,
    const int* const bl,
    const int* const bu,
    const int* const BCL,
    const int* const BCU,
    const int* const ss,
    const int* const nn,
    const double* phi );

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
            space_->getStencilWidths()->getDimNcbG(i),
            space_->getStencilWidths()->getNcbG(i),
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
    int dim = space_->dim();
    for( int i=0; i<dim; ++i) {
      x.exchange(i);

      OP_grad(
          i+1,
          space_->nLoc(),
          space_->bl(),
          space_->bu(),
          space_->gl(),
          space_->gu(),
          space_->sInd(i),
          space_->eInd(i),
          c_[i],
          x.s_,
          y.vec(i) );
      // necessary?
      OP_SetBCZero(
          space_->nLoc(),
          space_->bl(),
          space_->bu(),
          space_->getDomain()->getBCLocal()->getBCL(),
          space_->getDomain()->getBCLocal()->getBCU(),
          space_->sIndB(i),
          space_->eIndB(i),
          y.vec(i) );
      // necessary?
      // OP_bc_extrapolation( i+1, y.vec_[i] ); // doesnot work with Schurcomplement, not cleary what it does anyway
    }
    y.changed();
  }

  void assignField( const RangeFieldT& mv ) {};
  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }


  void print( std::ostream& out=std::cout ) const {
    out << " --- stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ndir: " << i << "\n";
      Ordinal nTemp1 = ( space_->nLoc(i) + 1 );
      Ordinal nTemp2 = ( space_->gu(i) - space_->gl(i) + 1 );
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

}; // end of class GradOp



/// \relates GradOp
template<class S=double, class O=int, int d=3>
Teuchos::RCP< GradOp<S,O,d> > createGradOp( const Teuchos::RCP< const Space<S,O,d> >& space ) {
  return( Teuchos::rcp( new GradOp<S,O,d>(space) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRADOP_HPP
