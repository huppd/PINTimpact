#pragma once
#ifndef PIMPACT_GRADOP_HPP
#define PIMPACT_GRADOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

namespace Pimpact{


extern "C" {
  void OP_grad( const int& m, double* phi, double *grad );
  void OP_bc_extrapolation( const int& m, double* phi );
}


/// \ingroup BaseOperator
template<class Scalar, class Ordinal, int dimension=3>
class Grad {

public:

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal,dimension>  RangeFieldT;

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    int dim = x.dim();
    for( int i=0; i<dim; ++i) {
      x.exchange(i);
//      OP_grad( i+1, x.s_, y.vec_[i] );
      OP_grad( i+1, x.s_, y.vec(i) );
//      OP_bc_extrapolation( i+1, y.vec_[i] ); // doesnot work with Schurcomplement, not cleary what it does anyway
    }
    y.changed();
  }

  void assignField( const RangeFieldT& mv ) {};
  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class Grad



/// \relates Grad
template<class S=double, class O=int, int d=3>
Teuchos::RCP< Grad<S,O,d> > createGradOp() {
  return( Teuchos::rcp( new Grad<S,O,d>() ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRADOP_HPP
