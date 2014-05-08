#pragma once
#ifndef PIMPACT_DIVOP_HPP
#define PIMPACT_DIVOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

namespace Pimpact{


extern "C" {
  void OP_div( double* phiU, double* phiV, double* phiW, double* div);
}


/// \brief Divergence operator.
/// \ingroup BaseOperator
template<class Scalar,class Ordinal>
class Div {
public:
  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
//  typedef NonModeOp OpType;

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    for( int dir=0; dir<x.dim(); ++dir )
      x.exchange( dir, dir );
    OP_div(x.vec_[0],x.vec_[1],x.vec_[2],y.s_);
    y.changed();

  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }
};


/// \relates Div
template< class S, class O >
Teuchos::RCP< Div<S,O> > createDivOp() {
  return( Teuchos::rcp( new Div<S,O>() ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVOP_HPP
