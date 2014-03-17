#pragma once
#ifndef PIMPACT_DIVGRADOP_HPP
#define PIMPACT_DIVGRADOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"

namespace Pimpact{


extern "C" {
  void OP_div_grad( const bool& corner_yes, double* phi, double* lap );
}

/// \brief "laplace" for pressure.
/// \ingroup BaseOperator
/// \todo not workin properly?
/// \todo add temporary variable
/// \warning does not hold test.
template<class Scalar,class Ordinal>
class DivGradOp {
  Teuchos::RCP<VectorField<Scalar,Ordinal> > temp_;
  Teuchos::RCP<Div<Scalar,Ordinal> > div_;
  Teuchos::RCP<Grad<Scalar,Ordinal> > grad_;
public:
  typedef ScalarField<Scalar,Ordinal>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
  typedef NonModeOp OpType;

  DivGradOp():
    temp_(Teuchos::null),
    div_(Teuchos::rcp(new Div<Scalar,Ordinal>() )),
    grad_(Teuchos::rcp(new Grad<Scalar,Ordinal>() )) {};
  DivGradOp( const Teuchos::RCP<VectorField<Scalar,Ordinal> >& temp):
    temp_(temp),
    div_(Teuchos::rcp(new Div<Scalar,Ordinal>() )),
    grad_(Teuchos::rcp(new Grad<Scalar,Ordinal>() )) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    grad_->apply( x, *temp_ );
    div_->apply( *temp_, y );

//    OP_div_grad( true, x.s_, y.s_ );
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }
};

template<class S, class O>
Teuchos::RCP<DivGradOp<S,O> > createDivGradOp( const Teuchos::RCP<VectorField<S,O> >& temp ) {
  return(
      Teuchos::rcp( new DivGradOp<S,O>( temp ) )
  );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVGRADOP_HPP
