#pragma once
#ifndef PIMPACT_DIVOPGRAD_HPP
#define PIMPACT_DIVOPGRAD_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_ModeField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"

#include "Pimpact_LinearProblem.hpp"


namespace Pimpact {

/// \ingroup ModeOperator
template<class Scalar,class Ordinal>
class DivOpGrad {
public:
  typedef MultiField<ModeField<VectorField<Scalar,Ordinal> > > MVF;
  typedef typename MVF::Scalar S;
  typedef typename MVF::Ordinal O;
  typedef ModeField<ScalarField<S,O> > DomainFieldT;
  typedef ModeField<ScalarField<S,O> > RangeFieldT;
//  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
  typedef OperatorBase<MVF> OP;
  typedef ModeOp OpType;
private:
  Teuchos::RCP<MVF> temp0_;
  Teuchos::RCP<MVF> temp1_;
  Teuchos::RCP<Div<S,O>> div_;
  Teuchos::RCP<Grad<S,O>> grad_;
  Teuchos::RCP<LinearProblem<MVF> > op_;

public:
//  DivOpGrad():
//        temp0_(Teuchos::null), temp1_(Teuchos::null),
//        div_(Teuchos::rcp( new Div<S,O> ) ),
//        grad_(Teuchos::rcp( new Grad<S,O> ) ),
//        op_(Teuchos::null) {};
  DivOpGrad(
      const Teuchos::RCP<MVF>& temp=Teuchos::null,
      const Teuchos::RCP< LinearProblem<MVF> >& op=Teuchos::null ):
        temp0_(temp->clone(1)), temp1_(temp->clone(1)),
        div_(Teuchos::rcp( new Div<S,O> ) ),
        grad_(Teuchos::rcp( new Grad<S,O> ) ),
        op_(op) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    grad_->apply( x.getConstCField(), temp0_->getFieldPtr(0)->getCField() );
    grad_->apply( x.getConstSField(), temp0_->getFieldPtr(0)->getSField() );
    op_->solve( temp1_, temp0_);
    div_->apply( temp1_->getFieldPtr(0)->getConstCField(), y.getCField() );
    div_->apply( temp1_->getFieldPtr(0)->getConstSField(), y.getSField() );
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DivOpGrad



/// \relates DivOpGrad
template<class S, class O>
Teuchos::RCP<DivOpGrad<S,O> > createDivOpGrad(
    const Teuchos::RCP<MultiField<ModeField<VectorField<S,O> > > > & temp,
    const Teuchos::RCP< LinearProblem<MultiField<ModeField<VectorField<S,O> > > > > op ) {
  return( Teuchos::rcp( new DivOpGrad<S,O>( temp, op ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVOPGRAD_HPP
