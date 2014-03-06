#pragma once
#ifndef PIMPACT_DIVOPGRAD_HPP
#define PIMPACT_DIVOPGRAD_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_ScalarField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"


namespace Pimpact {


template<class Scalar,class Ordinal>
class DivOpGrad {
public:
  typedef MultiField<ModeField<VectorField<Scalar,Ordinal> > > MVF;
  typedef typename MVF::Scalar S;
  typedef typename MVF::Ordinal O;
  typedef ScalarField<S,O>  DomainFieldT;
  typedef ScalarField<S,O>  RangeFieldT;
//  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
  typedef OperatorBase<MVF> OP;
  typedef ModeOp OpType;
private:
  Teuchos::RCP<MVF> temp0_;
  Teuchos::RCP<MVF> temp1_;
  Teuchos::RCP<Div<S,O>> div_;
  Teuchos::RCP<Grad<S,O>> grad_;
  Teuchos::RCP<LinearProblem<S, MVF, OP > > op_;

public:
  DivOpGrad():
        temp0_(Teuchos::null), temp1_(Teuchos::null),
        div_(Teuchos::rcp( new Div<S,O> ) ),
        grad_(Teuchos::rcp( new Grad<S,O> ) ),
        op_(Teuchos::null) {};
  DivOpGrad(
      const Teuchos::RCP<MVF>& temp,
      const Teuchos::RCP< LinearProblem<S, MVF, OP > >& op ):
        temp0_(temp->clone(1)), temp1_(temp->clone(1)),
        div_(Teuchos::rcp( new Div<S,O> ) ),
        grad_(Teuchos::rcp( new Grad<S,O> ) ),
        op_(op) {};

  void apply(const ModeField<DomainFieldT>& x, ModeField<RangeFieldT>& y) const {
    grad_->apply( *x.getConstFieldC(), *temp0_->getField(0).getFieldC() );
    grad_->apply( *x.getConstFieldS(), *temp0_->getField(0).getFieldS() );
    op_->solve( temp1_, temp0_);
    div_->apply( *temp1_->getField(0).getConstFieldC(), *y.getFieldC() );
    div_->apply( *temp1_->getField(0).getConstFieldS(), *y.getFieldS() );
  }

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DivOpGrad


template<class S, class O>
Teuchos::RCP<DivOpGrad<S,O> > createDivOpGrad(
    const Teuchos::RCP<MultiField<ModeField<VectorField<S,O> > > > & temp,
    const Teuchos::RCP< LinearProblem<S,MultiField<ModeField<VectorField<S,O> > >,
      OperatorBase<MultiField<ModeField<VectorField<S,O> > > > > > op ) {
  return( Teuchos::rcp( new DivOpGrad<S,O>( temp, op ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVOPGRAD_HPP
