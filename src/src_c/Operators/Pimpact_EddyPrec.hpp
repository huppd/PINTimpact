#pragma once
#ifndef PIMPACT_EDDYPREC_HPP
#define PIMPACT_EDDYPREC_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"



namespace Pimpact {


/// \ingroup ModeOperator
template<class Scalar,class Ordinal>
class EddyPrec {

public:
  typedef MultiField<ModeField<VectorField<Scalar,Ordinal> > > MVF;
  typedef typename MVF::Scalar S;
  typedef typename MVF::Ordinal O;
  typedef VectorField<S,O>  DomainFieldT;
  typedef VectorField<S,O>  RangeFieldT;
//  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
  typedef OperatorBase<MVF> OP;
  typedef ModeOp OpType;

private:
  Teuchos::RCP<MVF> temp_;
  Teuchos::RCP<LinearProblem<S, MVF, OP > > op_;

public:
  EddyPrec():
        temp_(Teuchos::null),
        op_(Teuchos::null) {};
  EddyPrec(
      const Teuchos::RCP<MVF>& temp,
      const Teuchos::RCP< LinearProblem<S, MVF, OP > >& op ):
        temp_(temp->clone(1)),
        op_(op) {};

  void apply(const ModeField<DomainFieldT>& x, ModeField<RangeFieldT>& y) const {
    temp_->getFieldPtr(0)->getCFieldPtr()->add( 1., x.getConstCField(),  1., x.getConstSField() );
    temp_->getFieldPtr(0)->getSFieldPtr()->add( 1., x.getConstCField(), -1., x.getConstSField() );

    op_->solve( createMultiField<ModeField<DomainFieldT> >(Teuchos::rcpFromRef(y)), temp_ );
    y.scale(0.5);
  }

  bool hasApplyTranspose() const { return( false ); }

}; // end of class EddyPrec


template<class S, class O>
Teuchos::RCP<EddyPrec<S,O> > createEddyPrec(
    const Teuchos::RCP<MultiField<ModeField<VectorField<S,O> > > > & temp,
    const Teuchos::RCP< LinearProblem<S,MultiField<ModeField<VectorField<S,O> > >,
      OperatorBase<MultiField<ModeField<VectorField<S,O> > > > > > op ) {
  return( Teuchos::rcp( new EddyPrec<S,O>( temp, op ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_EDDYPREC_HPP
