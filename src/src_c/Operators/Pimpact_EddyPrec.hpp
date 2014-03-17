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
  typedef ModeField<VectorField<S,O> > DomainFieldT;
  typedef ModeField<VectorField<S,O> > RangeFieldT;
//  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
//  typedef OperatorBase<MVF> OP;
  typedef ModeOp OpType;

private:
  Teuchos::RCP<MVF> temp_;
  Teuchos::RCP<LinearProblem<MVF> > op_;

public:
  EddyPrec(
      const Teuchos::RCP<MVF>& temp=Teuchos::null,
      const Teuchos::RCP< LinearProblem<MVF> >& op=Teuchos::null ):
        temp_(temp->clone(1)),
        op_(op) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    temp_->getFieldPtr(0)->getCFieldPtr()->add( 1., x.getConstCField(),  1., x.getConstSField() );
    temp_->getFieldPtr(0)->getSFieldPtr()->add( 1., x.getConstCField(), -1., x.getConstSField() );

    op_->solve( createMultiField<DomainFieldT>(Teuchos::rcpFromRef(y)), temp_ );
    y.scale(0.5);
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class EddyPrec


template<class S, class O>
Teuchos::RCP<EddyPrec<S,O> > createEddyPrec(
    const Teuchos::RCP<MultiField<ModeField<VectorField<S,O> > > > & temp,
    const Teuchos::RCP< LinearProblem<MultiField<ModeField<VectorField<S,O> > > > > op ) {
  return( Teuchos::rcp( new EddyPrec<S,O>( temp, op ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_EDDYPREC_HPP
