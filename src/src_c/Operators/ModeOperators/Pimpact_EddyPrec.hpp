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

  typedef Scalar S;
  typedef Ordinal O;

public:

  typedef ModeField< VectorField<S,O> > DomainFieldT;
  typedef ModeField< VectorField<S,O> > RangeFieldT;

  typedef MultiField<DomainFieldT> MVF;

  typedef OperatorBase<MVF> Op;

private:

  Teuchos::RCP<MVF> temp_;
  Teuchos::RCP<Op> op_;

public:

  EddyPrec(
      const Teuchos::RCP<MVF>& temp=Teuchos::null,
      const Teuchos::RCP<Op>& op=Teuchos::null ):
        temp_(temp),
        op_(op) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) {
    if( temp_.is_null() ) {
      temp_ = createMultiField( x.clone() );
    }

    temp_->getFieldPtr(0)->getCFieldPtr()->add( 1., x.getConstCField(),  1., x.getConstSField() );
    temp_->getFieldPtr(0)->getSFieldPtr()->add( 1., x.getConstCField(), -1., x.getConstSField() );

    op_->apply( *temp_, *createMultiField<DomainFieldT>(Teuchos::rcpFromRef(y)) );

    y.scale(0.5);

  }

  void assignField( const DomainFieldT& mv ) {
//    op_->assignField( createMultiField(Teuchos::rcpFromRef(mv))  );
  };

  bool hasApplyTranspose() const { return( false ); }

}; // end of class EddyPrec



/// \relates EddyPrec
template<class S, class O>
Teuchos::RCP<EddyPrec<S,O> > createEddyPrec(
    const Teuchos::RCP<typename EddyPrec<S,O>::MVF> & temp,
    const Teuchos::RCP<typename EddyPrec<S,O>::Op> op ) {
  return( Teuchos::rcp( new EddyPrec<S,O>( temp, op ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_EDDYPREC_HPP
