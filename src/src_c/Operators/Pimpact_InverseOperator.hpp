#pragma once
#ifndef PIMPACT_INVERSEOPERATOR_HPP
#define PIMPACT_INVERSEOPERATOR_HPP

#include "Pimpact_Types.hpp"

#include "Teuchos_RCP.hpp"
//#include "Pimpact_ScalarField.hpp"
//#include "Pimpact_VectorField.hpp"
//#include "Pimpact_MultiField.hpp"

//#include "Pimpact_HelmholtzOp.hpp"
//#include "Pimpact_OperatorMV.hpp"

#include "Pimpact_LinearProblem.hpp"




namespace Pimpact{


/// \ingroup Operator
/// \tparam MultiField has to be of type \c Pimpact::MultiField
template< class MultiField >
class InverseOperator {

  typedef MultiField           MF;
  typedef typename MF::Scalar  S;
  typedef OperatorBase<MF>     Op;

public:

  typedef MultiField  DomainFieldT;
  typedef MultiField  RangeFieldT;

protected:

  Teuchos::RCP< LinearProblem<MF> > linprob_;

public:

  InverseOperator( const Teuchos::RCP< LinearProblem<MF> >& linprob=Teuchos::null ):
    linprob_(linprob) {};

  void apply(const MF& x, MF& y, Belos::ETrans trans=Belos::NOTRANS ) const {
    linprob_->solve( Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x) );
  }

  void assignField( const DomainFieldT& mv ) {
    auto prob = linprob_->getProblem();
    Teuchos::rcp_const_cast<Op>( prob->getOperator() )->assignField( mv );
    if( prob->isLeftPrec() ) {
      auto opPrec = Teuchos::rcp_const_cast<Op>( prob->getLeftPrec() );
      opPrec->assignField( mv );
    }
  };

  bool hasApplyTranspose() const { return( false ); }

}; // end of class InverseOperator



/// \relates InverseOperator
template< class MF>
Teuchos::RCP< InverseOperator<MF> >
createInverseOperator(
    const Teuchos::RCP<LinearProblem<MF> >& linprob ) {

    return( Teuchos::rcp( new InverseOperator<MF>( linprob ) ) );

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_INVERSEOPERATOR_HPP
