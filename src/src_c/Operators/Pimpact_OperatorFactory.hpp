#pragma once
#ifndef PIMPACT_OPERATORFACTORY_HPP
#define PIMPACT_OPERATORFACTORY_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Space.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_Operator.hpp"




namespace Pimpact {



/// \relates MultiOpWrap
/// \relates ModeOpWrap
template<class Op>
Teuchos::RCP< MultiOpWrap<ModeOpWrap<Op> > >
createMultiModeOpWrap( const Teuchos::RCP<Op>& op ) {

  return( createMultiOpWrap( createModeOpWrap( op ) ) );

}



/// \relates MultiOpWrap
/// \relates OperatorBase
/// \todo substitute tparam MF with MultiField<OP::DomainFIeld>...
template<class MF, class Op>
Teuchos::RCP< OperatorBase<MF> > createMultiOperatorBase( const Teuchos::RCP<Op>& op ) {

//  return(
//      createOperatorBase(
//          createMultiOpWrap<Op>( op )
//      )
//  );
  return(
      Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
          Teuchos::rcp( new OperatorPimpl< MultiOpWrap<Op> >( createMultiOpWrap<Op>(op) ) ) )
  );

}


/// \relates OperatorBase
/// \relates MultiOpWrap
/// \relates ModeOpWrap
template<class MF, class Op>
Teuchos::RCP<const OperatorBase<MF> > createMultiModeOperatorBase( const Teuchos::RCP<Op>& op ) {

  return(
      createOperatorBase(
          createMultiOpWrap(
              createModeOpWrap(op)
          )
      )
  );
//  return(
//      Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
//          Teuchos::rcp(
//              new OperatorPimpl< MultiOpWrap<ModeOpWrap<Op> > >(
//                  createMultiOpWrap(
//                      createModeOpWrap(op) )
//              ) ) ) );

}


/// \relates OperatorBase
/// \relates MultiOpWrap
/// \relates InverseOperator
template<class MF>
Teuchos::RCP< OperatorBase<MF> >
createInverseOperatorBase( const Teuchos::RCP< LinearProblem<MF> >& linProb ) {

  return(
      createOperatorBase(
          createInverseOperator<MF>(linProb)
      )
  );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORFACTORY_HPP
