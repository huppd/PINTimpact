#pragma once
#ifndef PIMPACT_OPERATORFACTORY_HPP
#define PIMPACT_OPERATORFACTORY_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_Operator.hpp"


namespace Pimpact {


/// \relates MultiOpWrap
/// \relates ModeOpWrap
template<class Op>
Teuchos::RCP< MultiOpWrap<ModeOpWrap<Op> > > createMultiModeOpWrap( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null(op) )
    return( createMultiOpWrap( createModeOpWrap( Teuchos::rcp( new Op() ) ) ) );
  else
    return( createMultiOpWrap( createModeOpWrap( op ) ) );
}



/// \relates MultiOpWrap
/// \relates OperatorBase
template<class MF, class Op>
Teuchos::RCP< OperatorBase<MF> > createMultiOperatorBase( const Teuchos::RCP<Op>& op=Teuchos::null ) {
//  if( op.is_null() )
//    return(
//        Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
//            Teuchos::rcp(
//                new OperatorPimpl< MF, MultiOpWrap<Op> >(
//                    createMultiOpWrap<Op>( Teuchos::rcp( new Op() ) ) ) ) )
//            );
//  else
    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
            Teuchos::rcp( new OperatorPimpl< MF, MultiOpWrap<Op> >( createMultiOpWrap<Op>(op) ) ) )
            );
}


/// \relates OperatorBase
/// \relates MultiOpWrap
/// \relates ModeOpWrap
template<class MF, class Op>
Teuchos::RCP<OperatorBase<MF> > createMultiModeOperatorBase( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null(op) )
    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
            Teuchos::rcp(
                new OperatorPimpl< MF, MultiOpWrap<ModeOpWrap<Op> > >(
                    createMultiOpWrap(
                        createModeOpWrap(
                            Teuchos::rcp( new Op() ) ) )
                        ) ) ) );
  else
    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
            Teuchos::rcp(
                new OperatorPimpl< MF, MultiOpWrap<ModeOpWrap<Op> > >(
                    createMultiOpWrap(
                        createModeOpWrap(op) )
                        ) ) ) );
}


/// \relates OperatorBase
/// \relates MultiOpWrap
/// \relates InverseOperator
template<class MF>
Teuchos::RCP< OperatorBase<MF> > createInverseOperatorBase( const Teuchos::RCP< LinearProblem<MF> >& linProb ) {
  return( createOperatorBase<MF,InverseOperator<MF>>( createInverseOperator<MF>(linProb) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORFACTORY_HPP
