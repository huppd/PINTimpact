#pragma once
#ifndef PIMPACT_OPERATORFACTORY_HPP
#define PIMPACT_OPERATORFACTORY_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_Operator.hpp"


namespace Pimpact {

template<class Op>
Teuchos::RCP< MultiOpWrap<ModeOpWrap<Op> > > createMultiModeOpWrap( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null(op) )
    return( createMultiOpWrap( createModeOpWrap( Teuchos::rcp( new Op() ) ) ) );
  else
    return( createMultiOpWrap( createModeOpWrap( op ) ) );
}

/// \deprecated
//template<class MV, class Op>
//Teuchos::RCP<OperatorBase<MV> > createOperatorBaseMV( const Teuchos::RCP<Op>& op ) {
//  return(
//      Teuchos::rcp_dynamic_cast< OperatorBase<MV> >(
//          Teuchos::rcp( new OperatorPimpl< MV, MultiOpWrap<Op> >( createMultiOpWrap(op) ) ) )
//      );
//}


template<class MV, class Op>
Teuchos::RCP<OperatorBase<MV> > createMultiOperatorBase( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null(op) )
    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<MV> >(
            Teuchos::rcp(
                new OperatorPimpl< MV, MultiOpWrap<Op> >(
                    createMultiOpWrap<Op>( Teuchos::rcp( new Op() ) ) ) ) )
            );
  else
    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<MV> >(
            Teuchos::rcp( new OperatorPimpl< MV, MultiOpWrap<Op> >( createMultiOpWrap<Op>(op) ) ) )
            );
}


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



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORFACTORY_HPP
