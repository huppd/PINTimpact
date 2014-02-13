#pragma once
#ifndef PIMPACT_OPERATORFACTORY_HPP
#define PIMPACT_OPERATORFACTORY_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_Operator.hpp"


namespace Pimpact {


template<class MV, class Op>
Teuchos::RCP<OperatorBase<MV> > createOperatorBase() {
  return(
      Teuchos::rcp_dynamic_cast<OperatorBase<MV> >(
          Teuchos::rcp( new OperatorPimpl<MV,Op>( createOperatorMV<Op>() ) )
          )
      );
}


template<class MV, class Op>
Teuchos::RCP<OperatorBase<MV> > createOperatorBase( const Teuchos::RCP<Op>& op ) {
  return(
      Teuchos::rcp_dynamic_cast<OperatorBase<MV> >(
          Teuchos::rcp( new OperatorPimpl<MV,Op>( createOperatorMV<Op>(op) ) ) )
      );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORFACTORY_HPP
