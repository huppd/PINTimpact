#pragma once
#ifndef PIMPACT_OPERATORFACTORY_HPP
#define PIMPACT_OPERATORFACTORY_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_Grid.hpp"




namespace Pimpact {



/// \relates MultiOpWrap
/// \relates ModeOpWrap
template<class Op>
Teuchos::RCP<MultiOpWrap<ModeOpWrap<Op> > >
createMultiModeOpWrap(const Teuchos::RCP<Op>& op) {

  return createMultiOpWrap(createModeOpWrap(op));
}



/// \relates MultiOpWrap
/// \relates OperatorBase
template<class Op>
Teuchos::RCP<OperatorBase<MultiField<typename Op::DomainFieldT>, MultiField<typename Op::RangeFieldT> > >
createMultiOperatorBase(const Teuchos::RCP<Op>& op) {

  return Teuchos::rcp_dynamic_cast<OperatorBase<MultiField<typename
    Op::DomainFieldT>, MultiField<typename Op::RangeFieldT> > >(Teuchos::rcp(new
          OperatorPimpl<MultiOpWrap<Op> >(createMultiOpWrap<Op>(op))));
}


/// \relates OperatorBase
/// \relates MultiOpWrap
/// \relates ModeOpWrap
template<class MF, class Op>
Teuchos::RCP<const OperatorBase<MF> > createMultiModeOperatorBase(const Teuchos::RCP<Op>& op) {

  return createOperatorBase(createMultiOpWrap(createModeOpWrap(op)));
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORFACTORY_HPP
