/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_TRANSFERMODEOP_HPP
#define PIMPACT_TRANSFERMODEOP_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_ModeField.hpp"




namespace Pimpact {



template<class InterT>
class TransferModeOp {

public:

  using FGridT = typename InterT::FGridT;
  using CGridT = typename InterT::CGridT;

  using GridT = typename InterT::GridT;

  using DomainFieldT = ModeField<typename InterT::DomainFieldT>;
  using RangeFieldT = ModeField<typename InterT::RangeFieldT>;

protected:

  using OT = typename GridT::Ordinal;

  Teuchos::RCP<InterT> op_;

public:

  TransferModeOp(
    const Teuchos::RCP<const FGridT>& gridC,
    const Teuchos::RCP<const CGridT>& gridF):
    op_(Teuchos::rcp(new InterT(gridC, gridF))) {}

  TransferModeOp(
    const Teuchos::RCP<const FGridT>& gridC,
    const Teuchos::RCP<const CGridT>& gridF,
    const Teuchos::Tuple<int, GridT::dimension>& nb):
    op_(Teuchos::rcp(new InterT(gridC, gridF, nb))) {}



  template<class DT, class RT>
  void apply(const DT& x, RT& y) const {

    op_->apply(x.getCField(), y.getCField());
    op_->apply(x.getSField(), y.getSField());
  }


  void print(std::ostream& out=std::cout) const {

    out << "=== TransferModeOP ===\n";
    op_->print(out);
  }


}; // end of class TransferModeOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRANSFEROP_HPP
