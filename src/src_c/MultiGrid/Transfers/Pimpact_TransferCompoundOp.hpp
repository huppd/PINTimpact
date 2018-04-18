/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_TRANSFERCOMPOUNDOP_HPP
#define PIMPACT_TRANSFERCOMPOUNDOP_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_CompoundField.hpp"




namespace Pimpact {



template<class TransVT, class TransST>
class TransferCompoundOp {

public:

  using FGridT = typename TransVT::FGridT;
  using CGridT = typename TransVT::CGridT;

  using GridT = typename TransVT::GridT;

  using DomainFieldT = CompoundField<typename TransVT::DomainFieldT, typename TransST::DomainFieldT >;
  using RangeFieldT = CompoundField<typename TransVT::RangeFieldT,  typename TransST::RangeFieldT  >;

protected:

  Teuchos::RCP<TransVT> opV_;
  Teuchos::RCP<TransST> opS_;

public:

  TransferCompoundOp(
    const Teuchos::RCP<const FGridT>& gridC,
    const Teuchos::RCP<const CGridT>& gridF):
    opV_(Teuchos::rcp(new TransVT(gridC, gridF))),
    opS_(Teuchos::rcp(new TransST(gridC, gridF))) {}

  TransferCompoundOp(
    const Teuchos::RCP<const FGridT>& gridC,
    const Teuchos::RCP<const CGridT>& gridF,
    const Teuchos::Tuple<int, GridT::dimension>& nb):
    opV_(Teuchos::rcp(new TransVT(gridC, gridF))),
    opS_(Teuchos::rcp(new TransST(gridC, gridF))) {}


  template<class DT, class RT>
  void apply(const DT& x, RT& y) const {

    opV_->apply(x.getVField(), y.getVField());
    opS_->apply(x.getSField(), y.getSField());
  }

  void print(std::ostream& out=std::cout) const {

    out << "=== TransferCompoundOP ===\n";
    opV_->print(out);
    opS_->print(out);
  }

}; // end of class TransferCompoundOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRANSFERCOMPOUNDOP_HPP
