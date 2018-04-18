/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_INTRESCOMPOUNDOP_HPP
#define PIMPACT_INTRESCOMPOUNDOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Teuchos_TestForException.hpp"

#include "Pimpact_CompoundField.hpp"




namespace Pimpact {



template<class TransVT, class TransST>
class IntResCompoundOp {

public:

  using FGridT = typename TransVT::FGridT;
  using CGridT = typename TransVT::CGridT;

  using GridT = typename TransVT::GridT;

  using DomainFieldT = CompoundField<typename TransVT::DomainFieldT, typename TransST::DomainFieldT >;
  using RangeFieldT = CompoundField<typename TransVT::RangeFieldT, typename TransST::RangeFieldT>;

protected:

  Teuchos::RCP<TransVT> opV_;
  Teuchos::RCP<TransST> opS_;

public:

  IntResCompoundOp(
    const Teuchos::RCP<const FGridT>& gridC,
    const Teuchos::RCP<const CGridT>& gridF):
    opV_(Teuchos::rcp(new TransVT(gridC, gridF))),
    opS_(Teuchos::rcp(new TransST(gridC, gridF))) {}

  IntResCompoundOp(
    const Teuchos::RCP<const FGridT>& gridC,
    const Teuchos::RCP<const CGridT>& gridF,
    const Teuchos::Tuple<int, GridT::dimension>& nb):
    opV_(Teuchos::rcp(new TransVT(gridC, gridF))),
    opS_(Teuchos::rcp(new TransST(gridC, gridF))) {}


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    opV_->apply(x.getVField(), y.getVField());
    opS_->apply(x.getSField(), y.getSField());
  }


  void print(std::ostream& out=std::cout) const {

    out << "=== IntResCompoundOP ===\n";
    opV_->print(out);
    opS_->print(out);
  }

}; // end of class IntResCompoundOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTRESCOMPOUNDOP_HPP
