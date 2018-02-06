#pragma once
#ifndef PIMPACT_TRANSFERTIMEOP_HPP
#define PIMPACT_TRANSFERTIMEOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Teuchos_TestForException.hpp"



namespace Pimpact {


template<class TransVT>
class TransferTimeOp {

public:

  using FGridT = typename TransVT::FGridT;
  using CGridT = typename TransVT::CGridT;

  using GridT = typename TransVT::GridT;

  using DomainFieldT = TimeField<typename TransVT::DomainFieldT>;
  using RangeFieldT = TimeField<typename TransVT::RangeFieldT>;

protected:

  Teuchos::RCP<TransVT> opV_;

public:

  TransferTimeOp(
    const Teuchos::RCP<const FGridT>& gridC,
    const Teuchos::RCP<const CGridT>& gridF):
    opV_(Teuchos::rcp(new TransVT(gridC, gridF))) {}


  template<class DT, class RT>
  void apply(const DT& x, RT& y) const {

    for (int i=x.grid()->si(F::S, 3); i<x.grid()->ei(F::S, 3); ++i)
      opV_->apply(x(i), y(i));

    y.changed();
  }


  void print(std::ostream& out=std::cout) const {
    out << "=== TransferTimeOP ===\n";
    opV_->print(out);
  }

}; // end of class TransferTime


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRANSFERCOMPOUNDOP_HPP
