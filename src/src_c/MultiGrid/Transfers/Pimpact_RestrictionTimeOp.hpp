/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_RESTRICTIONTIMEOP_HPP
#define PIMPACT_RESTRICTIONTIMEOP_HPP

#include "Teuchos_RCP.hpp"
#include <BelosTypes.hpp>
#include "Pimpact_Utils.hpp"
#include "Pimpact_TimeField.hpp"



namespace Pimpact {
template<class OperatorT>
class RestrictionTimeOp  {

public:

  using DomainFieldT = TimeField<typename OperatorT::DomainFieldT>;
  using RangeFieldT = TimeField<typename OperatorT::RangeFieldT>;

  using GridT = typename DomainFieldT::GridT;

  using FGridT = typename OperatorT::FGridT;
  using CGridT = typename OperatorT::CGridT;

  using Ordinal = typename GridT::Ordinal;

protected:

  Teuchos::RCP<OperatorT> op_; // restriction in grid

public:

  RestrictionTimeOp(
    const Teuchos::RCP<const GridT>& gridF,
    const Teuchos::RCP<const GridT>& gridC):
    op_(Teuchos::rcp(new OperatorT(gridF, gridC))) {};


  // x is fn in this case
  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    typename OperatorT::RangeFieldT temp(gridC());
    Ordinal d = gridF()->nLoc(3)/gridC()->nLoc(3);

    x.exchange();

    if (d > 1)
      op_->apply(x(0), temp);

    for(Ordinal i=gridF()->si(F::S, 3); i<gridF()->ei(F::S, 3); ++i)  {

      if ((i-gridF()->si(F::S, 3))%d==0) {

        Ordinal iC = (i-gridF()->si(F::S, 3))/d + gridC()->si(F::S, 3);
        op_->apply(x(i), y(iC));

        if (d > 1)
          y(iC).add(0.25, temp, 0.5, y(iC));
      } else if (d > 1) {
        op_->apply(x(i), temp);
        Ordinal iC = (i-gridF()->si(F::S, 3) - 1)/d + gridC()->si(F::S, 3);
        y(iC).add(1., y(iC), 0.25, temp);
      }
    }
    y.changed();
  }

  Teuchos::RCP<const GridT> gridC() const {
    return op_->gridC();
  };
  Teuchos::RCP<const GridT> gridF() const {
    return op_->gridF();
  };

  Teuchos::RCP<OperatorT> getOperatorPtr() {
    return op_;
  }

};

template<class OpT >
Teuchos::RCP<RestrictionTimeOp<OpT> > createRestrictionTimeOp(
    const Teuchos::RCP<const typename OpT::FGridT>& gridF, const Teuchos::RCP<const typename OpT::CGridT>& gridC) {

  return Teuchos::rcp(new RestrictionTimeOp<OpT>(gridF, gridC));
}

} // end of namespace Pimpact


#endif
