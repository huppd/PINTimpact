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

  using SpaceT = typename DomainFieldT::SpaceT;

  using FSpaceT = typename OperatorT::FSpaceT;
  using CSpaceT = typename OperatorT::CSpaceT;

  using Ordinal = typename SpaceT::Ordinal;

protected:

  Teuchos::RCP<OperatorT> op_; // restriction in space

public:

  RestrictionTimeOp(
    const Teuchos::RCP<const SpaceT>& spaceF,
    const Teuchos::RCP<const SpaceT>& spaceC):
    op_(Teuchos::rcp(new OperatorT(spaceF, spaceC))) {};


  // x is fn in this case
  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    typename OperatorT::RangeFieldT temp(spaceC());
    Ordinal d = spaceF()->nLoc(3)/spaceC()->nLoc(3);

    x.exchange();

    if (d > 1)
      op_->apply(x(0), temp);

    for(Ordinal i=spaceF()->si(F::S, 3); i<spaceF()->ei(F::S, 3); ++i)  {

      if ((i-spaceF()->si(F::S, 3))%d==0) {

        Ordinal iC = (i-spaceF()->si(F::S, 3))/d + spaceC()->si(F::S, 3);
        op_->apply(x(i), y(iC));

        if (d > 1)
          y(iC).add(0.25, temp, 0.5, y(iC));
      } else if (d > 1) {
        op_->apply(x(i), temp);
        Ordinal iC = (i-spaceF()->si(F::S, 3) - 1)/d + spaceC()->si(F::S, 3);
        y(iC).add(1., y(iC), 0.25, temp);
      }
    }
    y.changed();
  }

  Teuchos::RCP<const SpaceT> spaceC() const {
    return op_->spaceC();
  };
  Teuchos::RCP<const SpaceT> spaceF() const {
    return op_->spaceF();
  };

  Teuchos::RCP<OperatorT> getOperatorPtr() {
    return op_;
  }

};

template<class OpT >
Teuchos::RCP<RestrictionTimeOp<OpT> > createRestrictionTimeOp(
    const Teuchos::RCP<const typename OpT::FSpaceT>& spaceF, const Teuchos::RCP<const typename OpT::CSpaceT>& spaceC) {

  return Teuchos::rcp(new RestrictionTimeOp<OpT>(spaceF, spaceC));
}

} // end of namespace Pimpact


#endif
