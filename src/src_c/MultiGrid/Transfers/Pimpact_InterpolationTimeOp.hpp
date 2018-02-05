#pragma once
#ifndef PIMPACT_INTERPOLATIONTIMEOP_HPP
#define PIMPACT_INTERPOLATIONTIMEOP_HPP

#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Utils.hpp"

#include "Pimpact_TimeField.hpp"



namespace Pimpact {


/// \brief Operator wrapper for \c TimeField's.
/// \ingroup TimeOperator
/// wraps and \c OperatorT and adds the functionality of handling \c TimeField.
template<class OperatorT>
class InterpolationTimeOp  {

public:

  using DomainFieldT = TimeField<typename OperatorT::DomainFieldT>;
  using RangeFieldT = TimeField<typename OperatorT::RangeFieldT>;

  using GridT = typename DomainFieldT::GridT;

  using FGridT = GridT;
  using CGridT = GridT;

  using Ordinal = typename GridT::Ordinal;

protected:

  Teuchos::RCP<OperatorT> op_; //  interpolation operator in grid

public:

  InterpolationTimeOp(
    const Teuchos::RCP<const GridT>& gridC,
    const Teuchos::RCP<const GridT>& gridF):
    op_(Teuchos::rcp(new OperatorT(gridC, gridF))) {
  };

  /// \brief default apply
  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    Ordinal d = (gridF()->nLoc(3)) / (gridC()->nLoc(3));

    x.exchange();

    for(Ordinal i=gridC()->si(F::S, 3); i <= gridC()->ei(F::S, 3); ++i) {
      Ordinal iF = d*(i-gridC()->si(F::S, 3)) + gridF()->si(F::S, 3);
      op_->apply(x(i), y(iF));
      if (gridC()->nLoc(3)==1 && d>1)
        op_->apply(x(1), y(2));
    }

    if (d > 1 && gridC()->nLoc(3)>1) {

      for(int i=gridF()->si(F::S, 3) + 1; i <gridF()->ei(F::S, 3); i=i+d) {
        y(i).add(0.5, y(i-1), 0.5, y(i+1));
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


}; // end of class InterpolationTimeOp



/// \relates
template<class OpT >
Teuchos::RCP<InterpolationTimeOp<OpT> > createInterpolationTimeOp(
  const Teuchos::RCP<const typename OpT::CGridT>& gridC,
  const Teuchos::RCP<const typename OpT::FGridT>& gridF) {

  return Teuchos::rcp(new InterpolationTimeOp<OpT>(gridC, gridF));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTERPOLATIONTIMEOP_HPP
