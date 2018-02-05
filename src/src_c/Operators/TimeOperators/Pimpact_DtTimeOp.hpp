#pragma once
#ifndef PIMPACT_DTTIMEOP_HPP
#define PIMPACT_DTTIMEOP_HPP


#include <cmath>

#include "Pimpact_TimeField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {


/// \ingroup TimeOperator
template<class ST>
class DtTimeOp {

public:

  using GridT = ST;

  using DomainFieldT = TimeField<VectorField<GridT> >;
  using RangeFieldT = TimeField<VectorField<GridT> >;

protected:

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  Teuchos::RCP<const GridT> grid_;

  Scalar mulI_;

public:

  DtTimeOp(const Teuchos::RCP<const GridT>& grid):
    grid_(grid) {

    Scalar pi = 4.*std::atan(1.);
    Scalar idt = (static_cast<Scalar>(grid_->nGlo(3)))/2./pi;
    mulI_ =
      grid_->getDomainSize()->getAlpha2()*idt/grid_->getDomainSize()->getRe();
  };


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    x.exchange();

    for(Ordinal i=grid_->si(F::S, 3); i<=grid_->ei(F::S, 3); ++i) {
      y(i).add(mulI_, x(i), -mulI_, x(i-1));
    }

    y.changed();
  }

  void assignField(const DomainFieldT& mv) {};

  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return grid_;
  };

  const std::string getLabel() const {
    return "DtTimeOp ";
  };

}; // end of class DtTimeOp



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DTTIMEOP_HPP
