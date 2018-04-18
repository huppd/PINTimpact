/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_TIMEOPWRAP_HPP
#define PIMPACT_TIMEOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_TimeField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief Operator wrapper for \c TimeField's.
/// \ingroup TimeOperator
/// wraps and \c OperatorT and adds the functionality of handling \c TimeField.
template<class OperatorT, bool CNyes=false>
class TimeOpWrap  {

public:

  using DomainFieldT = TimeField<typename OperatorT::DomainFieldT>;
  using RangeFieldT = TimeField<typename OperatorT::RangeFieldT>;

  using GridT = typename DomainFieldT::GridT;

protected:

  using Ordinal = typename GridT::Ordinal;

  Teuchos::RCP<OperatorT> op_;

public:

  TimeOpWrap(const Teuchos::RCP<const GridT>& grid):
    op_(create<OperatorT>(grid)) {};

  TimeOpWrap(const Teuchos::RCP<OperatorT>& op): op_(op) {};

  /// \brief default apply
  void apply(const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N) const {

    if(true==CNyes) {

      typename OperatorT::RangeFieldT temp(grid());

      x.exchange();
      for(int i=0; i <=grid()->ei(F::S, 3); ++i) {

        op_->apply(x(i), temp);

        if(i>=grid()->si(F::S, 3))
          y(i).add(1., y(i), 0.5, temp);

        if(i+1<=grid()->ei(F::S, 3))
          y(i+1).add(0., y(i+1), 0.5, temp);
      }
    } else {
      for(Ordinal i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
        op_->apply(x(i) , y(i), add);
    }
    y.changed();
  }

  void assignField(const DomainFieldT& mv) {};

  bool hasApplyTranspose() const {
    return op_->hasApplyTranspose();
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op_->grid();
  };

  Teuchos::RCP<OperatorT> getOperatorPtr() {
    return op_;
  }
  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    op_->setParameter(para);
  }

  void print(std::ostream& out=std::cout) const {
    op_->print(out);
  }

  const std::string getLabel() const {
    return "TimeOpWrap (";
  };

}; // end of class TimeOpWrap



/// \relates TimeOpWrap
template<class OpT, bool CNY=false >
Teuchos::RCP<TimeOpWrap<OpT, CNY> > createTimeOpWrap(
  const Teuchos::RCP<OpT>& op) {

  return Teuchos::rcp(new TimeOpWrap<OpT, CNY>(op));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMEOPWRAP_HPP
