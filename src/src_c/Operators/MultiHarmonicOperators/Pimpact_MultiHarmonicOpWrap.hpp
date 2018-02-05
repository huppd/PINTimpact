#pragma once
#ifndef PIMPACT_MULTIHARMONICOPWRAP_HPP
#define PIMPACT_MULTIHARMONICOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

//#include "Pimpact_MultiField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_Utils.hpp"



namespace Pimpact {


/// \brief Operator wrapper.
/// \ingroup MultiHarmonicOperator
/// wraps a \ref BaseOperator "base operator" and adds the functionality of handling \c MultiHarmonicField's.
template<class OpT>
class MultiHarmonicOpWrap  {


public:

  using DomainFieldT = MultiHarmonicField<typename OpT::DomainFieldT>;
  using RangeFieldT  = MultiHarmonicField<typename OpT::RangeFieldT>;

  using GridT = typename DomainFieldT::GridT;

protected:

  Teuchos::RCP<OpT> op_;

public:

  MultiHarmonicOpWrap(const Teuchos::RCP<const GridT>& grid):
    op_(Teuchos::rcp(new OpT(grid))) {};

  MultiHarmonicOpWrap(const Teuchos::RCP<OpT>& op): op_(op) {};


  void apply(const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N) const {

    if(0==grid()->si(F::U, 3))
      op_->apply(x.get0Field(), y.get0Field(), add);

    for(typename GridT::Ordinal i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i) {
      op_->apply(x.getCField(i), y.getCField(i), add);
      op_->apply(x.getSField(i), y.getSField(i), add);
    }

    y.changed();
  };


  void assignField(const DomainFieldT& mv) {
    op_->assignField(mv.get0Field());
  };

  bool hasApplyTranspose() const {
    return op_->hasApplyTranspose();
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op_->grid();
  };

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    op_->setParameter(para);
  }

  Teuchos::RCP<OpT> getOperatorPtr() {
    return op_;
  }

  const std::string getLabel() const {
    return "MH_"+op_->getLabel();
  };

  void print(std::ostream& out=std::cout) const {
    out << getLabel() <<":\n";
    op_->print(out);
  }

}; // end of class MultiHarmonicOpWrap



/// \relates MultiHarmonicOpWrap
template<class OpT>
Teuchos::RCP<MultiHarmonicOpWrap<OpT> >
createMultiHarmonicOpWrap(const Teuchos::RCP<OpT>& op) {
  return Teuchos::rcp(new MultiHarmonicOpWrap<OpT>(op));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIHARMONICOPWRAP_HPP
