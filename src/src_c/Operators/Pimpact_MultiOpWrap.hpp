/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_MULTIOPWRAP_HPP
#define PIMPACT_MULTIOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_MultiField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief Operator wrapper for \c MultiField's.
/// \ingroup MultiOperator
/// wraps and \c Operator and adds the functionality of handling \c MultiField.
template<class OT>
class MultiOpWrap  {

public:

  using OperatorT = OT;

  using GridT = typename OperatorT::GridT;

  using DomainFieldT = MultiField<typename OperatorT::DomainFieldT>;
  using RangeFieldT = MultiField<typename OperatorT::RangeFieldT>;

protected:

  Teuchos::RCP<OperatorT> op_;

public:

  MultiOpWrap(const Teuchos::RCP<const GridT>& grid):
    op_(Teuchos::rcp(new OperatorT(grid))) {}

  MultiOpWrap(const Teuchos::RCP<OperatorT>& op):op_(op) {}

//	template<class IOperatorT>
//	MultiOpWrap(const Teuchos::RCP<MultiOpWrap<IOperatorT> >& op):
//		op_(Teuchos::rcp(new OperatorT(op->getOperatorPtr()))) {}


  void apply(const DomainFieldT& x, RangeFieldT& y, const Belos::ETrans
      trans=Belos::NOTRANS) const {

    assert(x.getNumberVecs()==y.getNumberVecs());

    for(int i=0; i<x.getNumberVecs(); ++i)
      op_->apply(x.getField(i), y.getField(i));
  }


  void assignField(const DomainFieldT& mv) {

    op_->assignField(mv.getField(0));
  };


  bool hasApplyTranspose() const {
    return op_->hasApplyTranspose();
  }

  Teuchos::RCP<const OperatorT> getOperatorPtr() const {
    return op_;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op_->grid();
  };

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    op_->setParameter(para);
  }

  const std::string getLabel() const {
    return op_->getLabel();
  };

  void print(std::ostream& out=std::cout) const {
    op_->print(out);
  }

}; // end of class MultiOpWrap



/// \relates MultiOpWrap
/// \deprecated use create<MultiOpWrap>(op) instead
template<class OperatorT>
Teuchos::RCP<MultiOpWrap<OperatorT> > createMultiOpWrap(const Teuchos::RCP<OperatorT>& op) {

  return Teuchos::rcp(new MultiOpWrap<OperatorT>(op));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIOPWRAP_HPP
