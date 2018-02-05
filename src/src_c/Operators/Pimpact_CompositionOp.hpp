#pragma once
#ifndef PIMPACT_COMPOSITIONOP_HPP
#define PIMPACT_COMPOSITIONOP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Grid.hpp"




namespace Pimpact {



/// \brief make the composition of two operators.
///
/// Both operators are applied sequentially.
/// the \c DomainFieldT \c OP2 has to equal to the \c RangeFieldT \c OP1.
/// \ingroup Operator
template<class OP1, class OP2 >
class CompositionOp {

public:

  using DomainFieldT = typename OP2::DomainFieldT;
  using RangeFieldT = typename OP1::RangeFieldT;

  using TempFieldT = typename OP2::RangeFieldT;

  using GridT = typename DomainFieldT::GridT;

protected:

  Teuchos::RCP<OP1> op1_;
  Teuchos::RCP<OP2> op2_;

public:

  CompositionOp(
    const Teuchos::RCP<OP1>& op1,
    const Teuchos::RCP<OP2>& op2):
    op1_(op1),
    op2_(op2) {};

  void apply(const DomainFieldT& x, RangeFieldT& y, const Belos::ETrans
      trans=Belos::NOTRANS) const {

    TempFieldT temp(grid());

    if(op1_.is_null()) {
//      op2_->apply(x, y);
    } else {
      op2_->apply(x, temp);
      op1_->apply(temp, y);
    }
  }

  void assignField(const DomainFieldT& field) {
    if(!op1_.is_null())
//      op1_->assignField(field);
      if(!op2_.is_null())
        op2_->assignField(field);
  };

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op1_->grid();
  };

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    if(!op1_.is_null()) op1_->setParameter(para);
    if(!op2_.is_null()) op2_->setParameter(para);
  }

  bool hasApplyTranspose() const {
    return op1_->hasApplyTranspose() && op2_->hasApplyTranspose() /*&& op3_->hasApplyTranspose()*/;
  }

  const std::string getLabel() const {
    return op1_->getLabel() + std::string("*") + op2_->getLabel();
  };

  void print(std::ostream& out=std::cout) const {
    out <<getLabel() <<":\n";
    op1_->print(out);
    op2_->print(out);
  }

}; // end of class CompositionOp



/// \relates CompositionOp
template<class OP1, class OP2>
Teuchos::RCP<CompositionOp<OP1, OP2> > createCompositionOp(
    const Teuchos::RCP<OP1>& op1=Teuchos::null,
    const Teuchos::RCP<OP2>& op2=Teuchos::null) {

  return Teuchos::rcp(new CompositionOp<OP1, OP2>(op1, op2));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOSITIONOP_HPP
