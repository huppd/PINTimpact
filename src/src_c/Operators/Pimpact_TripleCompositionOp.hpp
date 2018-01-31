#pragma once
#ifndef PIMPACT_TRIPLECOMPOSITIONOP_HPP
#define PIMPACT_TRIPLECOMPOSITIONOP_HPP


#include "Teuchos_RCP.hpp"

#include "BelosTypes.hpp"



namespace Pimpact {


/// \brief make the composition of two operators.
///
/// \f[ op1 op2 op3 \mathbf{x} = \mathbf{y} \f]
/// Both operators are applied sequentially.
/// the \c DomainFieldT \c OP1 has to equal to the \c RangeFieldT \c OP2.
/// \ingroup Operator
template<class OP1, class OP2, class OP3=OP2 >
class TripleCompositionOp {

public:

  using DomainFieldT = typename OP3::DomainFieldT;
  using RangeFieldT = typename OP1::RangeFieldT;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  Teuchos::RCP<OP1> op1_;
  Teuchos::RCP<OP2> op2_;
  Teuchos::RCP<OP3> op3_;

public:

  TripleCompositionOp(
    const Teuchos::RCP<OP1>& op1,
    const Teuchos::RCP<OP2>& op2,
    const Teuchos::RCP<OP3>& op3):
    op1_(op1),
    op2_(op2),
    op3_(op3) {};

  void apply(const DomainFieldT& x, RangeFieldT& y, const Belos::ETrans
      trans=Belos::NOTRANS) const {

    typename OP2::RangeFieldT temp1(space()); // has to be equal to OP2::DomainFieldT

    op3_->apply(x, temp1);

    typename OP2::DomainFieldT temp2(space()); // has to be equal to OP3::DomainFieldT

    //temp1.extrapolateBC();

    op2_->apply(temp1, temp2);

    temp2.extrapolateBC();

    op1_->apply(temp2, y);
  }

  /// \note here nothing happens, because it is assumed to be done somewhere else
  void assignField(const RangeFieldT& field) { };

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op1_->space();
  };

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    if(!op1_.is_null()) op1_->setParameter(para);
    if(!op2_.is_null()) op2_->setParameter(para);
    if(!op3_.is_null()) op3_->setParameter(para);
  }

  bool hasApplyTranspose() const {
    return op1_->hasApplyTranspose() && op2_->hasApplyTranspose() && op3_->hasApplyTranspose();
  }

  const std::string getLabel() const {
    return std::string("(") + op1_->getLabel() + std::string(" * ") + op2_->getLabel() + std::string(" * ") + op3_->getLabel() + std::string(")");
  };


  void print(std::ostream& out=std::cout) const {
    out <<"TripleComposition: " <<getLabel() <<"\n";
    op1_->print(out);
    op2_->print(out);
    op3_->print(out);
  }

}; // end of class TripleCompositionOp



/// \relates TripleCompositionOp
template<class OP1, class OP2, class OP3>
Teuchos::RCP<TripleCompositionOp<OP1, OP2, OP3> > createTripleCompositionOp(
  const Teuchos::RCP<OP1>& op1,
  const Teuchos::RCP<OP2>& op2,
  const Teuchos::RCP<OP3>& op3) {

  return Teuchos::rcp(new TripleCompositionOp<OP1, OP2, OP3>(op1, op2, op3));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRIPLECOMPOSITIONOP_HPP
