#pragma once
#ifndef PIMPACT_PRECINVERSEOP_HPP
#define PIMPACT_PRECINVERSEOP_HPP


#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief hides all the linear solver typeerasiure stuff
///
/// \todo inherit from InverseOp
/// \ingroup Operator
template<class OpT, template<class> class PT >
class PrecInverseOp {

public:

  using OperatorT = OpT;
  using PreconditionerT = PT<OpT>;

  using SpaceT = typename OperatorT::SpaceT;

  using DomainFieldT = typename OperatorT::DomainFieldT;
  using RangeFieldT = typename OperatorT::RangeFieldT;

  using MF = typename Pimpact::MultiField<DomainFieldT>;

  using MOpT = OperatorBase<MF>;

protected:

  Teuchos::RCP<LinearProblem<MF> > linprob_;

public:

  PrecInverseOp(const Teuchos::RCP<OperatorT>& op) {

    auto para =
      createLinSolverParameter("GMRES", 1.e-12, -1, Teuchos::rcp(&std::cout, false), 200);

    linprob_ = createLinearProblem<MF>(
                 createMultiOperatorBase(op),
                 create<MF>(op->space()),
                 create<MF>(op->space()),
                 para,
                 "GMRES");

    linprob_->setRightPrec(
      createMultiOperatorBase(
        create<PreconditionerT>(op)));
  }


  PrecInverseOp(const Teuchos::RCP<OperatorT>& op,
                 const Teuchos::RCP<Teuchos::ParameterList>& pl):
    linprob_(createLinearProblem<MF>(
                createMultiOperatorBase(op),
                create<MF>(op->space()),
                create<MF>(op->space()),
                Teuchos::sublist(pl, "Solver"),
                pl->get<std::string>("Solver name", "GMRES"))) {

    if(pl->get<bool>("LeftPrec", false))
      linprob_->setLeftPrec(
        createMultiOperatorBase(
          create<PreconditionerT>(op, Teuchos::sublist(pl, "Preconditioner"))));
    else
      linprob_->setRightPrec(
        createMultiOperatorBase(
          create<PreconditionerT>(op, Teuchos::sublist(pl, "Preconditioner"))));
  }


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    linprob_->solve(
      wrapMultiField(Teuchos::rcpFromRef(y)),
      wrapMultiField(Teuchos::rcpFromRef(const_cast<DomainFieldT&>(x))));
  }


  void assignField(const DomainFieldT& mvr) {

    auto mv = wrapMultiField(
                Teuchos::rcpFromRef<DomainFieldT>(const_cast<DomainFieldT&>(mvr)));

    auto prob = linprob_->getProblem();

    Teuchos::rcp_const_cast<MOpT>(prob->getOperator())->assignField(*mv);

    if(prob->isLeftPrec()) {
      auto opPrec = Teuchos::rcp_const_cast<MOpT>(prob->getLeftPrec());
      opPrec->assignField(*mv);
    }

    if(prob->isRightPrec()) {
      auto opPrec = Teuchos::rcp_const_cast<MOpT>(prob->getRightPrec());
      opPrec->assignField(*mv);
    }
  };


  Teuchos::RCP<const SpaceT> space() const {
    return linprob_->space();
  };

  Teuchos::RCP<LinearProblem<MF> > getLinearProblem() {
    return linprob_;
  }


  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    auto prob = linprob_->getProblem();

    Teuchos::rcp_const_cast<MOpT>(prob->getOperator())->setParameter(para);

    if(prob->isLeftPrec()) {
      auto opPrec = Teuchos::rcp_const_cast<MOpT>(prob->getLeftPrec());
      opPrec->setParameter(para);
    }

    if(prob->isRightPrec()) {
      auto opPrec = Teuchos::rcp_const_cast<MOpT>(prob->getRightPrec());
      opPrec->setParameter(para);
    }
  }


  /// \brief Set left preconditioner (\c LP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setLeftPrec(const Teuchos::RCP<const OperatorBase<MF> > &LP) {
    linprob_->getProblem()->setLeftPrec(LP);
  }

  /// \brief Set right preconditioner (\c RP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setRightPrec(const Teuchos::RCP<const OperatorBase<MF> > &RP) {
    linprob_->getProblem()->setRightPrec(RP);
  }

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return linprob_->getProblem()->getOperator()->getLabel() + std::string("^-1(prec) ");
  };

  void print(std::ostream& out=std::cout) const {
    out <<"PrecInverse:\n";
    linprob_->print(out);
  }

}; // end of class PrecInverseOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_PRECINVERSEOP_HPP
