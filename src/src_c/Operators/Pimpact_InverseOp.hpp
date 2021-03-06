/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_INVERSEOP_HPP
#define PIMPACT_INVERSEOP_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosSolverManager.hpp"

#include "Pimpact_EmptyProjector.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_MultiOpWrap.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// hides all the linear solver typeerasiure stuff
///
/// \todo merge with Prec...
/// \ingroup Operator
template<class OpT, template<class> class  ProjectorT=EmptyProjector >
class InverseOp {

public:

  using OperatorT = OpT;
  using ProT = ProjectorT<OpT>;

  using GridT = typename OperatorT::GridT;

  using DomainFieldT = typename OperatorT::DomainFieldT;
  using RangeFieldT = typename OperatorT::RangeFieldT;

  using MF = typename Pimpact::MultiField<DomainFieldT>;

  using MOpT = OperatorBase<MF>;

protected:

  using ST = typename GridT::Scalar;

  bool level_;
  bool nullspaceOrtho_;
  bool initZero_;
  bool debug_;

  ST relTol_;

  std::string solverName_;

  Teuchos::RCP<Teuchos::ParameterList > solverParameter_;
  Teuchos::RCP<Belos::LinearProblem<ST, MF, MOpT> > problem_;

  Teuchos::RCP<const ProjectorT<OpT> > projector_;

public:

  /// should be avoided(used in peri_navierSMGXML)
  /// does not work
  InverseOp(const Teuchos::RCP<const GridT>& grid):
    level_(false),
    nullspaceOrtho_(false),
    initZero_(false),
    debug_(false),
    relTol_(1.),
    solverName_("GMRES"),
    solverParameter_(
        createLinSolverParameter("GMRES", 1.e-1, -1, Teuchos::rcp(new Teuchos::oblackholestream), 200)),
    problem_(Teuchos::rcp(new Belos::LinearProblem<ST, MF, MOpT> ())) {}


  InverseOp(const Teuchos::RCP<OpT>& op):
    level_(false),
    nullspaceOrtho_(false),
    initZero_(false),
    debug_(false),
    relTol_(1.),
    solverName_("GMRES"),
    solverParameter_(
        createLinSolverParameter("GMRES", 1.e-1, -1, Teuchos::rcp(new Teuchos::oblackholestream), 10)),
    problem_(
      Teuchos::rcp(
        new Belos::LinearProblem<ST, MF, MOpT>(
          createMultiOperatorBase(op),
          Teuchos::null,
          Teuchos::null))) {
      if(nullspaceOrtho_)
        projector_ = Teuchos::rcp(new ProjectorT<OpT>(op));
    }


  //template<class IOperatorT>
  InverseOp(const Teuchos::RCP<OperatorT>& op,
             const Teuchos::RCP<Teuchos::ParameterList>& pl):
    level_(pl->get<bool>("level", false)),
    nullspaceOrtho_(pl->get<bool>("nullspace ortho", false)),
    initZero_(pl->get<bool>("initZero", false)),
    debug_(pl->get<bool>("debug", false)),
    relTol_(pl->get<ST>("relative Tol", 1.)),
    solverName_(pl->get<std::string>("Solver name", "GMRES")),
    solverParameter_(Teuchos::sublist(pl, "Solver")),
    problem_(
      Teuchos::rcp(new Belos::LinearProblem<ST, MF, MOpT>(
                      createMultiOperatorBase(op),
                      Teuchos::null,
                      Teuchos::null))) {
      if(nullspaceOrtho_)
        projector_ = Teuchos::rcp(new ProjectorT<OpT>(op));
    }


  void apply(const DomainFieldT& rhs, RangeFieldT& x, const Add add=Add::N) const {
    apply(
        wrapMultiField(Teuchos::rcpFromRef<DomainFieldT>(const_cast<DomainFieldT&>(rhs))),
        wrapMultiField(Teuchos::rcpFromRef<RangeFieldT>(x)));
  }


  /// \brief MultiField helper (useful for MH ops)
  void apply(const MF& rhs, MF& x, const Add add=Add::N) const {

    apply(Teuchos::rcpFromRef<const MF>(rhs), Teuchos::rcpFromRef<MF>(x));
  }

private:

  void apply(const Teuchos::RCP<const MF>& rhs, const Teuchos::RCP<MF>& x, const Add add=Add::N) const {

    if(initZero_) x->init();

    if(nullspaceOrtho_) {
      for(int i=0; i<rhs->getNumberVecs(); ++i)
        (*projector_)(const_cast<RangeFieldT&>(rhs->getField(i)));
    }

    problem_->setProblem(x, rhs);
    Belos::SolverFactory<ST, MF, MOpT> factory;

    Teuchos::RCP<Belos::SolverManager<ST, MF, MOpT> > solver =
      factory.create(solverName_, solverParameter_);

    solver->setProblem(problem_);
    Belos::ReturnType succes = solver->solve();

    if(debug_) {
      x->write();
      rhs->write(100);
      std::cout << getLabel() << ": succes? " << succes << "\n";
      assert(Belos::ReturnType::Converged==succes);
    }

    if(level_) x->level();
    problem_->setProblem();
  }

public:

  void assignField(const DomainFieldT& mvr) {

    auto mv = wrapMultiField(
                Teuchos::rcpFromRef<DomainFieldT>(const_cast<DomainFieldT&>(mvr)));

    Teuchos::rcp_const_cast<MOpT>(problem_->getOperator())->assignField(*mv);

    if(problem_->isLeftPrec()) {
      auto opPrec = Teuchos::rcp_const_cast<MOpT>(problem_->getLeftPrec());
      opPrec->assignField(*mv);
    }

    if(problem_->isRightPrec()) {
      auto opPrec = Teuchos::rcp_const_cast<MOpT>(problem_->getRightPrec());
      opPrec->assignField(*mv);
    }
  };


  constexpr const Teuchos::RCP<const GridT>& grid() {
    return problem_->getOperator()->grid();
  };


  constexpr Teuchos::RCP<const MOpT> getOperator() {
    return problem_->getOperator();
  };


  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    if(para->name()=="Linear Solver")
      solverParameter_->set<ST>("Convergence Tolerance",
                                 para->get<double>("Convergence Tolerance")*relTol_);

    Teuchos::rcp_const_cast<MOpT>(problem_->getOperator())->setParameter(para);

    if(problem_->isLeftPrec())
      Teuchos::rcp_const_cast<MOpT>(problem_->getLeftPrec())->setParameter(para);

    if(problem_->isRightPrec())
      Teuchos::rcp_const_cast<MOpT>(problem_->getRightPrec())->setParameter(para);
  }


  /// \brief Set left preconditioner (\c LP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setLeftPrec(const Teuchos::RCP<const MOpT>& LP) {
    problem_->setLeftPrec(LP);
  }

  /// \brief Set right preconditioner (\c RP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setRightPrec(const Teuchos::RCP<const MOpT>& RP) {
    problem_->setRightPrec(RP);
  }

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return problem_->getOperator()->getLabel() + std::string("^-1 ");
  };

  void print(std::ostream& out=std::cout) const {
    out << "Inverse:\n";
    problem_->getOperator()->print(out);
  }


}; // end of class InverseOp



/// \relates InverseOp
template<class OpT>
Teuchos::RCP<InverseOp<OpT> >
createInverseOp(
  const Teuchos::RCP<OpT>& op,
  const Teuchos::RCP<Teuchos::ParameterList>& pl) {

  return Teuchos::rcp(new InverseOp<OpT>(op, pl));
}



/// \relates InverseOp
template<template<class> class  ProjectorT, class OpT  >
Teuchos::RCP<InverseOp<OpT, ProjectorT> >
createInverseOp(
  const Teuchos::RCP<OpT>& op,
  const Teuchos::RCP<Teuchos::ParameterList>& pl) {

  return Teuchos::rcp(new InverseOp<OpT, ProjectorT>(op, pl));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSEOP_HPP
