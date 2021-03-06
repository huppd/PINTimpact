/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef NOX_PIMPACT_INTERFACE_HPP
#define NOX_PIMPACT_INTERFACE_HPP


#include <mpi.h>
#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_OperatorBase.hpp"




namespace NOX {
namespace Pimpact {

/// \brief Provides a set of interfaces for users to provide information about the nonlinear problem to NOX.
/// Contains interfaces for the user to supply
/// (1) the evaluation of the nonlinear equations,
/// (2) the Jacobian, and
/// (3) any preconditioning if required.

/// \brief Supplies NOX with the set nonlinear equations.
///
/// This is the minimum required information to solve a nonlinear
/// problem using the NOX::Epetra objects for the linear algebra
/// implementation.  Used by NOX::Epetra::Group to provide a link
/// to the external code for residual fills.
/// \tparam FieldT hast to be of type \c Pimpact::MultiField.
/// \note it is assumed that op_ is getting assigned if jopInv is assigned
template<class FT, class OpT=::Pimpact::OperatorBase<FT>, class IOpT=::Pimpact::OperatorBase<FT> >
class Interface {

public:

  using FieldT = FT;

  using OperatorT = OpT;

protected:

  Teuchos::RCP<FieldT> fu_;
  Teuchos::RCP<OpT>    op_;
  Teuchos::RCP<IOpT>   jopInv_;
  Teuchos::RCP<FieldT> resScaling_;

public:

  /// Constructor
  Interface(
    Teuchos::RCP<FieldT> fu=Teuchos::null,
    Teuchos::RCP<OpT>    op=Teuchos::null,
    Teuchos::RCP<IOpT>   jop=Teuchos::null,
    Teuchos::RCP<FieldT> resScaling=Teuchos::null):
    fu_(fu),
    op_(op),
    jopInv_(jop),
    resScaling_(resScaling) {};


  /// Compute the function, F, given the specified input vector x.
  NOX::Abstract::Group::ReturnType computeF(const FieldT& x, FieldT& f) {

    op_->apply(x, f);
    f.add(1., f, -1., *fu_);

    return NOX::Abstract::Group::Ok;
  }


  /// \brief Compute the Jacobian Operator, given the specified input vector x.
  NOX::Abstract::Group::ReturnType computeJacobian(const FieldT& x) {

    jopInv_->assignField(x);
    return NOX::Abstract::Group::Ok;
  }


  NOX::Abstract::Group::ReturnType applyJacobian(const FieldT& x, FieldT& y,
      const Belos::ETrans type=Belos::NOTRANS) {

    return NOX::Abstract::Group::NotDefined;
  }


  NOX::Abstract::Group::ReturnType applyJacobianInverse(Teuchos::ParameterList &params,
      const FieldT& x, FieldT& y) {

    double atol = params.get<double>("Tolerance");
    if(atol>0.) {
      double res = x.norm();

      Teuchos::RCP<Teuchos::ParameterList> para = Teuchos::parameterList("Linear Solver");
      para->set<double>("Convergence Tolerance", atol*res);

      jopInv_->setParameter(para);
    }

    jopInv_->apply(x, y);
    return NOX::Abstract::Group::Ok;
  }


  NOX::Abstract::Group::ReturnType applyPreconditioner(const FieldT& x, FieldT& y) {
    return NOX::Abstract::Group::NotDefined;
  }


  Teuchos::RCP<const OpT> getOperatorPtr() const {
    return op_;
  }

  Teuchos::RCP<const FieldT> getScalingPtr() const {
    return resScaling_;
  }

}; // end of class Interface



/// \relates Interface
template<class FT, class OpT=::Pimpact::OperatorBase<FT>, class IOpT=::Pimpact::OperatorBase<FT> >
Teuchos::RCP<Interface<FT, OpT, IOpT> > createInterface(
  Teuchos::RCP<FT> fu=Teuchos::null,
  Teuchos::RCP<OpT>  op=Teuchos::null,
  Teuchos::RCP<IOpT> jop=Teuchos::null,
  Teuchos::RCP<FT> scaleField=Teuchos::null) {

  return Teuchos::rcp(new Interface<FT, OpT, IOpT>(fu, op, jop, scaleField));
}


} // end of namespace Pimpact
} // end of namespace NOX


#endif // end of #ifndef NOX_PIMPACT_INTERFACE_HPP
