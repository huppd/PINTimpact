/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef NOX_PIMPACT_PREPOSTWRITERESTART_HPP
#define NOX_PIMPACT_PREPOSTWRITERESTART_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "NOX_Common.H"  // for NOX_Config.h
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Solver_Generic.H"



namespace NOX {
namespace Pimpact {


template<class FieldT>
class PrePostWriteRestart : public NOX::Abstract::PrePostOperator {

  int writeIterPost_;

public:

  PrePostWriteRestart():
    writeIterPost_(-1) {};

  PrePostWriteRestart(const Teuchos::RCP<Teuchos::ParameterList>& pl):
    writeIterPost_(pl->get<int>("restart", -1)) {};


  PrePostWriteRestart(const NOX::Abstract::PrePostOperator&) = delete;

  PrePostWriteRestart(NOX::Abstract::PrePostOperator&& that) {
    NOX::Pimpact::PrePostWriteRestart<FieldT>& that_ (that);

    writeIterPost_ = that_.writeIterPost_;
  }

  virtual ~PrePostWriteRestart() {};

  virtual void runPreIterate(const NOX::Solver::Generic& solver) {}

  virtual void runPostIterate(const NOX::Solver::Generic& solver) {

    if(writeIterPost_!=-1 && (solver.getNumIterations()%writeIterPost_==0)) {
      const NOX::Abstract::Group& group = solver.getSolutionGroup();

      Teuchos::rcp_const_cast<FieldT>(Teuchos::rcp_dynamic_cast<const FieldT>(
            //group.getXPtr()))->getField().write(0, true);
            group.getXPtr()))->getField().write(100*solver.getNumIterations()/writeIterPost_, true);
    }
  }

  virtual void runPreSolve(const NOX::Solver::Generic& solver) {}

  virtual void runPostSolve(const NOX::Solver::Generic& solver) {}

}; // class PrePostWriteRestart

} // namespace Abstract
} // namespace NOX


#endif // end of #ifndef NOX_PIMPACT_PREPOSTWRITERESTART_HPP
