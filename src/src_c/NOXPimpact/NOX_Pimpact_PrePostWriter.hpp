#pragma once
#ifndef NOX_PIMPACT_PREPOSTWRITER_HPP
#define NOX_PIMPACT_PREPOSTWRITER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "NOX_Common.H"  // for NOX_Config.h
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Solver_Generic.H"



namespace NOX {
namespace Pimpact {


template<class FieldT>
class PrePostWriter : public NOX::Abstract::PrePostOperator {

  void writeSol(const NOX::Solver::Generic& solver) {

    const NOX::Abstract::Group& group = solver.getSolutionGroup();

    Teuchos::rcp_const_cast<FieldT>( Teuchos::rcp_dynamic_cast<const FieldT>(
                                       group.getXPtr() ))->getField().write();
  };

  void writeRes(const NOX::Solver::Generic& solver) {
    //Teuchos::RCP<NOX::Abstract::Group> group = solver.getSolutionGroup();

    //Teuchos::rcp_const_cast<FieldT>(Teuchos::rcp_dynamic_cast<const FieldT>(
    //group.getFPtr() ))->getField().write();
  };

  bool writeSolIterPost_;
  bool writeSolIterPre_;
  bool writeResIterPost_;
  bool writeResIterPre_;

  bool writeSolSolPost_;
  bool writeSolSolPre_;
  bool writeResSolPost_;
  bool writeResSolPre_;

public:

  PrePostWriter():
    writeSolIterPost_(false),
    writeSolIterPre_(false),
    writeResIterPost_(false),
    writeResIterPre_(false),
    writeSolSolPost_(false),
    writeSolSolPre_(false),
    writeResSolPost_(false),
    writeResSolPre_(false) {};

  PrePostWriter(const Teuchos::RCP<Teuchos::ParameterList>& pl):
    writeSolIterPost_( pl->sublist("sol").sublist("iter").get<bool>("post", false) ),
    writeSolIterPre_(  pl->sublist("sol").sublist("iter").get<bool>("pre", false) ),
    writeResIterPost_( pl->sublist("res").sublist("iter").get<bool>("post", false) ),
    writeResIterPre_(  pl->sublist("res").sublist("iter").get<bool>("pre", false) ),
    writeSolSolPost_(  pl->sublist("sol").sublist("solv").get<bool>("post", false) ),
    writeSolSolPre_(   pl->sublist("sol").sublist("solv").get<bool>("pre", false) ),
    writeResSolPost_(  pl->sublist("res").sublist("solv").get<bool>("post", false) ),
    writeResSolPre_(   pl->sublist("res").sublist("solv").get<bool>("pre", false) ) {};

  PrePostWriter(const NOX::Abstract::PrePostOperator& ) = delete;

  PrePostWriter(NOX::Abstract::PrePostOperator&& that) {
    NOX::Pimpact::PrePostWriter<FieldT>& that_ (that);

    writeSolIterPost_ = that_.writeSolIterPost_;
    writeSolIterPre_  = that_.writeSolIterPre_;
    writeResIterPost_ = that_.writeResIterPost_;
    writeResIterPre_  = that_.writeResIterPre_;

    writeSolSolPost_ =  that_.writeSolSolPost_;
    writeSolSolPre_  =  that_.writeSolSolPre_;
    writeResSolPost_ =  that_.writeResSolPost_;
    writeResSolPre_  =  that_.writeResSolPre_;
  }


  virtual ~PrePostWriter() {};

  virtual void runPreIterate(const NOX::Solver::Generic& solver) {
    if( writeSolIterPre_ ) writeSol(solver);
    if( writeResIterPre_ ) writeRes(solver);
  }

  virtual void runPostIterate(const NOX::Solver::Generic& solver) {
    if( writeSolIterPost_ ) writeSol(solver);
    if( writeResIterPost_ ) writeRes(solver);
  }

  virtual void runPreSolve(const NOX::Solver::Generic& solver) {
    if( writeSolSolPre_ ) writeSol(solver);
    if( writeResSolPre_ ) writeRes(solver);
  }

  virtual void runPostSolve(const NOX::Solver::Generic& solver) {
    if( writeSolSolPost_ ) writeSol(solver);
    if( writeResSolPost_ ) writeRes(solver);
  }

}; // class PrePostWriter

} // namespace Abstract
} // namespace NOX


#endif // end of #ifndef NOX_PIMPACT_PREPOSTWRITER_HPP
