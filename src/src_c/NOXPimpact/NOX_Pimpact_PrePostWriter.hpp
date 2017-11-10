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

    Teuchos::rcp_const_cast<FieldT>(
        Teuchos::rcp_dynamic_cast<const FieldT>(
          group.getXPtr() ))->getField().write(sol_);
  };


  void writeRes(const NOX::Solver::Generic& solver) {

    const NOX::Abstract::Group& group = solver.getSolutionGroup();

    Teuchos::rcp_const_cast<FieldT>(
        Teuchos::rcp_dynamic_cast<const FieldT>(
          group.getFPtr() ))->getField().write(res_);
  };


  void writeCor(const NOX::Solver::Generic& solver) {

    const NOX::Abstract::Group& group = solver.getSolutionGroup();

    Teuchos::rcp_const_cast<FieldT>(
        Teuchos::rcp_dynamic_cast<const FieldT>(
          group.getNewtonPtr() ))->getField().write(cor_);
  };


  int sol_;
  int res_;
  int cor_;

  bool writeIterPre_;
  bool writeIterPost_;

  bool writeSolvePre_;
  bool writeSolvePost_;

public:

  PrePostWriter():
    sol_(-1),
    res_(-1),
    cor_(-1),
    writeIterPre_(false),
    writeIterPost_(false),
    writeSolvePre_(false),
    writeSolvePost_(false) {};

  PrePostWriter(const Teuchos::RCP<Teuchos::ParameterList>& pl):
    sol_( pl->get<int>("solution", -1 )),
    res_( pl->get<int>("residual", -1 )),
    cor_( pl->get<int>("update", -1 )),
    writeIterPre_(  pl->sublist("iteration").get<bool>("pre", false) ),
    writeIterPost_( pl->sublist("iteration").get<bool>("post", false) ),
    writeSolvePre_(  pl->sublist("solve").get<bool>("pre", false) ),
    writeSolvePost_(  pl->sublist("solve").get<bool>("post", false) ) {

      //if( (sol_<0 or res_<0 or cor_<0) and ( !writeIterPre_ and !writeIterPost_ and
            //!writeSolvePre_ and !writeSolvePost_ ) ) 
        //std::cout << "Warning!!! no write happening!\n";
    };


  PrePostWriter(const NOX::Abstract::PrePostOperator& ) = delete;

  PrePostWriter(NOX::Abstract::PrePostOperator&& that) {
    NOX::Pimpact::PrePostWriter<FieldT>& that_ (that);

    sol_ = that_.sol_;
    res_ = that_.res_;
    cor_ = that_.cor_;

    writeIterPre_  = that_.writeIterPre_;
    writeIterPost_ = that_.writeIterPost_;

    writeSolvePre_  =  that_.writeSolvePre_;
    writeSolvePost_ =  that_.writeSolvePost_;
  }


  virtual ~PrePostWriter() {};

  virtual void runPreIterate(const NOX::Solver::Generic& solver) {
    if( writeIterPre_ ) {
      if( sol_>=0 ) writeSol(solver);
      if( res_>=0 ) writeRes(solver);
      if( cor_>=0 ) writeCor(solver);
    }
  }

  virtual void runPostIterate(const NOX::Solver::Generic& solver) {
    if( writeIterPost_ ) {
      if( sol_>=0 ) writeSol(solver);
      if( res_>=0 ) writeRes(solver);
      if( cor_>=0 ) writeCor(solver);
    }
  }

  virtual void runPreSolve(const NOX::Solver::Generic& solver) {
    if( writeSolvePre_ ) {
      if( sol_>=0 ) writeSol(solver);
      if( res_>=0 ) writeRes(solver);
      if( cor_>=0 ) writeCor(solver);
    }
  }

  virtual void runPostSolve(const NOX::Solver::Generic& solver) {
    if( writeSolvePost_ ) {
      if( sol_>=0 ) writeSol(solver);
      if( res_>=0 ) writeRes(solver);
      if( cor_>=0 ) writeCor(solver);
    }
  }

}; // class PrePostWriter

} // namespace Abstract
} // namespace NOX


#endif // end of #ifndef NOX_PIMPACT_PREPOSTWRITER_HPP
