#pragma once
#ifndef NOX_PIMPACT_PREPOSTSPECTRUM_HPP
#define NOX_PIMPACT_PREPOSTSPECTRUM_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "NOX_Common.H"  // for NOX_Config.h
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Solver_Generic.H"

#include "Pimpact_Utils.hpp"
#include "Pimpact_AnalysisTools.hpp"


namespace NOX {
namespace Pimpact {


template<class NV>
class PrePostSpectrum : public NOX::Abstract::PrePostOperator {

  using FieldT = typename NV::FieldT;

  using ST = typename FieldT::SpaceT::Scalar;
  using OT = typename FieldT::SpaceT::Ordinal;

  bool sol_;
  bool res_;
  bool cor_;

  bool iterPost_;
  bool iterPre_;
  bool solPost_;
	bool solPre_;

  int refinement_;


	void writeSpectrum(const NOX::Solver::Generic& solver) {

    int nIter = solver.getNumIterations();

    const NOX::Abstract::Group& group = solver.getSolutionGroup();


    if(sol_) {
      auto& x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>(
            group.getXPtr()))->getField();

      {
        auto out = ::Pimpact::createOstream(
            "xv_"+std::to_string(refinement_)+"_"+std::to_string(nIter)+".txt",
            x.space()->rankST());
        ::Pimpact::writeSpectrum(x.getField(0).getVField(), *out);
      }

      {
        auto out = ::Pimpact::createOstream(
            "xp_"+std::to_string(refinement_)+"_"+std::to_string(nIter)+".txt",
            x.space()->rankST());
        ::Pimpact::writeSpectrum(x.getField(0).getSField(), *out);
      }
    }

    if(res_) {
      auto& x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>(
            group.getFPtr()))->getField();
      {
        auto out = ::Pimpact::createOstream(
            "resv_"+std::to_string(refinement_)+"_"+std::to_string(nIter)+".txt",
            x.space()->rankST());
        ::Pimpact::writeSpectrum(x.getField(0).getVField(), *out);
      }
      {
        auto out = ::Pimpact::createOstream(
            "resp_"+std::to_string(refinement_)+"_"+std::to_string(nIter)+".txt",
            x.space()->rankST());
        ::Pimpact::writeSpectrum(x.getField(0).getSField(), *out);
      }
    }

    if(cor_) {
      auto& x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>(
            group.getNewtonPtr()))->getField();
      {
        auto out = ::Pimpact::createOstream(
            "corv_"+std::to_string(refinement_)+"_"+std::to_string(nIter)+".txt",
            x.space()->rankST());
        ::Pimpact::writeSpectrum(x.getField(0).getVField(), *out);
      }
      {
        auto out = ::Pimpact::createOstream(
            "corp_"+std::to_string(refinement_)+"_"+std::to_string(nIter)+".txt",
            x.space()->rankST());
        ::Pimpact::writeSpectrum(x.getField(0).getSField(), *out);
      }
    }
  };

public:

	PrePostSpectrum():
    sol_(false),
    res_(false),
    cor_(false),
		iterPost_(false),
		iterPre_(false),
		solPost_(false),
		solPre_(false),
		refinement_(0) {};

	PrePostSpectrum(const Teuchos::RCP<Teuchos::ParameterList>& pl):
    sol_(pl->get<bool>("solution", false)),
    res_(pl->get<bool>("residual", false)),
    cor_(pl->get<bool>("update", false)),
    iterPost_(pl->sublist("iteration").get<bool>("post", false)),
		iterPre_(pl->sublist("iteration").get<bool>("pre", false)),
		solPost_(pl->sublist("solve").get<bool>("post", false)),
		solPre_(pl->sublist("solve").get<bool>("pre", false)),
		refinement_(pl->get<int>("refinement", 0)) {

	};

	PrePostSpectrum(const NOX::Abstract::PrePostOperator&) = delete;

  PrePostSpectrum(NOX::Abstract::PrePostOperator&& that) : PrePostSpectrum() {

    NOX::Pimpact::PrePostSpectrum<NV>& that_ (that);

    iterPost_ = that_.iterPost_;
    iterPre_  = that_.iterPre_;

    solPost_ = that_.solPost_;
    solPre_  = that_.solPre_;

    refinement_ = that_.refinement_;
  }


  virtual ~PrePostSpectrum() {};


	virtual void runPreIterate(const NOX::Solver::Generic& solver) {
		if(iterPre_) writeSpectrum(solver);
	}

	virtual void runPostIterate(const NOX::Solver::Generic& solver) {
		if(iterPost_) writeSpectrum(solver);
	}

	virtual void runPreSolve(const NOX::Solver::Generic& solver) {
		if(solPre_) writeSpectrum(solver);
	}

	virtual void runPostSolve(const NOX::Solver::Generic& solver) {
		if(solPost_) writeSpectrum(solver);
	}

}; // class PrePostSpectrum

} // namespace Abstract
} // namespace NOX


#endif // end of #ifndef NOX_PIMPACT_PREPOSTSPECTRUM_HPP
