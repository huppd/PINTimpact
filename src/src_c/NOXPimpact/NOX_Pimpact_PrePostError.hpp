#pragma once
#ifndef NOX_PIMPACT_PREPOSTERRORCOMPUTE_HPP
#define NOX_PIMPACT_PREPOSTERRORCOMPUTE_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "NOX_Common.H"  // for NOX_Config.h
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Solver_Generic.H"



namespace NOX {
namespace Pimpact {


template<class NV>
class PrePostErrorCompute : public NOX::Abstract::PrePostOperator {

  using FieldT = typename NV::FieldT;

  bool computeErrorIterPost_;
  bool computeErrorIterPre_;
  bool computeErrorSolPost_;
	bool computeErrorSolPre_;

	Teuchos::RCP<FieldT> sol_;

	Teuchos::RCP<std::ostream> eStream_;


	void computeError(const NOX::Solver::Generic& solver) {

		const NOX::Abstract::Group& group = solver.getSolutionGroup();

    auto& x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>(
						group.getXPtr()))->getField();

		if(!sol_.is_null()) {
			auto temp = x.getField(0).getVField().clone();

      double solNorm = sol_->getField(0).getVField().norm();

      temp->add(1., sol_->getField(0).getVField(), -1., x.getField(0).getVField());

      if(solNorm!=0.)
        *eStream_ << solver.getNumIterations() << "\t" << temp->norm()/solNorm << "\n";
		}
	};

public:

	PrePostErrorCompute():
		computeErrorIterPost_(false),
		computeErrorIterPre_(false),
		computeErrorSolPost_(false),
		computeErrorSolPre_(false),
		sol_(Teuchos::null),
		eStream_(Teuchos::null) {};

	PrePostErrorCompute(const Teuchos::RCP<Teuchos::ParameterList>& pl,
			Teuchos::RCP<FieldT> sol=Teuchos::null):

    computeErrorIterPost_(pl->sublist("iter").get<bool>("post", false)),
		computeErrorIterPre_(pl->sublist("iter").get<bool>("pre", false)),
		computeErrorSolPost_(pl->sublist("solv").get<bool>("post", false)),
		computeErrorSolPre_(pl->sublist("solv").get<bool>("pre", false)),
		sol_(sol),
		eStream_(Teuchos::null) {

			if(sol_!=Teuchos::null) {
				int world_rank = sol->grid()->rankST();
        if(0==world_rank) 
          eStream_ = Teuchos::rcp(new std::ofstream("errorIter.txt"));
        else
          eStream_ = Teuchos::rcp(new Teuchos::oblackholestream);
			}
	};

	PrePostErrorCompute(const NOX::Abstract::PrePostOperator&) = delete;

  PrePostErrorCompute(NOX::Abstract::PrePostOperator&& that) : PrePostErrorCompute() {

    NOX::Pimpact::PrePostErrorCompute<FieldT>& that_ (that);

    computeErrorIterPost_ = that_.computeErrorIterPost_;
    computeErrorIterPre_  = that_.computeErrorIterPre_;

    computeErrorSolPost_ = that_.computeErrorSolPost_;
    computeErrorSolPre_  = that_.computeErrorSolPre_;

    sol_.swap(that_.sol_);
    eStream_.swap(that_.eStream_);
  }


	virtual ~PrePostErrorCompute() {
    sol_= Teuchos::null;
    eStream_ = Teuchos::null;
  };


	virtual void runPreIterate(const NOX::Solver::Generic& solver) {
		if(computeErrorIterPre_) computeError(solver);
	}

	virtual void runPostIterate(const NOX::Solver::Generic& solver) {
		if(computeErrorIterPost_) computeError(solver);
	}

	virtual void runPreSolve(const NOX::Solver::Generic& solver) {
		if(computeErrorSolPre_) computeError(solver);
	}

	virtual void runPostSolve(const NOX::Solver::Generic& solver) {
		if(computeErrorSolPost_) computeError(solver);
	}

}; // class PrePostErrorCompute

} // namespace Abstract
} // namespace NOX


#endif // end of #ifndef NOX_PIMPACT_PREPOSTERRORCOMPUTE_HPP
