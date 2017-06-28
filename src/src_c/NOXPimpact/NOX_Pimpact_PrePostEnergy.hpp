#pragma once
#ifndef NOX_PIMPACT_PREPOSTENERGYCOMPUTE_HPP
#define NOX_PIMPACT_PREPOSTENERGYCOMPUTE_HPP

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
class PrePostEnergyCompute : public NOX::Abstract::PrePostOperator {

  //using FieldT = typename NV::FieldT;
  using FieldT = ::Pimpact::VectorField<typename NV::FieldT::SpaceT>;


  using ST = typename FieldT::SpaceT::Scalar;
  using OT = typename FieldT::SpaceT::Ordinal;

  bool computeEnergyIterPost_;
  bool computeEnergyIterPre_;
  bool computeEnergySolPost_;
	bool computeEnergySolPre_;

  int refinement_;

  ST gamma_;
  ::Pimpact::ECoord dir_;

	Teuchos::RCP<FieldT> base_;

	Teuchos::RCP<std::ostream> eStream_;


	void computeEnergy(const NOX::Solver::Generic& solver) {

    int nIter = solver.getNumIterations();

    const NOX::Abstract::Group& group = solver.getSolutionGroup();

    auto& x = Teuchos::rcp_const_cast<NV>( Teuchos::rcp_dynamic_cast<const NV>(
            group.getXPtr() ))->getField();
    auto space = x.space();

    std::string prefix = "energy_"+ ::Pimpact::toString(dir_) +"_r" +
      std::to_string(refinement_) + "_i" + std::to_string(nIter) + "_";

    // compute glob energy in y-dir
    if( 0==space->si(::Pimpact::F::U,3) ) {
      auto vel =  x.getField(0).getVField().get0Field().clone( ::Pimpact::ECopy::Deep );
      vel->add( 1., *vel, -1., *base_ );

      auto out = ::Pimpact::createOstream(
           prefix + "0.txt",
          space->getProcGrid()->getRankBar(dir_));

      ::Pimpact::computeEnergyDir( *vel, *out, dir_, gamma_ );
    }

    for( OT i=std::max(space->si(::Pimpact::F::U,3),1); i<=space->ei(::Pimpact::F::U,3); ++i ) {
      {
        auto out = ::Pimpact::createOstream( prefix + "C"+std::to_string(i) + ".txt",
            space->getProcGrid()->getRankBar(dir_) );

        ::Pimpact::computeEnergyDir(
            x.getField(0).getVField().getCField(i), *out, dir_, gamma_ );
      }
      {
        auto out = ::Pimpact::createOstream( prefix + "S"+std::to_string(i) + ".txt",
            space->getProcGrid()->getRankBar(dir_) );

        ::Pimpact::computeEnergyDir(
            x.getField(0).getVField().getSField(i), *out, dir_, gamma_ );
      }
    }
  };

public:

	PrePostEnergyCompute():
		computeEnergyIterPost_(false),
		computeEnergyIterPre_(false),
		computeEnergySolPost_(false),
		computeEnergySolPre_(false),
    refinement_(0),
    gamma_(0.),
    dir_(::Pimpact::ECoord::X),
		base_(Teuchos::null),
		eStream_(Teuchos::null) {};

	PrePostEnergyCompute(const Teuchos::RCP<Teuchos::ParameterList>& pl,
			Teuchos::RCP<FieldT> sol=Teuchos::null ):
    computeEnergyIterPost_( pl->sublist("iter").get<bool>("post", false) ),
		computeEnergyIterPre_(  pl->sublist("iter").get<bool>("pre", false) ),
		computeEnergySolPost_(  pl->sublist("solv").get<bool>("post", false) ),
		computeEnergySolPre_(   pl->sublist("solv").get<bool>("pre", false) ),
    refinement_(pl->get<int>("refinement", 0)),
    gamma_( pl->get<ST>("gamma", 0.) ),
    dir_( static_cast<::Pimpact::ECoord>( pl->get<int>("dir", 1) ) ),
		base_( sol ),
		eStream_(Teuchos::null) {

			if( base_!=Teuchos::null ) {
				int world_rank = sol->space()->rankST();
        if( 0==world_rank ) 
          eStream_ = Teuchos::rcp( new std::ofstream( "errorIter.txt" ) );
        else
          eStream_ = Teuchos::rcp( new Teuchos::oblackholestream );
			}
	};

	PrePostEnergyCompute(const NOX::Abstract::PrePostOperator& ) = delete;

  PrePostEnergyCompute(NOX::Abstract::PrePostOperator&& that) : PrePostEnergyCompute() {

    NOX::Pimpact::PrePostEnergyCompute<NV>& that_ (that);

    computeEnergyIterPost_ = that_.computeEnergyIterPost_;
    computeEnergyIterPre_  = that_.computeEnergyIterPre_;

    computeEnergySolPost_ = that_.computeEnergySolPost_;
    computeEnergySolPre_  = that_.computeEnergySolPre_;

    refinement_ = that_.refinement_;

    gamma_ = that_.gamma_;
    dir_ = that_.dir_;

    base_.swap( that_.base_);
    eStream_.swap( that_.eStream_ );
  }


	virtual ~PrePostEnergyCompute() {
    base_= Teuchos::null;
    eStream_ = Teuchos::null;
  };


	virtual void runPreIterate(const NOX::Solver::Generic& solver) {
		if( computeEnergyIterPre_ ) computeEnergy(solver);
	}

	virtual void runPostIterate(const NOX::Solver::Generic& solver) {
		if( computeEnergyIterPost_ ) computeEnergy(solver);
	}

	virtual void runPreSolve(const NOX::Solver::Generic& solver) {
		if( computeEnergySolPre_ ) computeEnergy(solver);
	}

	virtual void runPostSolve(const NOX::Solver::Generic& solver) {
		if( computeEnergySolPost_ ) computeEnergy(solver);
	}

}; // class PrePostEnergyCompute

} // namespace Abstract
} // namespace NOX


#endif // end of #ifndef NOX_PIMPACT_PREPOSTENERGYCOMPUTE_HPP
