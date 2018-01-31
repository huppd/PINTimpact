#include <iostream>
#include <cmath>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_AnalysisTools.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_VectorField.hpp"

#include "Pimpact_Test.hpp"




namespace {





TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(VectorField, initField, SpaceT) {

	setParameter(SpaceT::sdim);

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(pl);

	space->getInterpolateV2S()->print();
	auto vel = Pimpact::create<Pimpact::VectorField>(space);
	auto divVec = Pimpact::create<Pimpact::ScalarField>(space);

	auto divOp = Pimpact::create<Pimpact::DivOp>(space);
	auto pl = Teuchos::parameterList();
	pl->set<std::string>("Type", "disturbance");

	vel->write();
	divOp->apply(*vel, *divVec);
	auto bla = divVec->norm(Pimpact::ENorm::Inf);
	if(0==space->rankST())
		std::cout <<"F: " << "\tmax div: " <<bla <<"\n";
	bla = divVec->norm(Pimpact::ENorm::Two);
	if(0==space->rankST())
		std::cout <<"F: " <<"\t||div||: " <<bla <<"\n";
	divVec->write(1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(VectorField, initField, D2) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(VectorField, initField, D3)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(VectorField, computeEnergy, SpaceT) {

	setParameter(SpaceT::sdim);

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(pl);

	Pimpact::VectorField<SpaceT> vel(space);
	vel.init(1.);

	std::cout <<"energy: " <<vel.norm(Pimpact::ENorm::L2) <<"\n";
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(VectorField, computeEnergy, D2) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(VectorField, computeEnergy, D3)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(VectorField, computeEnergyY, SpaceT) {

	setParameter(SpaceT::sdim);

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(pl);

	Pimpact::VectorField<SpaceT> vel(space);
	vel.init(1.);

  Pimpact::computeEnergyDir(vel);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(VectorField, computeEnergyY, D2) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(VectorField, computeEnergyY, D3)

} // end of namespace
