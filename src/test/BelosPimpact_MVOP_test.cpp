#include <iostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "BelosOutputManager.hpp"
#include "BelosMVOPTester.hpp"

#include "BelosPimpactAdapter.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"




namespace {


using ST = double;
using OT = int;
const int d = 3;
const int dNC = 4;

bool testMpi = true;
double eps = 1.e-6;
int domain = 0;


using GridT = typename Pimpact::Grid<ST, OT, 3, d, dNC>;

using SF = typename Pimpact::ScalarField<GridT>;
using VF = typename Pimpact::VectorField<GridT>;
using MSF = typename Pimpact::ModeField<SF>;
using MVF = typename Pimpact::ModeField<VF>;

using CF = Pimpact::CompoundField<VF, SF>;
using CMF = Pimpact::CompoundField<MVF, MSF>;


auto pl = Teuchos::parameterList();


TEUCHOS_STATIC_SETUP() {
	Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
	clp.addOutputSetupOptions(true);
	clp.setOption(
			"test-mpi", "test-serial", &testMpi,
			"Test MPI (if available) or force test of serial.  In a serial build, "
			" this option is ignored and a serial comm is always used.");
	clp.setOption(
			"error-tol-slack", &eps,
			"Slack off of machine epsilon used to check test results");
	clp.setOption(
			"domain", &domain,
			"domain");

	Pimpact::setBoundaryConditions(pl, domain);
	pl->set("nx", 33);
	pl->set("ny", 17);
	pl->set("nz", 9);

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TempField, BelosMVTest, FType) {

	auto grid = Pimpact::create<GridT>(pl);

	auto mx = Teuchos::rcp(new Pimpact::MultiField<FType>(grid, 5));

	// Create an output manager to handle the I/O from the solver
	Teuchos::RCP<Belos::OutputManager<ST> > MyOM =
		Teuchos::rcp(new Belos::OutputManager<ST>(Belos::Warnings, rcp(&out, false)));

	bool res = Belos::TestMultiVecTraits<ST, Pimpact::MultiField<FType> >(MyOM, mx);
	TEST_EQUALITY_CONST(res, true);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, BelosMVTest, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, BelosMVTest, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, BelosMVTest, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, BelosMVTest, MVF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, BelosMVTest, CF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, BelosMVTest, CMF)




TEUCHOS_UNIT_TEST(BelosOperatorMV, HelmholtzMV) {

	using VF = Pimpact::VectorField<GridT>;
	using MVF = Pimpact::MultiField<VF>;

	using OpBase = Pimpact::OperatorBase<MVF>;


	auto grid = Pimpact::create<GridT>(pl);

	auto mv = Teuchos::rcp(new Pimpact::MultiField<VF>(grid, 10));

	auto opm = Pimpact::createMultiOpWrap(Pimpact::create<Pimpact::DiffusionOp>(grid));

	auto op = Pimpact::createOperatorBase(opm);

	Teuchos::RCP<Belos::OutputManager<ST> > MyOM = Teuchos::rcp(
			new Belos::OutputManager<ST>(Belos::Warnings, rcp(&out, false)));

	bool res = Belos::TestOperatorTraits<ST, MVF, OpBase >(MyOM, mv, op);

	TEST_EQUALITY(res, true);
}




TEUCHOS_UNIT_TEST(BelosOperatorMV, DivGrad) {

	using SF = Pimpact::ScalarField<GridT>;
	using BSF = Pimpact::MultiField<SF>;

	using OpBase = Pimpact::OperatorBase<BSF>;

	auto grid = Pimpact::create<GridT>(pl);

	auto mv = Teuchos::rcp(new Pimpact::MultiField<Pimpact::ScalarField<GridT> >(grid , 10));
	auto mv2 = mv->clone();

	auto lap = Pimpact::createOperatorBase(
			Pimpact::createMultiOpWrap(
				Pimpact::createDivGradOp(
					Pimpact::create<Pimpact::DivOp>(grid),
					Pimpact::create<Pimpact::GradOp>(grid)
				)
			)
		);

	Teuchos::RCP<Belos::OutputManager<ST> > MyOM =
		Teuchos::rcp(new Belos::OutputManager<ST>(Belos::Errors + Belos::Warnings + Belos::IterationDetails +
					Belos::OrthoDetails + Belos::FinalSummary + Belos::TimingDetails +
					Belos::StatusTestDetails + Belos::Debug, rcp(&out, false)));


	bool res =// true;
	Belos::TestOperatorTraits<ST, BSF, OpBase > (MyOM, mv, lap);

	TEST_EQUALITY(res, true);
}



} // namespace

