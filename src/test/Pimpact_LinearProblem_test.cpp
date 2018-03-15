#include <iostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"

#include "BelosPimpactAdapter.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"




namespace {

using ST = double;
using OT = int;
const int sd = 3;

using GridT = Pimpact::Grid<ST, OT, sd, 3, 4>;

bool testMpi = true;
double eps = 1e-6;

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

	Pimpact::setBoundaryConditions(pl, 0);

  pl->set("nx", 25);
  pl->set("ny", 17);
  pl->set("nz",  9);
}



TEUCHOS_UNIT_TEST(BelosSolver, HelmholtzMV) {

  using VF = Pimpact::VectorField<GridT>;
  using MVF = Pimpact::MultiField<VF>;

  using BOp = Pimpact::OperatorBase<MVF>;

  auto grid = Pimpact::create<GridT>(pl);

  auto x = Pimpact::wrapMultiField(Pimpact::create<Pimpact::VectorField>(grid));
  auto b = Pimpact::wrapMultiField(Pimpact::create<Pimpact::VectorField>(grid));

  b->init(1.);

  auto op = Pimpact::createMultiOperatorBase(Pimpact::create<Pimpact::DiffusionOp>(grid));

  auto para = Pimpact::createLinSolverParameter("GMRES", eps);

  Belos::SolverFactory<ST, MVF, BOp> factory;

  // Create the GMRES solver.
  Teuchos::RCP<Belos::SolverManager<ST, MVF, BOp > > solver =
      factory.create("GMRES", para);

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  Teuchos::RCP<Belos::LinearProblem<ST, MVF, BOp > > problem =
      Teuchos::rcp (new Belos::LinearProblem<ST, MVF, BOp > (op, x, b));

  problem->setProblem(x, b);

  // Tell the solver what problem you want to solve.
  solver->setProblem(problem);

  solver->solve();
  TEST_EQUALITY(solver->achievedTol()<eps, true);

	if(solver->achievedTol()>=eps)
		x->write(111);

}



TEUCHOS_UNIT_TEST(BelosSolver, DivGrad) {

  using SF = Pimpact::ScalarField<GridT>;
  using BSF = Pimpact::MultiField<SF>;

  using BOp = Pimpact::OperatorBase<BSF>;


  auto grid = Pimpact::create<GridT>(pl);


  Pimpact::VectorField<GridT> temp(grid);
	temp(Pimpact::F::U).initField(Pimpact::Poiseuille2D_inX);

  auto x = Teuchos::rcp(new Pimpact::MultiField<Pimpact::ScalarField<GridT> >(grid, 1));
  auto b = x->clone();


	temp.init();
  auto op = Pimpact::createOperatorBase(
      Pimpact::createMultiOpWrap(
				Pimpact::create<Pimpact::DivGradOp<GridT> >(grid)
          //Pimpact::createDivGradOp(grid)
     )
 );

  //  x->init(1.);
  x->random();
  op->apply(*x, *b);
  x->init(0.);
  b->write(99);

  auto para = Pimpact::createLinSolverParameter("GMRES", eps);
  para->set("Num Blocks", 500);

  //  auto para = Teuchos::parameterList();


  Belos::SolverFactory<ST, BSF, BOp > factory;
  // Make an empty new parameter list.

  // Create the GMRES solver.
  Teuchos::RCP<Belos::SolverManager<ST, BSF, BOp > > solver =
      factory.create("GMRES", para);

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  Teuchos::RCP<Belos::LinearProblem<ST, BSF, BOp > > problem =
      Teuchos::rcp (new Belos::LinearProblem<ST, BSF, BOp > (op, x, b));

  //   std::cout <<"param\n" <<*solver->getValidParameters();
  problem->setProblem(x, b);

  // Tell the solver what problem you want to solve.
  solver->setProblem(problem);
  //       TEST_EQUALITY(solver->getProblem(), *problem)
  //       std::cout <<*problem;

  solver->solve();
  TEST_EQUALITY(solver->achievedTol()<eps, true);

	if(solver->achievedTol()>=eps) {
		x->write(999);
		temp.write(999);
	}

}


TEUCHOS_UNIT_TEST(LinearProblem, HelmholtzMV) {

	using MF = Pimpact::MultiField<Pimpact::VectorField<GridT> >;


	auto grid = Pimpact::create<GridT>(pl);

	auto x = Teuchos::rcp(new Pimpact::MultiField<Pimpact::VectorField<GridT> >(grid, 1));
	auto b = Teuchos::rcp(new Pimpact::MultiField<Pimpact::VectorField<GridT> >(grid, 1));

	x->init(0.);
	b->init(1.);

	auto A = Pimpact::createMultiOperatorBase(Pimpact::create<Pimpact::DiffusionOp>(grid));

	auto param = Pimpact::createLinSolverParameter("GMRES", 1.e-4);

	auto linprob = Pimpact::createLinearProblem<MF>(A, x, b, param, "GMRES");

	linprob->solve(x, b);

  TEST_EQUALITY(linprob->getSolver()->achievedTol()<1.e-4, true);
	if(linprob->getSolver()->achievedTol()>=1.e-4) 
		x->write();

}


} // end of namespace

