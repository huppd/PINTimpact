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

typedef double  S;
typedef int     O;

typedef Pimpact::Space<S,O,3,4> SpaceT;

bool testMpi = true;
double eps = 1e-6;

auto pl = Teuchos::parameterList();

TEUCHOS_STATIC_SETUP() {
	Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
	clp.addOutputSetupOptions(true);
	clp.setOption(
			"test-mpi", "test-serial", &testMpi,
			"Test MPI (if available) or force test of serial.  In a serial build,"
			" this option is ignored and a serial comm is always used." );
	clp.setOption(
			"error-tol-slack", &eps,
			"Slack off of machine epsilon used to check test results" );

  pl->set( "dim", 3);
  pl->set( "domain", 0);

  pl->set( "nx", 25 );
  pl->set( "ny", 17 );
  pl->set( "nz",  9 );

}



TEUCHOS_UNIT_TEST( BelosSolver, HelmholtzMV ) {

  typedef Pimpact::VectorField<SpaceT> VF;
  typedef Pimpact::MultiField<VF> MVF;

  typedef Pimpact::OperatorBase<MVF> BOp;

  auto space = Pimpact::createSpace( pl );

  auto x = Pimpact::createMultiField( Pimpact::create<Pimpact::VectorField>(space) );
  auto b = Pimpact::createMultiField( Pimpact::create<Pimpact::VectorField>(space) );

  b->init( 1. );

  auto op = Pimpact::createMultiOperatorBase( Pimpact::create<Pimpact::HelmholtzOp>( space ) );

  auto para = Pimpact::createLinSolverParameter("GMRES",eps);

  Belos::SolverFactory<S, MVF, BOp> factory;

  // Create the GMRES solver.
  Teuchos::RCP<Belos::SolverManager<S, MVF, BOp > > solver =
      factory.create( "GMRES", para );

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  Teuchos::RCP<Belos::LinearProblem<S, MVF, BOp > > problem =
      Teuchos::rcp (new Belos::LinearProblem<S, MVF, BOp > (op, x, b));

  problem->setProblem(x,b);

  // Tell the solver what problem you want to solve.
  solver->setProblem( problem );

  solver->solve();
  TEST_EQUALITY( solver->achievedTol()<eps, true );

	if( solver->achievedTol()>=eps )
		x->write(111);

}



//TEUCHOS_UNIT_TEST( BelosSolver, DtLapOp ) {
//
//
//  typedef Pimpact::VectorField<SpaceT> VF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField<MVF> BVF;
//
//  typedef Pimpact::OperatorBase<BVF> OpBase;
//
//  typedef OpBase  BOp;
//
//  auto space = Pimpact::createSpace( pl );
//
//  auto b = Pimpact::createInitMVF( Pimpact::Streaming2DFlow, space );
//  auto x = Pimpact::createInitMVF( Pimpact::Zero2DFlow, space );
//
//
//  auto op = Pimpact::createMultiOperatorBase( Pimpact::createDtLapOp(space) );
//
//  auto para = Pimpact::createLinSolverParameter("GMRES",eps);
//
//  Belos::SolverFactory<S, BVF, BOp > factory;
//  // Make an empty new parameter list.
//
//  // Create the GMRES solver.
//  Teuchos::RCP<Belos::SolverManager<S, BVF, BOp > > solver =
//      factory.create( "GMRES", para );
//
//  // Create a LinearProblem struct with the problem to solve.
//  // A, X, B, and M are passed by (smart) pointer, not copied.
//  Teuchos::RCP<Belos::LinearProblem<S, BVF, BOp > > problem =
//      Teuchos::rcp (new Belos::LinearProblem<S, BVF, BOp > (op, x, b));
//
//  problem->setProblem(x,b);
//
//  // Tell the solver what problem you want to solve.
//  solver->setProblem( problem );
//
//  solver->solve();
//  TEST_EQUALITY( solver->achievedTol()<eps, true );
//
//  x->write(20);
//
//}



TEUCHOS_UNIT_TEST( BelosSolver, DivGrad ) {

  typedef Pimpact::ScalarField<SpaceT> SF;
  typedef Pimpact::MultiField<SF> BSF;

  typedef Pimpact::OperatorBase<BSF>  BOp;


  auto space = Pimpact::createSpace( pl );


  auto temp = Pimpact::create<Pimpact::VectorField>( space );
  temp->initField(Pimpact::PoiseuilleFlow2D_inX);
  //  temp->init(0.);

  auto p = Pimpact::createScalarField(space);

  auto x = Pimpact::createMultiField( *p, 1 );
  auto b = x->clone();


  temp->initField(Pimpact::ZeroFlow);
  auto op = Pimpact::createOperatorBase(
      Pimpact::createMultiOpWrap(
				Pimpact::create< Pimpact::DivGradOp<SpaceT> >(space)
          //Pimpact::createDivGradOp(space)
      )
  );

  //  x->init( 1. );
  x->random();
  op->apply( *x, *b );
  x->init( 0. );
  b->write(99);

  auto para = Pimpact::createLinSolverParameter("GMRES",eps);
  para->set( "Num Blocks", 500 );

  //  auto para = Teuchos::parameterList();


  Belos::SolverFactory<S, BSF, BOp > factory;
  // Make an empty new parameter list.

  // Create the GMRES solver.
  Teuchos::RCP<Belos::SolverManager<S, BSF, BOp > > solver =
      factory.create( "GMRES", para );

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  Teuchos::RCP<Belos::LinearProblem<S, BSF, BOp > > problem =
      Teuchos::rcp (new Belos::LinearProblem<S, BSF, BOp > (op, x, b));

  //   std::cout << "param\n" << *solver->getValidParameters();
  problem->setProblem(x,b);

  // Tell the solver what problem you want to solve.
  solver->setProblem( problem );
  //       TEST_EQUALITY( solver->getProblem(), *problem)
  //       std::cout << *problem;

  solver->solve();
  TEST_EQUALITY( solver->achievedTol()<eps, true );

	if( solver->achievedTol()>=eps ) {
		x->write(999);
		temp->write(999);
	}

}


TEUCHOS_UNIT_TEST( LinearProblem, HelmholtzMV ) {

	typedef Pimpact::MultiField< Pimpact::VectorField<SpaceT> > MF;


	auto space = Pimpact::createSpace( pl );

	auto vel = Pimpact::create<Pimpact::VectorField>(space);

	auto x = Pimpact::createMultiField(*vel,1);
	auto b = Pimpact::createMultiField(*vel,1);

	x->init(0.);
	b->init(1.);

	auto A = Pimpact::createMultiOperatorBase( Pimpact::create<Pimpact::HelmholtzOp>( space ) );

	auto param = Pimpact::createLinSolverParameter("GMRES",1.e-4);

	auto linprob = Pimpact::createLinearProblem<MF>(A,x,b,param,"GMRES");

	linprob->solve(x,b);

  TEST_EQUALITY( linprob->getSolver()->achievedTol()<1.e-4, true );
	if( linprob->getSolver()->achievedTol()>=1.e-4 ) 
		x->write();

}


} // end of namespace

