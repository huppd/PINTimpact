#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinSolverParameter.hpp"

#include "BelosPimpactAdapter.hpp"





namespace {


typedef double S;
typedef int O;

typedef Pimpact::Space<S,O,3> SpaceT;

bool testMpi = true;
double errorTolSlack = 1e-4;

bool isImpactInit = false;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used." );
  clp.setOption(
      "error-tol-slack", &errorTolSlack,
      "Slack off of machine epsilon used to check test results" );
}


TEUCHOS_UNIT_TEST( BelosSolver, HelmholtzMV ) {

  typedef Pimpact::VectorField<SpaceT> VF;
  typedef Pimpact::MultiField<VF> MVF;

  typedef Pimpact::OperatorBase<MVF> BOp;

  auto pl = Teuchos::parameterList();

  pl->set( "domain", 1);

  auto space = Pimpact::createSpace( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto x = Pimpact::createMultiField( Pimpact::createVectorField(space) );
  auto b = Pimpact::createMultiField( Pimpact::createVectorField(space) );

  b->init( 1. );

  auto op = Pimpact::createMultiOperatorBase<MVF>( Pimpact::createHelmholtzOp(space) );

  auto para = Pimpact::createLinSolverParameter("GMRES",errorTolSlack);

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
  TEST_EQUALITY( solver->achievedTol()<errorTolSlack, true );

  x->write(111);

}




TEUCHOS_UNIT_TEST( BelosSolver, DtLapOp ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  typedef Pimpact::VectorField<SpaceT> VF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField<MVF> BVF;

  typedef Pimpact::OperatorBase<BVF> OpBase;

  typedef OpBase  BOp;


  auto pl = Teuchos::parameterList();

  pl->set( "domain", 1);

  auto space = Pimpact::createSpace( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto b = Pimpact::createInitMVF( Pimpact::Streaming2DFlow, space );
  auto x = Pimpact::createInitMVF( Pimpact::Zero2DFlow, space );


  auto op = Pimpact::createMultiOperatorBase<BVF>( Pimpact::createDtLapOp(space) );

  auto para = Pimpact::createLinSolverParameter("GMRES",errorTolSlack);

  Belos::SolverFactory<S, BVF, BOp > factory;
  // Make an empty new parameter list.

  // Create the GMRES solver.
  Teuchos::RCP<Belos::SolverManager<S, BVF, BOp > > solver =
      factory.create( "GMRES", para );

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  Teuchos::RCP<Belos::LinearProblem<S, BVF, BOp > > problem =
      Teuchos::rcp (new Belos::LinearProblem<S, BVF, BOp > (op, x, b));

  problem->setProblem(x,b);

  // Tell the solver what problem you want to solve.
  solver->setProblem( problem );

  solver->solve();
  TEST_EQUALITY( solver->achievedTol()<errorTolSlack, true );

  x->write(20);

}



TEUCHOS_UNIT_TEST( BelosSolver, DivGrad ) {

  typedef Pimpact::ScalarField<SpaceT> SF;
  typedef Pimpact::MultiField<SF> BSF;

  typedef Pimpact::DivGradOp<SpaceT> Op;
  typedef Pimpact::MultiOpWrap<Op> MuOp;
  typedef Pimpact::OperatorBase<BSF> OpBase;
  typedef OpBase  BOp;

  auto pl = Teuchos::parameterList();

  pl->set( "domain", 1);

  auto space = Pimpact::createSpace( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto temp = Pimpact::createVectorField( space );
  temp->initField(Pimpact::PoiseuilleFlow2D_inX);
  //  temp->init(0.);

  auto p = Pimpact::createScalarField(space);

  auto x = Pimpact::createMultiField<SF>(*p,1);
  auto b = x->clone();

//  auto div = Pimpact::createDivOp( space );

  temp->initField(Pimpact::ZeroFlow);
  auto op = Pimpact::createOperatorBase<BSF,MuOp>(
      Pimpact::createMultiOpWrap(
          Pimpact::createDivGradOp(
              temp,
              Pimpact::createDivOp( space ),
              Pimpact::createGradOp( space )
          )
      )
  );

  //  x->init( 1. );
  x->random();
  op->apply( *x, *b );
  x->init( 0. );
  b->write(99);

  auto para = Pimpact::createLinSolverParameter("GMRES",errorTolSlack);
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
  TEST_EQUALITY( solver->achievedTol()<errorTolSlack, true );

  x->write(999);
  temp->write(999);

}

} // end of namespace

