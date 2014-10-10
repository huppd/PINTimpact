#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include "Teuchos_Range1D.hpp"
#include <Teuchos_CommHelpers.hpp>
#include "BelosOutputManager.hpp"

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_Space.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorMV.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "BelosPimpactAdapter.hpp"
#include "Pimpact_LinSolverParameter.hpp"

#include "BelosSolverFactory.hpp"

#include <iostream>

namespace {

typedef double S;
typedef int O;

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

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::MultiField<VF> MVF;

//  typedef Pimpact::HelmholtzOp<S,O> Op;
//  typedef Pimpact::MultiOpWrap<Op> MOp;
  typedef Pimpact::OperatorBase<MVF> BOp;

	auto space = Pimpact::createSpace();

  auto x = Pimpact::createMultiField( Pimpact::createVectorField<S,O>(space) );
  auto b = Pimpact::createMultiField( Pimpact::createVectorField<S,O>(space) );

  b->init( 1. );

  auto op = Pimpact::createMultiOperatorBase<MVF>( Pimpact::createHelmholtzOp(space) );

  auto para = Pimpact::createLinSolverParameter("GMRES",1.e-3);

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



TEUCHOS_UNIT_TEST( BelosSolver, PrecDivGrad ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

//  typedef Pimpact::VectorField<S,O> VF;

  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::MultiField<SF> MSF;

//  typedef Pimpact::HelmholtzOp<S,O> Op;
  typedef Pimpact::OperatorBase<MSF> BOp;


	auto space = Pimpact::createSpace();

  auto u   = Pimpact::createVectorField<S,O>( space );
  auto temp= Pimpact::createVectorField<S,O>( space );

  auto rhs = Pimpact::createScalarField<S,O>( space );
  auto sol = Pimpact::createScalarField<S,O>( space );

  auto lap = Pimpact::createHelmholtzOp<S,O>( space, 0., 1. );

  auto div = Pimpact::createDivOp<S,O>( space );

  u->initField( Pimpact::PoiseuilleFlow2D_inX );
  u->write();


  lap->apply( *u, *temp );
  div->apply( *temp, *rhs );

  rhs->write(0);
  auto b = Pimpact::createMultiField( rhs );
  auto x = Pimpact::createMultiField( sol );
//  b->random();
  b->init(1.);
  x->init(0.);

  auto op = Pimpact::createMultiOperatorBase<MSF>(
      Pimpact::createDivGradOp<S,O>(
          u->clone(),
          Pimpact::createDivOp( space ),
          Pimpact::createGradOp( space )
      )
  );

  auto prec = Pimpact::createMultiOperatorBase<MSF>(
      Pimpact::createMGVDivGradOp<S,O>(true) );

  auto para = Pimpact::createLinSolverParameter("GMRES",1.e-6,1);
//  auto para = Teuchos::parameterList();
  para->set( "Maximum Iterations", 10000 );


  Belos::SolverFactory<S, MSF, BOp> factory;

  // Create the GMRES solver.
  Teuchos::RCP<Belos::SolverManager<S, MSF, BOp > > solver =
//      factory.create( "LSQR", para );
      factory.create( "GMRES", para );
//      factory.create( "TFQMR", para );
//      factory.create( "CG", para ); < not working

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  auto problem = Teuchos::rcp( new Belos::LinearProblem<S, MSF, BOp >( op, x, b) );

  problem->setProblem(x,b);
//  problem->setLeftPrec( prec );
  problem->setRightPrec( prec );

  // Tell the solver what problem you want to solve.
  solver->setProblem( problem );

//  Belos::ReturnType ret =
  solver->solve();
  x->write(222);

  TEST_EQUALITY( solver->achievedTol()<errorTolSlack, true );

}



TEUCHOS_UNIT_TEST( BelosSolver, DtLapOp ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField<MVF> BVF;

//  typedef Pimpact::DtLapOp<S,O> Op;
//  typedef Pimpact::MultiOpWrap<Op> OpMV;
  typedef Pimpact::OperatorBase<BVF> OpBase;
  //  typedef Pimpact::OperatorPimpldep<BVF,Op> OpPimpl;
  typedef OpBase  BOp;
  //  Teuchos::RCP<Pimpact::OperatorBase<Pimpact::MultiField<Field> > >
  //  typedef OpBase BOp;

	auto space = Pimpact::createSpace();

  auto b = Pimpact::createInitMVF<S,O>( Pimpact::Streaming2DFlow, space );
  auto x = Pimpact::createInitMVF<S,O>( Pimpact::Zero2DFlow, space );


  auto op = Pimpact::createMultiOperatorBase<BVF>( Pimpact::createDtLapOp(space) );

  auto para = Pimpact::createLinSolverParameter("GMRES",1.e-3);

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

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::MultiField<SF> BSF;

  typedef Pimpact::DivGradOp<S,O> Op;
  typedef Pimpact::MultiOpWrap<Op> MuOp;
  typedef Pimpact::OperatorBase<BSF> OpBase;
//  typedef Pimpact::OperatorPimpl<BSF,MuOp> OpPimpl;
  typedef OpBase  BOp;

	auto space = Pimpact::createSpace();

  auto temp = Pimpact::createVectorField<double,int>( space );
  temp->initField(Pimpact::PoiseuilleFlow2D_inX);
  //  temp->init(0.);

  auto p = Pimpact::createScalarField<S,O>(space);

  auto x = Pimpact::createMultiField<SF>(*p,1);
  auto b = x->clone();

  auto div = Teuchos::rcp( new Pimpact::DivOp<S,O>( space ) );

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

  auto para = Pimpact::createLinSolverParameter("GMRES",1.e-9);
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

