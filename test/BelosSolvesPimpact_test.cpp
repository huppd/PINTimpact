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
double errorTolSlack = 1e+1;

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

	init_impact(0,0);

  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::MultiField<VF> MVF;

  typedef Pimpact::Helmholtz<S,O> Op;
  typedef Pimpact::MultiOpWrap<Op> MOp;
  typedef Pimpact::OperatorBase<MVF> BOp;

  auto fS = Pimpact::createFieldSpace<O>();

  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();


  auto x = Pimpact::createMultiField( Pimpact::createVectorField<S,O>(fS, iIS, fIS) );
  auto b = Pimpact::createMultiField( Pimpact::createVectorField<S,O>(fS, iIS, fIS) );

  b->init( 1. );

  auto op = Pimpact::createOperatorBase<MVF,MOp>();

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
  TEST_EQUALITY( solver->achievedTol()<1.e-3, true );

  x->write(111);

}



TEUCHOS_UNIT_TEST( BelosSolver, DtLapOp ) {

  typedef int O;
  typedef double S;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField<MVF> BVF;

  typedef Pimpact::DtLapOp<S,O> Op;
  typedef Pimpact::MultiOpWrap<Op> OpMV;
  typedef Pimpact::OperatorBase<BVF> OpBase;
//  typedef Pimpact::OperatorPimpldep<BVF,Op> OpPimpl;
  typedef OpBase  BOp;
//  Teuchos::RCP<Pimpact::OperatorBase<Pimpact::MultiField<Field> > >
//  typedef OpBase BOp;

  auto fS = Pimpact::createFieldSpace<O>();

  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto b = Pimpact::createInitMVF<S,O>( Pimpact::Streaming2DFlow, fS, iIS, fIS );
  auto x = Pimpact::createInitMVF<S,O>( Pimpact::Zero2DFlow, fS, iIS, fIS );

//  auto A = Pimpact::createOperatorMV<Pimpact::Helmholtz<double,int> >();

  auto op = Pimpact::createOperatorBase<BVF,OpMV>();

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
  TEST_EQUALITY( solver->achievedTol()<1.e-3, true );

  x->write(20);

}



TEUCHOS_UNIT_TEST( BelosSolver, DivGrad ) {

  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::MultiField<SF> BSF;

  typedef Pimpact::DivGradOp<S,O> Op;
  typedef Pimpact::MultiOpWrap<Op> MuOp;
  typedef Pimpact::OperatorBase<BSF> OpBase;
  typedef Pimpact::OperatorPimpl<BSF,MuOp> OpPimpl;
  typedef OpBase  BOp;

  auto fS = Pimpact::createFieldSpace<O>();

  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto temp = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  temp->initField(Pimpact::Poiseuille2D_inX);
//  temp->init(0.);



  auto p = Pimpact::createScalarField<S,O>(fS);

  auto x = Pimpact::createMultiField<SF>(*p,1);
  auto b = x->clone();

  auto div = Teuchos::rcp( new Pimpact::Div<S,O>() );
//  div->apply(*temp,b->getField(0) );
//  b->scale(-1.);
//  b->init( 1. );
//  b->write(9999);

//  x->init(0.);
//  b->init(0.);
//  x->random();
//  x->scale(1000.);
//  b->random();

  temp->initField(Pimpact::ZeroProf);
  auto op = Pimpact::createOperatorBase<BSF,MuOp>(
      Pimpact::createMultiOpWrap(
          Pimpact::createDivGradOp(temp) ) );

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
   TEST_EQUALITY( solver->achievedTol()<1.e-9, true );

   x->write(999);
   temp->write(999);

}

} // end of namespace

