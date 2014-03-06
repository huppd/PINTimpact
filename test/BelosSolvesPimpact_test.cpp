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
	auto fS = Pimpact::createFieldSpace<int>();

	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto X = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,1);
	auto B = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,1);

	X->init(0.);
	B->init(1.);

	auto A = Pimpact::createOperatorMV<Pimpact::Helmholtz<double,int> >();

	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing
	typedef double Scalar;
	typedef Pimpact::MultiField<Pimpact::VectorField<double,int> > MV;
	typedef Pimpact::OperatorMV< Pimpact::Helmholtz<double,int> >  OP;

	// The ellipses represent the code you would normally use to create
	// the sparse matrix, preconditioner, right-hand side, and initial
	// guess for the linear system AX=B you want to solve.
	//			 RCP<OP> A = ...; // The sparse matrix / operator A
	////			 RCP<OP> M = ...; // The (right) preconditioner M
	//			 RCP<MV> B = ...; // Right-hand side of AX=B
	//			 RCP<MV> X = ...; // Initial guess for the solution

	Belos::SolverFactory<Scalar, MV, OP > factory;
	// Make an empty new parameter list.
	RCP<ParameterList> solverParams = parameterList();

	 // Set some GMRES parameters.
	 //
	 // "Num Blocks" = Maximum number of Krylov vectors to store.  This
	 // is also the restart length.  "Block" here refers to the ability
	 // of this particular solver (and many other Belos solvers) to solve
	 // multiple linear systems at a time, even though we are only solving
	 // one linear system in this example.
	 solverParams->set ("Num Blocks", 40);
	 solverParams->set ("Maximum Iterations", 400);
	 solverParams->set ("Convergence Tolerance", 1.0e-2);
	 solverParams->set ("Output Frequency", 50);
	 solverParams->set ("Output Style", 1);
	 solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
			 Belos::TimingDetails + Belos::StatusTestDetails);

	 // Create the GMRES solver.
	 RCP<Belos::SolverManager<Scalar, MV, OP > > solver =
		 factory.create ("CG", solverParams);

	 // Create a LinearProblem struct with the problem to solve.
	 // A, X, B, and M are passed by (smart) pointer, not copied.
	 RCP<Belos::LinearProblem<Scalar, MV, OP > > problem =
		 rcp (new Belos::LinearProblem<Scalar, MV, OP > (A, X, B));
//			 problem->setRightPrec (M);

	 std::cout << "param\n" << *solver->getValidParameters();
	 problem->setProblem(X,B);
//			 std::cout << "problem set? " << problem->isProblemSet() << "\n";

	 // Tell the solver what problem you want to solve.
	 solver->setProblem (problem);
//			 TEST_EQUALITY( solver->getProblem(), *problem)
//			 std::cout << *problem;

	 // Attempt to solve the linear system.  result == Belos::Converged
	 // means that it was solved to the desired tolerance.  This call
	 // overwrites X with the computed approximate solution.
	 Belos::ReturnType result = solver->solve();
	 TEST_EQUALITY( result,Belos::Converged);

	 X->write(2);

	 // Ask the solver how many iterations the last solve() took.
//			 const int numIters = solver->getNumIters();

//				TEST_EQUALITY( res, true );
}


TEUCHOS_UNIT_TEST( BelosSolver, HelmholtzMV2 ) {
  typedef int O;

  auto fS = Pimpact::createFieldSpace<O>();

  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

//				vel->init(0.);
//				vel->init_field();

	auto X = Pimpact::createMultiField<Pimpact::VectorField<double,int> >( *vel, 1, Pimpact::DeepCopy );
	auto B = Pimpact::createMultiField<Pimpact::VectorField<double,int> >( *vel, 1 );

//			 X->Init(0.);
	X->getField(0).initField();
	B->init(1.);

	auto A = Pimpact::createOperatorMV<Pimpact::Helmholtz<double,int> >();

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing
	typedef double Scalar;
	typedef Pimpact::MultiField<Pimpact::VectorField<double,int> > MV;
	typedef Pimpact::OperatorMV< Pimpact::Helmholtz<double,int> >  OP;

				 // The ellipses represent the code you would normally use to create
				 // the sparse matrix, preconditioner, right-hand side, and initial
				 // guess for the linear system AX=B you want to solve.
	//			 RCP<OP> A = ...; // The sparse matrix / operator A
	////			 RCP<OP> M = ...; // The (right) preconditioner M
	//			 RCP<MV> B = ...; // Right-hand side of AX=B
	//			 RCP<MV> X = ...; // Initial guess for the solution

				 Belos::SolverFactory<Scalar, MV, OP > factory;
				 // Make an empty new parameter list.
				 RCP<ParameterList> solverParams = parameterList();

				 // Set some GMRES parameters.
				 //
				 // "Num Blocks" = Maximum number of Krylov vectors to store.  This
				 // is also the restart length.  "Block" here refers to the ability
				 // of this particular solver (and many other Belos solvers) to solve
				 // multiple linear systems at a time, even though we are only solving
				 // one linear system in this example.
				 solverParams->set ("Num Blocks", 40);
				 solverParams->set ("Maximum Iterations", 400);
				 solverParams->set ("Convergence Tolerance", 1.0e-2);
				 solverParams->set ("Output Frequency", 50);
				 solverParams->set ("Output Style", 1);
				 solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
			       Belos::TimingDetails + Belos::StatusTestDetails);

				 // Create the GMRES solver.
				 RCP<Belos::SolverManager<Scalar, MV, OP > > solver =
				   factory.create ("CG", solverParams);

				 // Create a LinearProblem struct with the problem to solve.
				 // A, X, B, and M are passed by (smart) pointer, not copied.
				 RCP<Belos::LinearProblem<Scalar, MV, OP > > problem =
				   rcp (new Belos::LinearProblem<Scalar, MV, OP > (A, X, B));
	//			 problem->setRightPrec (M);

				 std::cout << "param\n" << *solver->getValidParameters();
				 problem->setProblem(X,B);
	//			 std::cout << "problem set? " << problem->isProblemSet() << "\n";

				 // Tell the solver what problem you want to solve.
				 solver->setProblem (problem);
	//			 TEST_EQUALITY( solver->getProblem(), *problem)
	//			 std::cout << *problem;

				 // Attempt to solve the linear system.  result == Belos::Converged
				 // means that it was solved to the desired tolerance.  This call
				 // overwrites X with the computed approximate solution.
				 Belos::ReturnType result = solver->solve();
				 TEST_EQUALITY( result,Belos::Converged);

				 X->write(20);

				 // Ask the solver how many iterations the last solve() took.
//				 const int numIters = solver->getNumIters();

	//				TEST_EQUALITY( res, true );
}


TEUCHOS_UNIT_TEST( BelosSolver, DtL ) {

  typedef int O;
  typedef double S;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField<MVF> BVF;

  typedef Pimpact::DtL<S,O> Op;
  typedef Pimpact::OperatorMV<Op> OpMV;
  typedef Pimpact::OperatorBase<BVF> OpBase;
  typedef Pimpact::OperatorPimpl<BVF,Op> OpPimpl;
  typedef OpBase  BOp;
//  Teuchos::RCP<Pimpact::OperatorBase<Pimpact::MultiField<Field> > >
//  typedef OpBase BOp;

  auto fS = Pimpact::createFieldSpace<O>();

  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto b = Pimpact::createInitMVF<S,O>( Pimpact::Streaming2DFlow, fS, iIS, fIS );
  auto x = Pimpact::createInitMVF<S,O>( Pimpact::Zero2DFlow, fS, iIS, fIS );

//  auto A = Pimpact::createOperatorMV<Pimpact::Helmholtz<double,int> >();

  auto op = Pimpact::createOperatorBase<BVF,Op>();

  auto para = Pimpact::createLinSolverParameter("GMRES",1.e-3);

  Belos::SolverFactory<S, BVF, BOp > factory;
  // Make an empty new parameter list.

  // Create the GMRES solver.
  Teuchos::RCP<Belos::SolverManager<S, BVF, BOp > > solver =
           factory.create( "GMRES", para->get() );

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  Teuchos::RCP<Belos::LinearProblem<S, BVF, BOp > > problem =
           Teuchos::rcp (new Belos::LinearProblem<S, BVF, BOp > (op, x, b));

  std::cout << "param\n" << *solver->getValidParameters();
  problem->setProblem(x,b);

  // Tell the solver what problem you want to solve.
  solver->setProblem( problem );
  //       TEST_EQUALITY( solver->getProblem(), *problem)
  //       std::cout << *problem;

  // Attempt to solve the linear system.  result == Belos::Converged
  // means that it was solved to the desired tolerance.  This call
  // overwrites X with the computed approximate solution.
  Belos::ReturnType result = solver->solve();
  TEST_EQUALITY( result,Belos::Converged);

  x->write(20);

         // Ask the solver how many iterations the last solve() took.
//         const int numIters = solver->getNumIters();
}

TEUCHOS_UNIT_TEST( BelosSolver, DivGrad ) {
  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::MultiField<SF> BSF;

  typedef Pimpact::Div_Grad<S,O> Op;
  typedef Pimpact::OperatorMV<Op> OpMV;
  typedef Pimpact::OperatorBase<BSF> OpBase;
  typedef Pimpact::OperatorPimpl<BSF,Op> OpPimpl;
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
  div->apply(*temp,b->getField(0) );
//  b->scale(-1.);
  b->write(9999);

  x->init(0.);
//  b->init(0.);
  x->random();
  x->scale(1000.);
//  b->random();

  temp->initField(Pimpact::ZeroProf);
  auto op = Pimpact::createOperatorBase<BSF,Op>( Teuchos::rcp( new Op(temp) ) );

  auto para = Pimpact::createLinSolverParameter("GMRES",1.e-9)->get();
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

   // Attempt to solve the linear system.  result == Belos::Converged
   // means that it was solved to the desired tolerance.  This call
   // overwrites X with the computed approximate solution.
   Belos::ReturnType result = solver->solve();
   TEST_EQUALITY( result, Belos::Converged);

   x->write(999);
   temp->write(999);

          // Ask the solver how many iterations the last solve() took.
 //         const int numIters = solver->getNumIters();
//  auto lap = Pimpact::createOperatorMV<Op>();
//
//  Teuchos::RCP<Belos::OutputManager<S> > MyOM =
//      Teuchos::rcp( new Belos::OutputManager<S>(Belos::Errors + Belos::Warnings + Belos::IterationDetails +
//          Belos::OrthoDetails + Belos::FinalSummary + Belos::TimingDetails +
//          Belos::StatusTestDetails + Belos::Debug, rcp(&out,false)) );
//
////  mv->random();
////  lap->apply( *mv, *mv2 );
////  mv2->write(20);
//
//  bool res =// true;
//      Belos::TestOperatorTraits< S, BSF, BOp > (MyOM,mv,lap);
//
//  TEST_EQUALITY( res, true );
}

} // end of namespace

