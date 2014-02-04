// Pimpact_SalarVectorSpace_test.cpp

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
#include "Pimpact_MultiField.hpp"
#include "Pimpact_OperatorMV.hpp"
#include "Pimpact_Operator.hpp"
#include "BelosPimpactAdapter.hpp"

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

	auto X = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,1);
	auto B = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,1);

	X->Init(0.);
	B->Init(1.);

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

				auto fS = Pimpact::createFieldSpace<int>();

				auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
				auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

				auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

//				vel->init(0.);
//				vel->init_field();

				auto X = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>( *vel, 1, Pimpact::DeepCopy );
				auto B = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>( *vel, 1 );

//			 X->Init(0.);
			 X->GetVec(0).init_field();
			 B->Init(1.);

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


//	TEUCHOS_UNIT_TEST( BelosSolver, Lap ) {
//
//				auto fS = Pimpact::createFieldSpace<int>();
//
////				auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
////				auto fIS = Pimpact::createFullFieldIndexSpaces<int>();
//
//				auto p = Pimpact::createScalarField<double,int>(fS);
//
//				auto X = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,1);
//				auto B = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,1);
//
//			 X->Init(0.);
////			 X->Random();
//			 B->Init(1.);
//
//				auto A = Pimpact::createOperatorMV<Pimpact::Div_Grad<double,int> >();
//
//			 using Teuchos::ParameterList;
//			 using Teuchos::parameterList;
//			 using Teuchos::RCP;
//			 using Teuchos::rcp; // Save some typing
//			 typedef double Scalar;
//			 typedef Pimpact::MultiField<Pimpact::ScalarField<double,int> > MV;
//			 typedef Pimpact::OperatorMV< Pimpact::Div_Grad<double,int> >  OP;
//
//				 // The ellipses represent the code you would normally use to create
//				 // the sparse matrix, preconditioner, right-hand side, and initial
//				 // guess for the linear system AX=B you want to solve.
//	//			 RCP<OP> A = ...; // The sparse matrix / operator A
//	////			 RCP<OP> M = ...; // The (right) preconditioner M
//	//			 RCP<MV> B = ...; // Right-hand side of AX=B
//	//			 RCP<MV> X = ...; // Initial guess for the solution
//
//				 Belos::SolverFactory<Scalar, MV, OP > factory;
//				 // Make an empty new parameter list.
//				 RCP<ParameterList> solverParams = parameterList();
//
//				 // Set some GMRES parameters.
//				 //
//				 // "Num Blocks" = Maximum number of Krylov vectors to store.  This
//				 // is also the restart length.  "Block" here refers to the ability
//				 // of this particular solver (and many other Belos solvers) to solve
//				 // multiple linear systems at a time, even though we are only solving
//				 // one linear system in this example.
//				 solverParams->set ("Num Blocks", 100);
//				 solverParams->set ("Maximum Iterations", 1800);
//				 solverParams->set ("Convergence Tolerance", 1.0e-8);
//				 solverParams->set ("Output Frequency", 50);
//				 solverParams->set ("Output Style", 1);
//				 solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
//			       Belos::TimingDetails + Belos::StatusTestDetails);
//
//				 // Create the GMRES solver.
//				 RCP<Belos::SolverManager<Scalar, MV, OP > > solver =
//				   factory.create ("GMRES", solverParams);
////				   factory.create ("CG", solverParams);
//
//				 // Create a LinearProblem struct with the problem to solve.
//				 // A, X, B, and M are passed by (smart) pointer, not copied.
//				 RCP<Belos::LinearProblem<Scalar, MV, OP > > problem =
//				   rcp (new Belos::LinearProblem<Scalar, MV, OP > (A, X, B));
//	//			 problem->setRightPrec (M);
//
//				 std::cout << "param\n" << *solver->getValidParameters();
//				 problem->setProblem(X,B);
//	//			 std::cout << "problem set? " << problem->isProblemSet() << "\n";
//
//				 // Tell the solver what problem you want to solve.
//				 solver->setProblem (problem);
//	//			 TEST_EQUALITY( solver->getProblem(), *problem)
//	//			 std::cout << *problem;
//
//				 // Attempt to solve the linear system.  result == Belos::Converged
//				 // means that it was solved to the desired tolerance.  This call
//				 // overwrites X with the computed approximate solution.
//				 Belos::ReturnType result = solver->solve();
//				 TEST_EQUALITY( result,Belos::Converged);
//
//				 X->write(3);
//
//				 // Ask the solver how many iterations the last solve() took.
//				 const int numIters = solver->getNumIters();
//
//	//				TEST_EQUALITY( res, true );
//		}
////		const int n = mv->GetVecLength();
////		std::vector<double> normval(m);
//
//	  // test different float values, assures that initial and norm work smoothly
//		for( double i=0.; i< 200.1; ++i ) {
//			mv->Init(i/2.);
//			mv->Norm(normval,Belos::TwoNorm);
////			mv->Assign(*mv);
////			mv->Dot(*mv,normval);
//			for( int j=0; j<m; ++j )
//				TEST_EQUALITY( std::pow(i/2.,2)*n, normval[j] );
//		}
//	}


} // namespace

