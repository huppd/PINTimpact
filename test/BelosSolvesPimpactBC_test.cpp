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
#include "Pimpact_BoundaryConditions.hpp"

#include "BelosSolverFactory.hpp"

#include <iostream>

namespace {

	bool testMpi = true;
	double errorTolSlack = 1e+1;

	TEUCHOS_STATIC_SETUP()
	{
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

			auto bc = Pimpact::createPeriodicChannelBC2D();
			init_impact(0,0,bc);
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

				 // Create the solver.
				 RCP<Belos::SolverManager<Scalar, MV, OP > > solver =
				   factory.create ("CG", solverParams);
//				   factory.create ("GMRES", solverParams);

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

				 X->write(40);

				 // Ask the solver how many iterations the last solve() took.
//				 const int numIters = solver->getNumIters();

	//				TEST_EQUALITY( res, true );
		}

	TEUCHOS_UNIT_TEST( BelosSolver, HelmholtzMV3 ) {

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

					auto A = Pimpact::createHelmholtz<double,int>( 0.,1.);

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

					 // Create the solver.
					 RCP<Belos::SolverManager<Scalar, MV, OP > > solver =
					   factory.create ("CG", solverParams);
	//				   factory.create ("GMRES", solverParams);

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

					 X->write(60);

					 // Ask the solver how many iterations the last solve() took.
//					 const int numIters = solver->getNumIters();

		//				TEST_EQUALITY( res, true );
			}
//
//
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
//
//
//	TEUCHOS_UNIT_TEST( MultiField, clone ) {
//
//		auto sVS = Pimpact::createFieldSpace<int>();
//		auto sIS = Pimpact::createScalarIndexSpace<int>();
//		auto p = Pimpact::createScalarField<double,int>(sVS);
//
//		auto mv = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,1);
//
//		auto mv2 = mv->Clone(10);
//
//		int n1(mv->GetNumberVecs());
//		int n2(mv2->GetNumberVecs());
//
//		TEST_EQUALITY( 1, n1 );
//		TEST_EQUALITY( 10, n2 );
//
//	}
//
//	TEUCHOS_UNIT_TEST( MultiField, CloneCopy ) {
//
//		auto sVS = Pimpact::createFieldSpace<int>();
//		auto sIS = Pimpact::createScalarIndexSpace<int>();
//		auto p = Pimpact::createScalarField<double,int>(sVS);
//
//		auto mv = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//		mv->Random();
//		auto mv2 = mv->CloneCopy();
//
//		int n1(mv->GetNumberVecs());
//		int n2(mv2->GetNumberVecs());
//
//		TEST_EQUALITY( 10, n1 );
//		TEST_EQUALITY( n1, n2 );
//
//		std::vector<double> norm1(n1);
//		std::vector<double> norm2(n2);
//
//		mv->Norm(norm1);
//		mv2->Norm(norm2);
//		for( int i=0; i<n1; ++i)
//			TEST_EQUALITY( norm1[i], norm2[i] );
//
//	}
//
//	TEUCHOS_UNIT_TEST( MultiField, CloneCopy2 ) {
//
//		auto sVS = Pimpact::createFieldSpace<int>();
//		auto sIS = Pimpact::createScalarIndexSpace<int>();
//		auto p = Pimpact::createScalarField<double,int>(sVS);
//
//		auto mv = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//		mv->Random();
//
//		std::vector<int> index(5);
//		for(int i=0; i<5; ++i)
//			index[i] = 2*i;
//
//		auto mv2 = mv->CloneCopy(index);
//
//		unsigned int n1 = (mv->GetNumberVecs());
//		unsigned int n2 = (mv2->GetNumberVecs());
//
//		TEST_EQUALITY( 10, n1 );
//		TEST_EQUALITY( 5, n2 );
//		TEST_EQUALITY( index.size(), n2 );
//
//		std::vector<double> norm1(n1);
//		std::vector<double> norm2(n2);
//
//		mv->Norm(norm1);
//		mv2->Norm(norm2);
//
//		for( unsigned int i=0; i<index.size(); ++i)
//			TEST_EQUALITY( norm1[index[i]], norm2[i] );
//	}
//
//	TEUCHOS_UNIT_TEST( MultiField, CloneCopy3 ) {
//
//		auto sVS = Pimpact::createFieldSpace<int>();
//		auto sIS = Pimpact::createScalarIndexSpace<int>();
//		auto p = Pimpact::createScalarField<double,int>(sVS);
//
//		auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//		mv1->Random();
//
//		Teuchos::Range1D index(2,7);
//
//		auto mv2 = mv1->CloneCopy(index);
//
//		unsigned int n1(mv1->GetNumberVecs());
//		unsigned int n2(mv2->GetNumberVecs());
//
//		TEST_EQUALITY( 10, n1 );
//		TEST_EQUALITY( index.size(), n2 );
//
//		std::vector<double> norm1(n1);
//		std::vector<double> norm2(n2);
//
//		mv1->Norm(norm1);
//		mv2->Norm(norm2);
//		for( int i=0; i<index.size(); ++i)
//			TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//
//	}
//
//	TEUCHOS_UNIT_TEST( MultiField, CloneViewNonConst1 ) {
//
//		auto sVS = Pimpact::createFieldSpace<int>();
//		auto sIS = Pimpact::createScalarIndexSpace<int>();
//		auto p = Pimpact::createScalarField<double,int>(sVS);
//
//		auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//		mv1->Init(0.);
//
//		std::vector<int> index(5);
//		for(int i=0; i<5; ++i)
//			index[i] = 2*i;
//
//		auto mv2 = mv1->CloneViewNonConst(index);
//
//		unsigned int n1 = (mv1->GetNumberVecs());
//		unsigned int n2 = (mv2->GetNumberVecs());
////
//		TEST_EQUALITY( 10, n1 );
//		TEST_EQUALITY( 5, n2 );
//		TEST_EQUALITY( index.size(), n2 );
//
//		std::vector<double> norm1(n1);
//		std::vector<double> norm2(n2);
//
//		mv2->Random();
//
//		mv1->Norm(norm1);
//		mv2->Norm(norm2);
//
//		for( unsigned int i=0; i<index.size(); ++i)
//			TEST_EQUALITY( norm1[index[i]], norm2[i] );
//	}
//
//
//	TEUCHOS_UNIT_TEST( MultiField, CloneViewNonConst2 ) {
//
//			auto sVS = Pimpact::createFieldSpace<int>();
//			auto sIS = Pimpact::createScalarIndexSpace<int>();
//			auto p = Pimpact::createScalarField<double,int>(sVS);
//
//			auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//			mv1->Init(0.);
//
//			Teuchos::Range1D index(2,7);
//
//			auto mv2 = mv1->CloneViewNonConst(index);
//
//			unsigned int n1 = (mv1->GetNumberVecs());
//			unsigned int n2 = (mv2->GetNumberVecs());
//
//			TEST_EQUALITY( 10, n1 );
//			TEST_EQUALITY( index.size(), n2 );
//
//			std::vector<double> norm1(n1);
//			std::vector<double> norm2(n2);
//
//			mv2->Random();
//
//			mv1->Norm(norm1);
//			mv2->Norm(norm2);
//
//			for( unsigned int i=0; i<index.size(); ++i)
//				TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//		}
//
//	TEUCHOS_UNIT_TEST( MultiField, CloneView1 ) {
//
//			auto sVS = Pimpact::createFieldSpace<int>();
//			auto sIS = Pimpact::createScalarIndexSpace<int>();
//			auto p = Pimpact::createScalarField<double,int>(sVS);
//
//			auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//			mv1->Init(0.);
//
//			std::vector<int> index(5);
//			for(int i=0; i<5; ++i)
//				index[i] = 2*i;
//
//			auto mv2 = mv1->CloneView(index);
//
//			unsigned int n1 = (mv1->GetNumberVecs());
//			unsigned int n2 = (mv2->GetNumberVecs());
//	//
//			TEST_EQUALITY( 10, n1 );
//			TEST_EQUALITY( 5, n2 );
//			TEST_EQUALITY( index.size(), n2 );
//
//			std::vector<double> norm1(n1);
//			std::vector<double> norm2(n2);
//
////			mv2->Random(); //< this should give compile error
//
//			mv1->Random();
//
//			mv1->Norm(norm1);
//			mv2->Norm(norm2);
//
//			for( unsigned int i=0; i<index.size(); ++i)
//				TEST_EQUALITY( norm1[index[i]], norm2[i] );
//		}
//
//	TEUCHOS_UNIT_TEST( MultiField, CloneViewt2 ) {
//
//			auto sVS = Pimpact::createFieldSpace<int>();
//			auto sIS = Pimpact::createScalarIndexSpace<int>();
//			auto p = Pimpact::createScalarField<double,int>(sVS);
//
//			auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//			mv1->Init(0.);
//
//			Teuchos::Range1D index(2,7);
//
//			auto mv2 = mv1->CloneView(index);
//
//			unsigned int n1 = (mv1->GetNumberVecs());
//			unsigned int n2 = (mv2->GetNumberVecs());
//
//			TEST_EQUALITY( 10, n1 );
//			TEST_EQUALITY( index.size(), n2 );
//
//			std::vector<double> norm1(n1);
//			std::vector<double> norm2(n2);
//
////			mv2->Random(); // has to give compilation error
//			mv1->Random();
//
//			mv1->Norm(norm1);
//			mv2->Norm(norm2);
//
//			for( unsigned int i=0; i<index.size(); ++i)
//				TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//		}
//
//	TEUCHOS_UNIT_TEST( MultiField, TimesMatAdd ) {
//
//				auto sVS = Pimpact::createFieldSpace<int>();
//				auto sIS = Pimpact::createScalarIndexSpace<int>();
//				auto p = Pimpact::createScalarField<double,int>(sVS);
//
//				auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//				mv1->Init(0.);
//
//				Teuchos::Range1D index1(0,9);
//				std::vector<int> index2(10);
//				for(int i=0; i<10; ++i)
//					index2[i] = i;
//
//				auto mv2 = mv1->CloneView(index1);
//				auto mv3 = mv1->Clone(mv1->GetNumberVecs() );
//
//
//				unsigned int n1 = (mv1->GetNumberVecs());
//				unsigned int n2 = (mv2->GetNumberVecs());
//				unsigned int n3 = (mv3->GetNumberVecs());
//
//				TEST_EQUALITY( n1, n2 );
//				TEST_EQUALITY( n3, n2 );
//				TEST_EQUALITY( n3, n1 );
//
//
//	//			mv2->Random(); // has to give compilation error
//				mv1->Init(1.);
////				mv2->Random();
////				mv3->Assign(2.);
//
//
//				Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//				for( unsigned int j=0; j<n1; ++j)
//					for( unsigned int i=0; i<n1; ++i)
//						B(j,i) = 1./n1;
//
//				mv1->TimesMatAdd( 0.5, *mv2, B, 0.5 );
//
//
//				std::vector<double> norm1(n1);
//				std::vector<double> norm2(n2);
//
//				mv1->Norm(norm1);
//				mv2->Norm(norm2);
//
//				for( unsigned int i=0; i<n1; ++i) {
//					TEST_EQUALITY( norm1[i], norm2[i] );
//					TEST_EQUALITY( mv1->GetVecLength(), norm2[i] );
//				}
//
//				std::vector<double> scales(n1);
//				for( unsigned int j=0; j<n1; ++j){
//					scales[j] = (j+1);
////					scales[j] = 0;
//					for( unsigned int i=0; i<n1; ++i)
//						B(j,i) = 1./n1/(j+1);
////						B(j,i) = 1./5.;
//				}
////				std::cout << B;
////				scales[5] = 5.;
//				mv1->Init(1.);
//				mv3->Init(1.);
//
//				mv3->Scale(scales);
//
//				mv1->TimesMatAdd( 1., *mv3, B, 0. );
//
//
//				mv1->Norm(norm1);
//				mv2->Norm(norm2);
//
//				for( unsigned int i=0; i<n1; ++i) {
////					TEST_EQUALITY( norm1[i], norm2[i] );
//					TEST_FLOATING_EQUALITY( (double)mv1->GetVecLength(), norm2[i], 0.01 );
//				}
//			}
//
//	TEUCHOS_UNIT_TEST( MultiField, Add ) {
//
//					auto sVS = Pimpact::createFieldSpace<int>();
//					auto sIS = Pimpact::createScalarIndexSpace<int>();
//					auto p = Pimpact::createScalarField<double,int>(sVS);
//
//					auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//					mv1->Init(0.);
//
//					Teuchos::Range1D index1(0,9);
//					std::vector<int> index2(10);
//					for(int i=0; i<10; ++i)
//						index2[i] = i;
//
////					auto mv2 = mv1->CloneCopy(index1);
//					auto mv2 = mv1->CloneViewNonConst(index1);
//					auto mv3 = mv1->CloneView(index1);
//
//					unsigned int n1 = (mv1->GetNumberVecs());
//					unsigned int n2 = (mv2->GetNumberVecs());
//					unsigned int n3 = (mv3->GetNumberVecs());
//
//					TEST_EQUALITY( n1, n2 );
//					TEST_EQUALITY( n3, n2 );
//					TEST_EQUALITY( n3, n1 );
//
//
//		//			mv2->Random(); // has to give compilation error
//					mv1->Init(1.);
//					mv2->Init(1.);
//	//				mv2->Random();
//	//				mv3->Assign(2.);
//
//
//
//					mv1->Add( 0.5, *mv2, 0.5, *mv3);
//
//
//					std::vector<double> norm1(n1);
//					std::vector<double> norm2(n2);
//
//					mv1->Norm(norm1);
//					mv2->Norm(norm2);
//
//					for( unsigned int i=0; i<n1; ++i) {
////						TEST_EQUALITY( norm1[i], norm2[i] );
//						TEST_EQUALITY( mv1->GetVecLength(), norm2[i] );
//					}
//
//					mv1->Init(1.);
//					mv2->Init(1.);
//
//					mv2->Scale(0.5);
//
//					mv1->Add( 1., *mv2, 1., *mv3 );
//
//
//					mv1->Norm(norm1);
//					mv2->Assign(*mv1);
//					mv2->Norm(norm2);
//
//					for( unsigned int i=0; i<n1; ++i) {
//						TEST_EQUALITY( norm1[i], norm2[i] );
//						TEST_FLOATING_EQUALITY( (double)mv1->GetVecLength(), norm1[i], 0.01 );
//					}
//				}
//
//
//	TEUCHOS_UNIT_TEST( MultiField, Dot ) {
//
//					auto sVS = Pimpact::createFieldSpace<int>();
//					auto sIS = Pimpact::createScalarIndexSpace<int>();
//					auto p = Pimpact::createScalarField<double,int>(sVS);
//
//					auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//					mv1->Init(0.);
//
//					Teuchos::Range1D index1(0,9);
//					std::vector<int> index2(10);
//					for(int i=0; i<10; ++i)
//						index2[i] = i;
//
////					auto mv2 = mv1->CloneCopy(index1);
//					auto mv2 = mv1->CloneViewNonConst(index1);
//					auto mv3 = mv1->CloneView(index1);
//
//					unsigned int n1 = (mv1->GetNumberVecs());
//					unsigned int n2 = (mv2->GetNumberVecs());
//					unsigned int n3 = (mv3->GetNumberVecs());
//
//					TEST_EQUALITY( n1, n2 );
//					TEST_EQUALITY( n3, n2 );
//					TEST_EQUALITY( n3, n1 );
//
//
//		//			mv2->Random(); // has to give compilation error
//					mv1->Init(1.);
//					mv2->Init(1.);
//	//				mv2->Random();
//	//				mv3->Assign(2.);
//
//
//					std::vector<double> dots(n1);
//
//					mv1->Dot( *mv2, dots );
//
//
//					for( unsigned int i=0; i<n1; ++i) {
////						TEST_EQUALITY( norm1[i], norm2[i] );
//						TEST_EQUALITY( mv1->GetVecLength(), dots[i] );
//					}
//
//					mv2->Init(2.);
//
//					mv3->Dot( *mv2, dots  );
//					for( unsigned int i=0; i<n1; ++i) {
//						TEST_EQUALITY( 4*mv1->GetVecLength(), dots[i] );
//					}
//				}
//
//	TEUCHOS_UNIT_TEST( MultiField, Trans ) {
//
//				auto sVS = Pimpact::createFieldSpace<int>();
//				auto sIS = Pimpact::createScalarIndexSpace<int>();
//				auto p = Pimpact::createScalarField<double,int>(sVS);
//
//				auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);
//
//				mv1->Init(0.);
//
//				Teuchos::Range1D index1(0,9);
//				std::vector<int> index2(10);
//				for(int i=0; i<10; ++i)
//					index2[i] = i;
//
//				auto mv2 = mv1->CloneView(index1);
//				auto mv3 = mv1->CloneView(index1);
////				auto mv3 = mv1->CloneV(mv1->GetNumberVecs() );
//
//
//				unsigned int n1 = (mv1->GetNumberVecs());
//				unsigned int n2 = (mv2->GetNumberVecs());
//				unsigned int n3 = (mv3->GetNumberVecs());
//
//				TEST_EQUALITY( n1, n2 );
//				TEST_EQUALITY( n3, n2 );
//				TEST_EQUALITY( n3, n1 );
//
//
//	//			mv2->Random(); // has to give compilation error
//				mv1->Init(1.);
////				mv2->Random();
////				mv3->Assign(2.);
//
//
//				Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//
//				mv1->Trans( 1., *mv2, B );
//
//				for( unsigned int j=0; j<n1; ++j){
//					for( unsigned int i=0; i<n1; ++i)
//						TEST_EQUALITY( mv1->GetVecLength(), B(j,i) );
//				}
//
////				std::vector<double> norm1(n1);
////				std::vector<double> norm2(n2);
////
////				mv1->Norm(norm1);
////				mv2->Norm(norm2);
////
////				for( unsigned int i=0; i<n1; ++i) {
////					TEST_EQUALITY( norm1[i], norm2[i] );
////					TEST_EQUALITY( mv1->GetVecLength(), norm2[i] );
////				}
//
//				std::vector<double> scales(n1);
//				for( unsigned int i=0; i<scales.size(); ++i)
//					scales[i] = i*2;
//				mv1->Scale(scales);
//
//				mv2->Trans( 1., *mv3, B );
//
//				for( unsigned int j=0; j<n1; ++j){
//					for( unsigned int i=0; i<n1; ++i)
//						TEST_EQUALITY( scales[i]*scales[j]*mv1->GetVecLength(), B(j,i) );
//				}
//
////					scales[j] = (j+1);
//////					scales[j] = 0;
////					for( unsigned int i=0; i<n1; ++i)
////						B(j,i) = 1./n1/(j+1);
//////						B(j,i) = 1./5.;
////				}
//////				std::cout << B;
//////				scales[5] = 5.;
////				mv1->Init(1.);
////				mv3->Init(1.);
////
////				mv3->Scale(scales);
////
////				mv1->TimesMatAdd( 1., *mv3, B, 0. );
////
////
////				mv1->Norm(norm1);
////				mv2->Norm(norm2);
////
////				for( unsigned int i=0; i<n1; ++i) {
//////					TEST_EQUALITY( norm1[i], norm2[i] );
////					TEST_FLOATING_EQUALITY( (double)mv1->GetVecLength(), norm2[i], 0.01 );
////				}
//			}


} // namespace

