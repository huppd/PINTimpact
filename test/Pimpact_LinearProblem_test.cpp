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
#include "Pimpact_FieldFactory.hpp"
#include "Pimpact_OperatorMV.hpp"
#include "Pimpact_Operator.hpp"
#include "BelosPimpactAdapter.hpp"
#include "Pimpact_LinearProblem.hpp"

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
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing
	typedef double Scalar;
	typedef Pimpact::MultiField<Pimpact::VectorField<double,int> > MV;
	typedef Pimpact::OperatorMV< Pimpact::Helmholtz<double,int> >  OP;

	init_impact(0,0);
	auto fS = Pimpact::createFieldSpace<int>();

	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto X = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,1);
	auto B = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,1);

	X->Init(0.);
	B->Random();
	B->Scale(100.);

	auto A = Pimpact::createOperatorMV<Pimpact::Helmholtz<double,int> >();

	// Make an empty new parameter list.
	RCP<ParameterList> solverParams = parameterList();

	// Set some Belos parameters.
	//
	// "Num Blocks" = Maximum number of Krylov vectors to store.  This
	// is also the restart length.  "Block" here refers to the ability
	// of this particular solver (and many other Belos solvers) to solve
	// multiple linear systems at a time, even though we are only solving
	// one linear system in this example.
	solverParams->set ("Num Blocks", 40);
	solverParams->set ("Maximum Iterations", 400);
	solverParams->set ("Convergence Tolerance", 1.0e-1);
	solverParams->set ("Output Frequency", 50);
	solverParams->set ("Output Style", 1);
	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
	Belos::TimingDetails + Belos::StatusTestDetails);

	// Create the Pimpact::LinearSolver solver.
	Belos::SolverFactory<Scalar, MV, OP > factory;

	RCP<Belos::SolverManager<Scalar, MV, OP > > solver =
				 factory.create ("CG", solverParams);

	RCP<Belos::LinearProblem<Scalar, MV, OP > > problem =
			 rcp (new Belos::LinearProblem<Scalar, MV, OP > (A, X, B));

	problem->setProblem(X,B);

	solver->setProblem (problem);


	auto H_prob = Teuchos::rcp( new Pimpact::LinearProblem<Scalar,MV,OP>(solver,problem) );

//	std::cout << H_prob->getSolver()->getCurrentParameters();

	Belos::ReturnType result = H_prob->solve(X,B);
	TEST_EQUALITY( result,Belos::Converged);

	X->write(200);

}


TEUCHOS_UNIT_TEST( BelosSolver, HelmholtzMV2 ) {
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing
	typedef double Scalar;
	typedef Pimpact::MultiField<Pimpact::VectorField<double,int> > MV;
	typedef Pimpact::OperatorMV< Pimpact::Helmholtz<double,int> >  OP;

	auto fS = Pimpact::createFieldSpace<int>();

	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto X = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,1);
	auto B = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,1);

	X->Init(0.);
	B->Random();
	B->Scale(100.);

	auto A = Pimpact::createOperatorMV<Pimpact::Helmholtz<double,int> >();

	// Make an empty new parameter list.
	RCP<ParameterList> solverParams = parameterList();

	// Set some Belos parameters.
	//
	// "Num Blocks" = Maximum number of Krylov vectors to store.  This
	// is also the restart length.  "Block" here refers to the ability
	// of this particular solver (and many other Belos solvers) to solve
	// multiple linear systems at a time, even though we are only solving
	// one linear system in this example.
	solverParams->set ("Num Blocks", 40);
	solverParams->set ("Maximum Iterations", 400);
	solverParams->set ("Convergence Tolerance", 1.0e-1);
	solverParams->set ("Output Frequency", 50);
	solverParams->set ("Output Style", 1);
	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
	Belos::TimingDetails + Belos::StatusTestDetails);

// Create the Pimpact::LinearSolver solver.
	auto H_prob = Pimpact::createLinearProblem<Scalar,MV,OP>( A, X, B, solverParams,"CG" );

	Belos::ReturnType result = H_prob->solve(X,B);
	TEST_EQUALITY( result,Belos::Converged);

	X->write(200);
}


TEUCHOS_UNIT_TEST( BelosSolver, Dt1L0 ) {
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing
	typedef double Scalar;
	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MV;
	typedef Pimpact::OperatorMV< Pimpact::DtL<double,int> >  OP;

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto X = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> >,double,int>( *vel,1 );
	auto B = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> >,double,int>( *vel,1 );

	X->Init(0.);

	X->Random();
	B->Random();
//	B->Scale(100.);

	auto A = Pimpact::createDtL<double,int>( 1., 0., 0. );

	// Make an empty new parameter list.
	RCP<ParameterList> solverParams = parameterList();

	// Set some Belos parameters.
	//
	// "Num Blocks" = Maximum number of Krylov vectors to store.  This
	// is also the restart length.  "Block" here refers to the ability
	// of this particular solver (and many other Belos solvers) to solve
	// multiple linear systems at a time, even though we are only solving
	// one linear system in this example.
	solverParams->set ("Num Blocks", 40);
	solverParams->set ("Maximum Iterations", 400);
	solverParams->set ("Convergence Tolerance", 1.0e-1);
	solverParams->set ("Output Frequency", 50);
	solverParams->set ("Output Style", 1);
	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
	Belos::TimingDetails + Belos::StatusTestDetails);

// Create the Pimpact::LinearSolver solver.
	auto H_prob = Pimpact::createLinearProblem<Scalar,MV,OP>( A, X, B, solverParams,"GMRES" );

	Belos::ReturnType result = H_prob->solve(X,B);
	TEST_EQUALITY( result,Belos::Converged);

	X->write(200);
}


TEUCHOS_UNIT_TEST( BelosSolver, Dt0L1 ) {
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing
	typedef double Scalar;
	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MV;
	typedef Pimpact::OperatorMV< Pimpact::DtL<double,int> >  OP;

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto X = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> >,double,int>( *vel,1 );
	auto B = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> >,double,int>( *vel,1 );

	X->Init(0.);
	B->Random();
	B->Scale(100.);

	auto A = Pimpact::createDtL<double,int>(0.,0.,10.);

	// Make an empty new parameter list.
	RCP<ParameterList> solverParams = parameterList();

	// Set some Belos parameters.
	//
	// "Num Blocks" = Maximum number of Krylov vectors to store.  This
	// is also the restart length.  "Block" here refers to the ability
	// of this particular solver (and many other Belos solvers) to solve
	// multiple linear systems at a time, even though we are only solving
	// one linear system in this example.
	solverParams->set ("Num Blocks", 40);
	solverParams->set ("Maximum Iterations", 400);
	solverParams->set ("Convergence Tolerance", 1.0e-1);
	solverParams->set ("Output Frequency", 50);
	solverParams->set ("Output Style", 1);
	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
	Belos::TimingDetails + Belos::StatusTestDetails);

// Create the Pimpact::LinearSolver solver.
	auto H_prob = Pimpact::createLinearProblem<Scalar,MV,OP>( A, X, B, solverParams,"GMRES" );

	Belos::ReturnType result = H_prob->solve(X,B);
	TEST_EQUALITY( result,Belos::Converged);

	X->write(200);
}


TEUCHOS_UNIT_TEST( BelosSolver, Schur ) {
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing
	typedef double Scalar;
	typedef Pimpact::MultiField<Pimpact::VectorField<double,int> > MV;
	typedef Pimpact::MultiField<Pimpact::ScalarField<double,int> > MVp;
	typedef Pimpact::OperatorMV< Pimpact::Helmholtz<double,int> >  OP;
	typedef Pimpact::OperatorMV< Pimpact::Div_Hinv_Grad<double,int> >  OPp;

	auto fS = Pimpact::createFieldSpace<int>();

	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto p = Pimpact::createScalarField<double,int>(fS);

	auto temp = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,1);

	auto Xp = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,1);
	auto Bp = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,1);
	auto X = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,1);
	auto B = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,1);

	X->Init(0.);
	B->Random();
	B->Scale(100.);

	Xp->Init(0.);
	Bp->Random();
	Bp->Scale(100.);

	// init operators
	auto lap = Pimpact::createHelmholtz<double,int>( 0.,1.);
	auto div = Teuchos::rcp( new Pimpact::Div<double,int>() );
	auto grad = Teuchos::rcp( new Pimpact::Grad<double,int>() );


	// Make an empty new parameter list.
	RCP<ParameterList> solverParams = parameterList();

	// Set some Belos parameters.
	//
	// "Num Blocks" = Maximum number of Krylov vectors to store.  This
	// is also the restart length.  "Block" here refers to the ability
	// of this particular solver (and many other Belos solvers) to solve
	// multiple linear systems at a time, even though we are only solving
	// one linear system in this example.
	solverParams->set ("Num Blocks", 40);
	solverParams->set ("Maximum Iterations", 400);
	solverParams->set ("Convergence Tolerance", 1.0e-1);
	solverParams->set ("Output Frequency", 50);
	solverParams->set ("Output Style", 1);

// Create the Pimpact::LinearSolver solver.
	auto H_prob = Pimpact::createLinearProblem<Scalar,MV,OP>( lap, X, B, solverParams,"CG" );

	auto schur = Pimpact::createDivHinvGrad<double,int>( X,
//			div, grad,
			H_prob );


	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
	Belos::TimingDetails + Belos::StatusTestDetails);

	auto schur_prob = Pimpact::createLinearProblem<Scalar,MVp,OPp>( schur, Xp,Bp, solverParams, "GMRES");
	Belos::ReturnType result = schur_prob->solve(Xp,Bp);
	TEST_EQUALITY( result,Belos::Converged);

	X->write(200);

}

TEUCHOS_UNIT_TEST( BelosSolver, Div_DtLinv_Grad ) {
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing

	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MVF;
	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > > MSF;
	typedef Pimpact::OperatorMV< Pimpact::DtL<double,int> >  OP;
	typedef Pimpact::OperatorMV< Pimpact::Div_DtLinv_Grad<double,int> >  OP2;

	auto temp = Pimpact::createMultiModeVectorField<double,int>();

	auto X = Pimpact::createMultiModeScalarField<double,int>();
	auto B = Pimpact::createMultiModeScalarField<double,int>();

	X->Init(0.);
	B->Random();

	auto A = Pimpact::createDtL<double,int>(0.,0.,10.);

	// Make an empty new parameter list.
	RCP<ParameterList> solverParams = parameterList();

	// Set some Belos parameters.
	//
	// "Num Blocks" = Maximum number of Krylov vectors to store.  This
	// is also the restart length.  "Block" here refers to the ability
	// of this particular solver (and many other Belos solvers) to solve
	// multiple linear systems at a time, even though we are only solving
	// one linear system in this example.
	solverParams->set ("Num Blocks", 40);
	solverParams->set ("Maximum Iterations", 400);
	solverParams->set ("Convergence Tolerance", 1.0e-1);
	solverParams->set ("Output Frequency", 50);
	solverParams->set ("Output Style", 1);

// Create the Pimpact::LinearSolver solver.
	auto H_prob = Pimpact::createLinearProblem<double,MVF,OP>( A, temp, temp, solverParams,"GMRES" );

	auto schur = Pimpact::createDivDtLinvGrad<double,int>( temp, H_prob );

	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
	Belos::TimingDetails + Belos::StatusTestDetails);

	auto schur_prob = Pimpact::createLinearProblem<double,MSF,OP2>( schur, X, B, solverParams, "GMRES" );

	schur_prob->solve(X,B);

}

} // namespace

