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
#include "Pimpact_OperatorFactory.hpp"
#include "BelosPimpactAdapter.hpp"
#include "Pimpact_LinearProblem.hpp"
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

	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing
	typedef double  S;
	typedef int     O;
	typedef Pimpact::MultiField< Pimpact::VectorField<S,O> > MF;
	typedef Pimpact::MultiOpWrap< Pimpact::HelmholtzOp<S,O> >  Op;


	auto space = Pimpact::createSpace();

	auto vel = Pimpact::createVectorField<S,O>(space);

	auto x = Pimpact::createMultiField<Pimpact::VectorField<S,O> >(*vel,1);
	auto b = Pimpact::createMultiField<Pimpact::VectorField<S,O> >(*vel,1);

	x->init(0.);
	b->init(1.);

	auto A = Pimpact::createMultiOperatorBase< MF, Pimpact::HelmholtzOp<S,O> >();


	RCP<ParameterList> param = Pimpact::createLinSolverParameter("GMRES",1.e-4);

	auto linprob = Pimpact::createLinearProblem<MF>(A,x,b,param,"GMRES");

	linprob->solve(x,b);
	x->write();

  TEST_EQUALITY( linprob->getSolver()->achievedTol()<1.e-4, true );

}

//
//TEUCHOS_UNIT_TEST( BelosSolver, HelmholtzMV2 ) {
//	using Teuchos::ParameterList;
//	using Teuchos::parameterList;
//	using Teuchos::RCP;
//	using Teuchos::rcp; // Save some typing
//	typedef double Scalar;
//	typedef Pimpact::MultiField<Pimpact::VectorField<double,int> > MV;
//	typedef Pimpact::OperatorMV< Pimpact::HelmholtzOp<double,int> >  OP;
//
//	auto fS = Pimpact::createFieldSpace<int>();
//
//	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
//	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();
//
//	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
//
//	auto X = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,1);
//	auto B = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,1);
//
//	X->init(0.);
//	B->init(1.);
//
//	auto A = Pimpact::createOperatorMV<Pimpact::HelmholtzOp<double,int> >();
//
//	// Make an empty new parameter list.
//	RCP<ParameterList> solverParams = parameterList();
//
//	// Set some Belos parameters.
//	//
//	// "Num Blocks" = Maximum number of Krylov vectors to store.  This
//	// is also the restart length.  "Block" here refers to the ability
//	// of this particular solver (and many other Belos solvers) to solve
//	// multiple linear systems at a time, even though we are only solving
//	// one linear system in this example.
//	solverParams->set ("Num Blocks", 40);
//	solverParams->set ("Maximum Iterations", 40);
//	solverParams->set ("Convergence Tolerance", 1.0e-1);
//	solverParams->set ("Output Frequency", 50);
//	solverParams->set ("Output Style", 1);
//	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
//	Belos::TimingDetails + Belos::StatusTestDetails);
//
//// Create the Pimpact::LinearSolver solver.
//	auto H_prob = Pimpact::createLinearProblem<Scalar,MV,OP>( A, X, B, solverParams,"CG" );
//
//	Belos::ReturnType result = H_prob->solve(X,B);
//	TEST_EQUALITY( result,Belos::Converged);
//
//	H_prob->apply( X, B );
//
//	X->write(200);
//}
//
//
//TEUCHOS_UNIT_TEST( BelosSolver, Dt1L0 ) {
//	using Teuchos::ParameterList;
//	using Teuchos::parameterList;
//	using Teuchos::RCP;
//	using Teuchos::rcp; // Save some typing
//	typedef double Scalar;
//	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MV;
//	typedef Pimpact::OperatorMV< Pimpact::DtL<double,int> >  OP;
//
//	auto fS = Pimpact::createFieldSpace<int>();
//	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
//	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();
//
//	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
//	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
//
//	auto vel = Pimpact::createModeField( velc, vels );
//
//	auto X = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >( *vel,1 );
//	auto B = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >( *vel,1 );
//
//
//	X->init(0.);
//	B->init(1.);
//
//
//	auto A = Pimpact::createOperatorMV( Pimpact::createDtL<double,int>( 1., 0., 0. ) );
//
//	// Make an empty new parameter list.
//	RCP<ParameterList> solverParams = parameterList();
//
//	// Set some Belos parameters.
//	//
//	// "Num Blocks" = Maximum number of Krylov vectors to store.  This
//	// is also the restart length.  "Block" here refers to the ability
//	// of this particular solver (and many other Belos solvers) to solve
//	// multiple linear systems at a time, even though we are only solving
//	// one linear system in this example.
//	solverParams->set ("Num Blocks", 40);
//	solverParams->set ("Maximum Iterations", 40);
//	solverParams->set ("Convergence Tolerance", 1.0e-1);
//	solverParams->set ("Output Frequency", 50);
//	solverParams->set ("Output Style", 1);
//	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
//	Belos::TimingDetails + Belos::StatusTestDetails);
//
//// Create the Pimpact::LinearSolver solver.
//	auto H_prob = Pimpact::createLinearProblem<Scalar,MV,OP>( A, X, B, solverParams,"GMRES" );
//
//	Belos::ReturnType result = H_prob->solve(X,B);
//
//	TEST_EQUALITY( result,Belos::Converged);
//
//	H_prob->apply( X, B );
//
//	X->write(200);
//}
//
//
//TEUCHOS_UNIT_TEST( BelosSolver, Dt0L1 ) {
//	using Teuchos::ParameterList;
//	using Teuchos::parameterList;
//	using Teuchos::RCP;
//	using Teuchos::rcp; // Save some typing
//	typedef double Scalar;
//	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MV;
//	typedef Pimpact::OperatorMV< Pimpact::DtL<double,int> >  OP;
//
//	auto fS = Pimpact::createFieldSpace<int>();
//	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
//	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();
//
//	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
//	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
//
//	auto vel = Pimpact::createModeField( velc, vels );
//
//	auto X = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >( *vel,1 );
//	auto B = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >( *vel,1 );
//
//	X->init(0.);
//	B->init(1.);
//
//	auto A = Pimpact::createDtL<double,int>(0.,0.,10.);
//
//	// Make an empty new parameter list.
//	RCP<ParameterList> solverParams = parameterList();
//
//	// Set some Belos parameters.
//	//
//	// "Num Blocks" = Maximum number of Krylov vectors to store.  This
//	// is also the restart length.  "Block" here refers to the ability
//	// of this particular solver (and many other Belos solvers) to solve
//	// multiple linear systems at a time, even though we are only solving
//	// one linear system in this example.
//	solverParams->set ("Num Blocks", 40);
//	solverParams->set ("Maximum Iterations", 40);
//	solverParams->set ("Convergence Tolerance", 1.0e-1);
//	solverParams->set ("Output Frequency", 50);
//	solverParams->set ("Output Style", 1);
//	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
//	Belos::TimingDetails + Belos::StatusTestDetails);
//
//// Create the Pimpact::LinearSolver solver.
//	auto H_prob = Pimpact::createLinearProblem<Scalar,MV,OP>( A, X, B, solverParams,"GMRES" );
//
//	Belos::ReturnType result = H_prob->solve(X,B);
//	TEST_EQUALITY( result,Belos::Converged);
//
//	H_prob->apply( X, B );
//
//	X->write(200);
//}
//
//
//TEUCHOS_UNIT_TEST( BelosSolver, Schur ) {
//	using Teuchos::ParameterList;
//	using Teuchos::parameterList;
//	using Teuchos::RCP;
//	using Teuchos::rcp; // Save some typing
//	typedef double Scalar;
//	typedef Pimpact::MultiField<Pimpact::VectorField<double,int> > MV;
//	typedef Pimpact::MultiField<Pimpact::ScalarField<double,int> > MVp;
//	typedef Pimpact::OperatorMV< Pimpact::HelmholtzOp<double,int> >  OP;
//	typedef Pimpact::OperatorMV< Pimpact::Div_Hinv_Grad<double,int> >  OPp;
//
//	auto fS = Pimpact::createFieldSpace<int>();
//
//	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
//	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();
//
//	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
//	auto p = Pimpact::createScalarField<double,int>(fS);
//
//
//	auto Xp = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,1);
//	auto Bp = Xp->clone();
//	auto X = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,1);
//	auto B = X->clone();
//
//	X->init(0.);
//	B->random();
//	B->scale(0.1);
//
//	Xp->init(0.);
//	Xp->random();
//	Bp->random();
//	Bp->scale(0.1);
//
//	// init operators
//	auto lap = Pimpact::createHelmholtzdep<double,int>( 0.,1.);
//
//	// Set some Belos parameters.
//  auto solverName = "CG";
//  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-3 );
//  solverParams->get()->set ("Verbosity", int( Belos::Errors) );
//	solverParams->get()->set ("Maximum Iterations", 40);
//
//// Create the Pimpact::LinearSolver solver.
//	auto H_prob = Pimpact::createLinearProblem<Scalar,MV,OP>( lap, X, B, solverParams->get(), solverName );
//
//	auto schur = Pimpact::createDivHinvGrad<double,int>( X->clone(),
////			div, grad,
//			H_prob );
//
//  solverName = "GMRES";
//  solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
//  solverParams->get()->set( "Maximum Iterations", 10 );
//
//
//	auto schur_prob = Pimpact::createLinearProblem<Scalar,MVp,OPp>( schur, Xp,Bp, solverParams->get(), solverName);
//
//	Belos::ReturnType result = schur_prob->solve(Xp,Bp);
//
//	TEST_EQUALITY( result,Belos::Converged);
//
//	schur_prob->apply( Xp, Bp );
//
//	Xp->write(200);
//
//}
//
//
//TEUCHOS_UNIT_TEST( BelosSolver, Div_DtLinv_Grad ) {
//	using Teuchos::ParameterList;
//	using Teuchos::parameterList;
//	using Teuchos::RCP;
//	using Teuchos::rcp; // Save some typing
//
//	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MVF;
//	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > > MSF;
//	typedef Pimpact::OperatorMV< Pimpact::DtL<double,int> >  OP;
//	typedef Pimpact::OperatorMV< Pimpact::Div_DtLinv_Grad<double,int> >  OP2;
//
//	auto temp = Pimpact::createMultiModeVectorField<double,int>();
//
//	auto X = Pimpact::createMultiModeScalarField<double,int>();
//	auto B = Pimpact::createMultiModeScalarField<double,int>();
//
//	X->init(0.);
//	B->init(0.);
//	B->random(0.);
//	B->scale(0.1);
//
//	auto A = Pimpact::createDtL<double,int>(1.,0.,1.);
//
//	// Make an empty new parameter list.
//	auto solverName = "GMRES";
//	auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
//  solverParams->get()->set ("Verbosity",  Belos::Errors );
//
//// Create the Pimpact::LinearSolver solver.
//	auto H_prob = Pimpact::createLinearProblem<double,MVF,OP>( A, temp, temp, solverParams->get(), solverName );
//
//	auto schur = Pimpact::createDivDtLinvGrad<double,int>( temp, H_prob );
//
//	solverParams->get()->set ("Verbosity",  Belos::Errors + Belos::Warnings +
//	Belos::TimingDetails + Belos::StatusTestDetails);
//	solverParams->get()->set ("Maximum Iterations", 10);
//
//	auto schur_prob = Pimpact::createLinearProblem<double,MSF,OP2>( schur, X, B, solverParams->get(), solverName );
//
//	schur_prob->solve( X, B );
//
//	schur_prob->apply( X, B );
//
//}
//
//
//TEUCHOS_UNIT_TEST( BelosSolver, Dt1L1 ) {
//  typedef int O;
//  typedef double S;
//
//  typedef Pimpact::VectorField<S,O> VF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField<MVF> BVF;
//
//  typedef Pimpact::DtL<S,O> Op;
//  typedef Pimpact::DtModeOp<S,O> Op2;
////  typedef Pimpact::OperatorMV<Op> OpMV;
////  typedef Pimpact::OperatorMV<Op2> OpMV2;
//  typedef Pimpact::OperatorBase<BVF> OpBase;
////  typedef Pimpact::OperatorPimpl<BVF,Op> OpPimpl;
////  typedef Pimpact::OperatorPimpl<BVF,Op> prec;
//
//  auto fS = Pimpact::createFieldSpace<O>();
//
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();
//
//  auto fieldc = Pimpact::createVectorField<S,O>(fS,iIS,fIS);
//  auto fields = Pimpact::createVectorField<S,O>(fS,iIS,fIS);
//
//  auto field = Pimpact::createModeField( fieldc, fields );
//
//  auto x = Pimpact::createInitMVF( Pimpact::Streaming2DFlow, fS, iIS, fIS,1.,1.,1.);
//  auto b = Pimpact::createInitMVF( Pimpact::Zero2DFlow, fS, iIS, fIS,1.,1.,1.);
//
//  auto op = Pimpact::createOperatorBaseMV<BVF,Op>( Pimpact::createDtL<S,O>( 100., 0., 1. ) );
////  auto prec = Pimpact::createOperatorBase<BVF,Op2>( Pimpact::createDtModeOp<S,O>( -1. ) );
//
//
//  // Make an empty new parameter list.
//  auto para = Pimpact::createLinSolverParameter("GMRES");
//
//  // Create the Pimpact::LinearSolver solver.
//  auto pimpact_problem = Pimpact::createLinearProblem<S,BVF,OpBase>( op, x, b, para->get(),"GMRES" );
////  pimpact_problem->setLeftPrec(prec);
//
//  pimpact_problem->solve(x,b);
//  TEST_EQUALITY( pimpact_problem->getSolver()->achievedTol()<1.e-8, true );
//
//
//  pimpact_problem->apply( x, b );
//
//
//  x->write(200);
//}


} // end of namespace

