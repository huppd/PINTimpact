// Pimpact_SalarVector_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "BelosTypes.hpp"

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_FieldFactory.hpp"
#include "Pimpact_OperatorMV.hpp"
#include "Pimpact_Operator.hpp"
#include "BelosPimpactAdapter.hpp"

#include <iostream>
#include <cmath>


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

	TEUCHOS_UNIT_TEST( Operator, Div ) {
		// init impact
    init_impact(0,0);

		auto fS = Pimpact::createFieldSpace<int>();
		auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
		auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

		auto p = Pimpact::createScalarField<double,int>(fS);
		auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

//		vel->init(rank);
//		vel->print();
		vel->init(1.);
		p->init(2.);

		Pimpact::Div<double,int> div;

		div.apply(*vel,*p);

		p->write(1);

//		TEST_EQUALITY( 0., p->norm() );
	}

	TEUCHOS_UNIT_TEST( Operator, Div2 ) {

		auto fS = Pimpact::createFieldSpace<int>();
		auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
		auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

		auto p = Pimpact::createScalarField<double,int>(fS);
		auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

		vel->init_field();
		p->init(2.);

		Pimpact::Div<double,int> div;

		div.apply(*vel,*p);

		p->write(2);

		TEST_INEQUALITY( 0., p->norm() );
	}


TEUCHOS_UNIT_TEST( Operator, Grad ) {

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto p = Pimpact::createScalarField<double,int>(fS);
	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	Pimpact::Grad<double,int> grad;

	grad.apply(*p,*vel);


	TEST_EQUALITY( 0, 0 );
}


	TEUCHOS_UNIT_TEST( Operator, Helmholtz ) {

			auto fS = Pimpact::createFieldSpace<int>();
			auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
			auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

//			auto p = Pimpact::createScalarField<double,int>(fS);
			auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
			auto res = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

			Pimpact::Helmholtz<double,int> helmholtz;

			helmholtz.apply(*vel,*res);


			TEST_EQUALITY( 0, 0 );
		}


TEUCHOS_UNIT_TEST( Operator, HelmholtzMV ) {

		auto fS = Pimpact::createFieldSpace<int>();
		auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
		auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

//			auto p = Pimpact::createScalarField<double,int>(fS);
		auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
		auto res = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

		auto mv = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,10);
//			auto mv2 = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,10);
		auto mv2 = mv->Clone(11);

		Pimpact::Helmholtz<double,int> helmholtz;
		auto helm = Pimpact::createOperatorMV<Pimpact::Helmholtz<double,int> >();

		helmholtz.apply(*vel,*res);
		helm->apply(*mv,*mv2);
		Belos::OperatorTraits<double,Pimpact::MultiField<Pimpact::VectorField<double,int> >, Pimpact::OperatorMV<Pimpact::Helmholtz<double,int> > >::Apply(*helm,*mv,*mv2);


		TEST_EQUALITY( 0, 0 );
	}


TEUCHOS_UNIT_TEST( Operator, HelmholtzMVMode ) {

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> >,double,int>(*vel,10);

	auto mv2 = mv->Clone(11);

	Pimpact::Helmholtz<double,int> helmholtz;
	auto helm = Pimpact::createOperatorMV<Pimpact::Helmholtz<double,int> >();

	helm->apply(*mv,*mv2);
	Belos::OperatorTraits<double,Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >, Pimpact::OperatorMV<Pimpact::Helmholtz<double,int> > >::Apply(*helm,*mv,*mv2);

	TEST_EQUALITY( 0, 0 );
}


TEUCHOS_UNIT_TEST( Operator, Dt ) {

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> >,double,int>(*vel,1);

	mv->GetVec(0).getFieldC()->init(1.);
	mv->GetVec(0).getFieldS()->init(0.);

	auto mv2 = mv->Clone(1);
	mv2->Init(0.);

	auto A = Pimpact::createOperatorMV<Pimpact::Dt<double,int> >();

	A->apply(*mv,*mv2);

	TEST_EQUALITY( mv->GetConstVec(0).getConstFieldC()->norm(), mv2->GetConstVec(0).getConstFieldS()->norm() );
	TEST_EQUALITY( mv->GetConstVec(0).getConstFieldS()->norm(), mv2->GetConstVec(0).getConstFieldC()->norm() );
	Belos::OperatorTraits<double,Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >, Pimpact::OperatorMV<Pimpact::Dt<double,int> > >::Apply(*A,*mv,*mv2);
}


TEUCHOS_UNIT_TEST( Operator, DtL ) {

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> >,double,int>(*vel,1);

	mv->Random();

	auto mv2 = mv->Clone(1);
	auto mv3 = mv->Clone(1);
	mv2->Init(0.);
	mv3->Init(0.);

	auto A = Pimpact::createDtL<double,int>( 1., 0., 0. );

	A->apply(*mv,*mv2);

	TEST_EQUALITY( mv->GetConstVec(0).getConstFieldC()->norm(), mv2->GetConstVec(0).getConstFieldS()->norm() );
	TEST_EQUALITY( mv->GetConstVec(0).getConstFieldS()->norm(), mv2->GetConstVec(0).getConstFieldC()->norm() );

	Belos::OperatorTraits<double,Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >, Pimpact::OperatorMV<Pimpact::DtL<double,int> > >::Apply(*A,*mv,*mv2);

	mv2->Init(0.);

	auto A2 = Pimpact::createDtL<double,int>( 0., 0., 1. );
	auto A3 = Pimpact::createHelmholtz<double,int>( 0., 1. );

	A2->apply(*mv,*mv2);
	A3->apply(*mv,*mv3);

	std::vector<double> norm2( 1, 0. );
	std::vector<double> norm3( 1, 0. );

	mv2->Norm(norm2);
	mv3->Norm(norm3);

	TEST_EQUALITY( norm2[0] , norm3[0] );

	mv2->Init(0.);

	A2 = Pimpact::createDtL<double,int>( 0., 1., 0. );
	A3 = Pimpact::createHelmholtz<double,int>( 1., 0. );

	A2->apply(*mv,*mv2);
	A3->apply(*mv,*mv3);

	mv2->Norm(norm2);
	mv3->Norm(norm3);

	TEST_EQUALITY( norm2[0] , norm3[0] );

}

TEUCHOS_UNIT_TEST( Operator, Div_DtLinv_Grad ) {
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing

	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MVF;
	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > > MSF;
	typedef Pimpact::OperatorMV< Pimpact::DtL<double,int> >  OP;

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
	solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
	Belos::TimingDetails + Belos::StatusTestDetails);

// Create the Pimpact::LinearSolver solver.
	auto H_prob = Pimpact::createLinearProblem<double,MVF,OP>( A, temp, temp, solverParams,"GMRES" );

	auto schur = Pimpact::createDivDtLinvGrad<double,int>( temp, H_prob );

	schur->apply( *B, *X );

}

} // namespace

