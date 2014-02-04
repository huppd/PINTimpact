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
#include "Pimpact_ModeField.hpp"

#include <iostream>
#include <vector>
#include <cmath>


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


TEUCHOS_UNIT_TEST( VectorModeField, create_init_print ) {
	// init impact
	int rank = init_impact(0,0);

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	vel->init(rank);
}


TEUCHOS_UNIT_TEST( VectorModeField, InfNorm_and_init ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	double norm;
	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		vel->init(i/2.);
		norm = vel->norm(Belos::InfNorm);
		TEST_EQUALITY( i/2., norm );
	}

	// one test with infty-norm
	int rank;
	double init;
	MPI_Comm_rank(sVS->comm_,&rank);
	for( double i = 0.; i<200.1; ++i) {
		init = 3*i-1.;
		init = (init<0)?-init:init;
		vel->init(rank*i-1.);
		norm = vel->norm(Belos::InfNorm);
		TEST_EQUALITY( init, norm );
	}
}



TEUCHOS_UNIT_TEST( VectorModeField, TwoNorm_and_init ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	double norm;
	int N = vel->getVecLength();

	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		vel->init(i/2.);
		norm = vel->norm(Belos::TwoNorm);
		TEST_EQUALITY( std::pow(i/2.,2)*N, norm );
	}
}


TEUCHOS_UNIT_TEST( VectorModeField, dot ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto p1c = Pimpact::createScalarField<double,int>(sVS);
	auto p1s = Pimpact::createScalarField<double,int>(sVS);

	auto p2c = Pimpact::createScalarField<double,int>(sVS);
	auto p2s = Pimpact::createScalarField<double,int>(sVS);

	auto vel1 = Pimpact::createModeField( p1c, p1s );
	auto vel2 = Pimpact::createModeField( p2c, p2s );

	double dot;

	TEST_EQUALITY( vel1->getVecLength(), vel2->getVecLength() )

	int N = vel1->getVecLength();

	vel1->init(0.);
	vel2->init(1.);
	dot = vel1->dot(*vel2);
	TEST_EQUALITY( 0, dot );

	vel1->init(1.);
	vel2->init(1.);
	dot = vel2->dot(*vel1);
	TEST_EQUALITY( N, dot );

	vel1->init(2.);
	vel2->init(1.);
	dot = vel1->dot(*vel2);
	TEST_EQUALITY( 2*N, dot );

	vel1->init(1.);
	vel2->init(2.);
	dot = vel1->dot(*vel2);
	TEST_EQUALITY( 2*N, dot );

}


TEUCHOS_UNIT_TEST( VectorModeField, scale ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	int N = vel->getVecLength();
	double norm;

	vel->init(1.);
	vel->scale(2.);
	norm = vel->norm(Belos::TwoNorm);
	TEST_EQUALITY( 4*N, norm)
}


TEUCHOS_UNIT_TEST( VectorModeField, random ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	int N = vel->getVecLength();
	double norm;

	vel->init(1.);
	vel->random();
	norm = vel->norm(Belos::TwoNorm);
	TEST_INEQUALITY( N, norm)
}


TEUCHOS_UNIT_TEST( VectorModeField, add ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto p1c = Pimpact::createScalarField<double,int>(sVS);
	auto p1s = Pimpact::createScalarField<double,int>(sVS);

	auto p2c = Pimpact::createScalarField<double,int>(sVS);
	auto p2s = Pimpact::createScalarField<double,int>(sVS);

	auto p3c = Pimpact::createScalarField<double,int>(sVS);
	auto p3s = Pimpact::createScalarField<double,int>(sVS);

	auto vel1 = Pimpact::createModeField( p1c, p1s );
	auto vel2 = Pimpact::createModeField( p2c, p2s );
	auto vel3 = Pimpact::createModeField( p3c, p3s );



	TEST_EQUALITY( vel1->getVecLength(), vel2->getVecLength() )
	TEST_EQUALITY( vel2->getVecLength(), vel3->getVecLength() )
	TEST_EQUALITY( vel1->getVecLength(), vel3->getVecLength() )

	double norm;
	int N = vel1->getVecLength();

	vel1->init(0.);
	vel2->init(1./2.);
	vel3->init(1./3.);

	vel1->add( 2., *vel2, 0., *vel3);
	norm = vel1->norm(Belos::TwoNorm);
	TEST_EQUALITY( N, norm)

	vel1->init(0.);
	vel2->init(1./2.);
	vel3->init(1./3.);

	vel1->add( 0., *vel2, 3., *vel3);
	norm = vel1->norm(Belos::TwoNorm);
	TEST_EQUALITY( N, norm)

	vel1->init(0.);
	vel2->init(1.);
	vel3->init(1.);

	vel1->add( 0.5, *vel2, 0.5, *vel3);
	norm = vel1->norm(Belos::TwoNorm);
	TEST_EQUALITY( N, norm )
}


TEUCHOS_UNIT_TEST( VectorModeField, write ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	vel->init( 1. );
	vel->write();

	vel->random();
	vel->write( 2 );

	TEST_EQUALITY( 0, 0 )
}

} // namespace

