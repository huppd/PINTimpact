#include <iostream>
#include <vector>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "BelosTypes.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

#include "Pimpact_CompoundField.hpp"



namespace {

bool testMpi = true;
double eps = 1e-6;
int domain = 1;

TEUCHOS_STATIC_SETUP() {
	Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
	clp.addOutputSetupOptions(true);
	clp.setOption(
			"test-mpi", "test-serial", &testMpi,
			"Test MPI (if available) or force test of serial.  In a serial build,"
			" this option is ignored and a serial comm is always used." );
	clp.setOption(
			"error-tol-slack", &eps,
			"Slack off of machine epsilon used to check test results" );
  clp.setOption(
      "domain", &domain,
      "domain" );
}



TEUCHOS_UNIT_TEST( VectorModeField, create_init_print ) {
  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

	auto vel = Pimpact::createVectorField( space );
	auto p = Pimpact::createScalarField( space );

	auto um = Pimpact::createCompoundField( vel, p );

	um->init( space->rankST() );
}



TEUCHOS_UNIT_TEST( VectorModeField, InfNorm_and_init ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

	auto vel = Pimpact::createVectorField( space );
	auto p = Pimpact::createScalarField( space );

	auto um = Pimpact::createCompoundField( vel, p );

	double norm;
	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		um->init(i/2.);
		norm = um->norm(Belos::InfNorm);
		TEST_EQUALITY( i/2., norm );
	}

	// one test with infty-norm
	int rank;
	double init;
	MPI_Comm_rank(space->comm(),&rank);
	for( double i = 0.; i<200.1; ++i) {
		init = 3*i-1.;
		init = (init<0)?-init:init;
		um->init(rank*i-1.);
		norm = um->norm(Belos::InfNorm);
		TEST_EQUALITY( init, norm );
	}
}



TEUCHOS_UNIT_TEST( VectorModeField, TwoNorm_and_init ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

	auto vel = Pimpact::createVectorField( space );
	auto p = Pimpact::createScalarField( space );

	auto q = Pimpact::createCompoundField( vel, p );

	double norm;
	int N = q->getLength();

	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		q->init(i/2.);
		norm = q->norm(Belos::TwoNorm);
    TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*N), norm, eps );
	}

}



TEUCHOS_UNIT_TEST( VectorModeField, dot ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

	auto vel = Pimpact::createVectorField( space );
	auto p = Pimpact::createScalarField( space );
	auto vel3 = Pimpact::createVectorField( space );
	auto p3 = Pimpact::createScalarField( space );

	auto vel1 = Pimpact::createCompoundField( vel, p );
	auto vel2 = Pimpact::createCompoundField( vel3, p3 );

	double dot;

	TEST_EQUALITY( vel1->getLength(), vel2->getLength() )

	int N = vel1->getLength();

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

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

	auto vel = Pimpact::createVectorField( space );
	auto p = Pimpact::createScalarField( space );

	auto q = Pimpact::createCompoundField( vel, p );

	int N = q->getLength();
	double norm;

	q->init(1.);
	q->scale(2.);
	norm = q->norm(Belos::TwoNorm);
	TEST_EQUALITY( std::sqrt(4*N), norm)

}



TEUCHOS_UNIT_TEST( VectorModeField, random ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

	auto vel = Pimpact::createVectorField( space );
	auto p = Pimpact::createScalarField( space );

	auto q = Pimpact::createCompoundField( vel, p );

	int N = q->getLength();
	double norm;

	q->init(1.);
	q->random();
	norm = q->norm(Belos::TwoNorm);
	TEST_INEQUALITY( N, norm)

}



TEUCHOS_UNIT_TEST( VectorModeField, add ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

	auto vel1c = Pimpact::createVectorField(space);
	auto vel2c = Pimpact::createVectorField(space);
	auto vel3c = Pimpact::createVectorField(space);

	auto vel1s = Pimpact::createScalarField( space );
	auto vel2s = Pimpact::createScalarField( space );
	auto vel3s = Pimpact::createScalarField( space );

	auto vel1 = Pimpact::createCompoundField( vel1c, vel1s );
	auto vel2 = Pimpact::createCompoundField( vel2c, vel2s );
	auto vel3 = Pimpact::createCompoundField( vel3c, vel3s );


	TEST_EQUALITY( vel1->getLength(), vel2->getLength() )
	TEST_EQUALITY( vel2->getLength(), vel3->getLength() )
	TEST_EQUALITY( vel1->getLength(), vel3->getLength() )

	double norm;
	int N = vel1->getLength();

	vel1->init(0.);
	vel2->init(1./2.);
	vel3->init(1./3.);

	vel1->add( 2., *vel2, 0., *vel3);
	norm = vel1->norm(Belos::TwoNorm);
	TEST_EQUALITY( std::sqrt(N), norm )

	vel1->init(0.);
	vel2->init(1./2.);
	vel3->init(1./3.);

	vel1->add( 0., *vel2, 3., *vel3);
	norm = vel1->norm(Belos::TwoNorm);
	TEST_EQUALITY( std::sqrt(N), norm )

	vel1->init(0.);
	vel2->init(1.);
	vel3->init(1.);

	vel1->add( 0.5, *vel2, 0.5, *vel3);
	norm = vel1->norm(Belos::TwoNorm);
	TEST_EQUALITY( std::sqrt(N), norm )

}



TEUCHOS_UNIT_TEST( VectorModeField, write ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

	auto vel = Pimpact::createVectorField( space );
	auto p = Pimpact::createScalarField( space );

	auto q = Pimpact::createCompoundField( vel, p );

	q->init( 1. );
	q->write();

	q->random();
	q->write( 2 );

}


} // end of namespace
