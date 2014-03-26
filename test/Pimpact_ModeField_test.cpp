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
#include "Pimpact_ModeField.hpp"

#include <iostream>
#include <vector>
#include <cmath>




namespace {


bool testMpi = true;
double errorTolSlack = 1e-6;



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



TEUCHOS_UNIT_TEST( ScalarModeField, create_init_print ) {
	// init impact
	int rank = init_impact(0,0);

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	vel->init(rank);
}



TEUCHOS_UNIT_TEST( ScalarModeField, InfNorm_and_init ) {

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



TEUCHOS_UNIT_TEST( ScalarModeField, TwoNorm_and_init ) {
	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	double norm;
	int N = vel->getLength();

	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		vel->init(i/2.);
		norm = vel->norm(Belos::TwoNorm);
		TEST_EQUALITY( std::sqrt(std::pow(i/2.,2)*N), norm );
	}
}



TEUCHOS_UNIT_TEST( ScalarModeField, dot ) {
	auto sVS = Pimpact::createFieldSpace<int>();

	auto p1c = Pimpact::createScalarField<double,int>(sVS);
	auto p1s = Pimpact::createScalarField<double,int>(sVS);

	auto p2c = Pimpact::createScalarField<double,int>(sVS);
	auto p2s = Pimpact::createScalarField<double,int>(sVS);

	auto vel1 = Pimpact::createModeField( p1c, p1s );
	auto vel2 = Pimpact::createModeField( p2c, p2s );

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



TEUCHOS_UNIT_TEST( ScalarModeField, scale ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	int N = vel->getLength();
	double norm;

	vel->init(1.);
	vel->scale(2.);
	norm = vel->norm(Belos::TwoNorm);
	TEST_EQUALITY( std::sqrt(4*N), norm)

}



TEUCHOS_UNIT_TEST( ScalarModeField, random ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	int N = vel->getLength();
	double norm;

	vel->init(1.);
	vel->random();
	norm = vel->norm(Belos::TwoNorm);
	TEST_INEQUALITY( std::sqrt(N), norm)

}



TEUCHOS_UNIT_TEST( ScalarModeField, add ) {

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



TEUCHOS_UNIT_TEST( ScalarModeField, write ) {

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



TEUCHOS_UNIT_TEST( VectorModeField, create_init_print ) {

  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto um = Pimpact::createModeField( velc, vels );

  um->init();
}



TEUCHOS_UNIT_TEST( VectorModeField, InfNorm_and_init ) {

  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto um = Pimpact::createModeField( velc, vels );

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
  MPI_Comm_rank(fS->comm_,&rank);
  for( double i = 0.; i<200.1; ++i) {
    init = 3*i-1.;
    init = (init<0)?-init:init;
    um->init(rank*i-1.);
    norm = um->norm(Belos::InfNorm);
    TEST_EQUALITY( init, norm );
  }
}



TEUCHOS_UNIT_TEST( VectorModeField, TwoNorm_and_init ) {

  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto vel = Pimpact::createModeField( velc, vels );

  double norm;
  int N = vel->getLength();

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    vel->init(i/2.);
    norm = vel->norm(Belos::TwoNorm);
    TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*N), norm, errorTolSlack );
  }
}



TEUCHOS_UNIT_TEST( VectorModeField, dot ) {

  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto vel1c = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vel2c = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto vel1s = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vel2s = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto vel1 = Pimpact::createModeField( vel1c, vel1s );
  auto vel2 = Pimpact::createModeField( vel2c, vel2s );

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

  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto vel = Pimpact::createModeField( velc, vels );

  int N = vel->getLength();
  double norm;

  vel->init(1.);
  vel->scale(2.);
  norm = vel->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(4*N), norm)

}


TEUCHOS_UNIT_TEST( VectorModeField, random ) {

  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto vel = Pimpact::createModeField( velc, vels );

  int N = vel->getLength();
  double norm;

  vel->init(1.);
  vel->random();
  norm = vel->norm(Belos::TwoNorm);
  TEST_INEQUALITY( std::sqrt(N), norm)

}


TEUCHOS_UNIT_TEST( VectorModeField, add ) {

  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto vel1c = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vel2c = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vel3c = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto vel1s = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vel2s = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vel3s = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto vel1 = Pimpact::createModeField( vel1c, vel1s );
  auto vel2 = Pimpact::createModeField( vel2c, vel2s );
  auto vel3 = Pimpact::createModeField( vel3c, vel3s );


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

  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  auto vel = Pimpact::createModeField( velc, vels );

  vel->init( 1. );
  vel->write();

  vel->random();
  vel->write( 2 );

  TEST_EQUALITY( 0, 0 )

}


} // end of namespace
