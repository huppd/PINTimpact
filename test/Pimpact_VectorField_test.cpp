#include <iostream>
#include <vector>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosTypes.hpp"

//#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"



namespace {

bool testMpi = true;
double eps = 1e-6;
int domain = 1;
int dim = 2;

auto pl = Teuchos::parameterList();

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
  clp.setOption(
      "dim", &dim,
      "dim" );

  pl->set( "domain", domain);
  pl->set( "dim", dim);
  pl->set( "lx", 2. );
  pl->set( "ly", 2. );
  pl->set( "lz", 1. );


  pl->set("nx", (48*3)+1 );
  pl->set("ny", 49 );
  pl->set("nz", 17 );

  //  // processor grid size
  //  pl->set("npx", 2 );
  //  pl->set("npy", 2 );
  //  pl->set("npz", 1 );

}


TEUCHOS_UNIT_TEST( VectorField, create_init_print ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  vel->init( space->rankST() );
  vel->print();

}



TEUCHOS_UNIT_TEST( VectorField, OneNorm_and_init ) {

  auto space = Pimpact::createSpace( pl, false );

  auto x = Pimpact::create<Pimpact::VectorField>(space);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    x->init(i/2.);
    TEST_EQUALITY( (i/2.)*x->getLength(), x->norm(Belos::OneNorm) );
  }
}




TEUCHOS_UNIT_TEST( VectorField, InfNorm_and_init ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

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
  MPI_Comm_rank( space->comm(), &rank );
  for( double i = 0.; i<200.1; ++i) {
    init = 3*i-1.;
    init = (init<0)?-init:init;
    vel->init(rank*i-1.);
    norm = vel->norm(Belos::InfNorm);
    TEST_EQUALITY( init, norm );
  }
}


TEUCHOS_UNIT_TEST( VectorField, InfNorm_and_initvec2d ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  double norm;


  auto alpha0 = Teuchos::tuple( 1., 0., 0. );
  vel->init(alpha0);
  norm = vel->norm(Belos::InfNorm);
  TEST_FLOATING_EQUALITY( 1., norm, eps );

  auto alpha1 = Teuchos::tuple( 0., 1., 0. );
  vel->init(alpha1);
  norm = vel->norm(Belos::InfNorm);
  TEST_FLOATING_EQUALITY( 1., norm, eps );

  auto alpha2 = Teuchos::tuple( 0., 0., 1. );
  vel->init(alpha2);
  norm = vel->norm(Belos::InfNorm);
  TEST_FLOATING_EQUALITY( (2==dim)?0.:1., norm, eps );
}


TEUCHOS_UNIT_TEST( VectorField, TwoNorm_and_init ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  double norm;
  int N = vel->getLength();

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    vel->init(i/2.);
    norm = vel->norm(Belos::TwoNorm);
    TEST_EQUALITY( std::sqrt( std::pow(i/2.,2)*N ), norm );
  }
}


TEUCHOS_UNIT_TEST( VectorField, dot ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel1 = Pimpact::create<Pimpact::VectorField>(space);
  auto vel2 = Pimpact::create<Pimpact::VectorField>(space);

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



TEUCHOS_UNIT_TEST( VectorField, scale ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  double N = vel->getLength();
  double norm;

  vel->init(1.);
  vel->scale(2.);
  norm = vel->norm(Belos::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(4*N), norm, eps )

}



TEUCHOS_UNIT_TEST( VectorField, random ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  double N = vel->getLength();
  double norm;

  vel->init(1.);
  vel->random();
  norm = vel->norm(Belos::TwoNorm);
  TEST_INEQUALITY( N, norm )

}



TEUCHOS_UNIT_TEST( VectorField, add ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel1 = Pimpact::create<Pimpact::VectorField>(space);
  auto vel2 = Pimpact::create<Pimpact::VectorField>(space);
  auto vel3 = Pimpact::create<Pimpact::VectorField>(space);

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



TEUCHOS_UNIT_TEST( VectorField, write ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel = Pimpact::create<Pimpact::VectorField>(space);

  vel->init( 1. );
  vel->write();

  vel->random();
  vel->write( 1 );

  TEST_EQUALITY( 0, 0 )
}

TEUCHOS_UNIT_TEST( VectorField, write_restart ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel = Pimpact::create<Pimpact::VectorField>(space);

  vel->init( 1. );
  vel->write();

  vel->random();
  vel->write( 99, true );

}



TEUCHOS_UNIT_TEST( VectorField, initField ) {

  auto space = Pimpact::createSpace( pl, false );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  for( int i=0; i<=19; ++i ) {
    vel->initField( Pimpact::EFlowField(i) );
    vel->write( i );
  }

}


} // end of namespace
