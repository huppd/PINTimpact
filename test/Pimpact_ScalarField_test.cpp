// Pimpact_SalarVector_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "BelosTypes.hpp"

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_ScalarField.hpp"

#include <iostream>
#include <cmath>

//#include "hdf5.h"


namespace {

bool testMpi = true;
double errorTolSlack = 1e-6;
int domain = 1;

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
  clp.setOption(
      "domain", &domain,
      "domain" );
}



TEUCHOS_UNIT_TEST( ScalarField, create_init_print ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

  space->print();

  auto p = Pimpact::createScalarField(space);

  p->init( space->rankST() );

}



TEUCHOS_UNIT_TEST( ScalarField, InfNorm_and_init ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

  auto p = Pimpact::createScalarField(space);

  double norm;

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    p->init(i/2.);
    norm = p->norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( i/2., norm, errorTolSlack );

  }

  // one test with infty-norm
  int rank;
  double init;
  MPI_Comm_rank(space->comm(),&rank);
  for( double i = 0.; i<200.1; ++i) {
    init = 3*i-1.;
    init = (init<0)?-init:init;
    p->init(rank*i-1.);
    norm = p->norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( init, norm, errorTolSlack );

  }
}



TEUCHOS_UNIT_TEST( ScalarField, OneNorm_and_init ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

  auto p = Pimpact::createScalarField(space);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    p->init(i/2.);
//    TEST_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm) );
    TEST_FLOATING_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm), errorTolSlack );

  }
}



TEUCHOS_UNIT_TEST( ScalarField, TwoNorm_and_init ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

  auto p = Pimpact::createScalarField(space);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    p->init(i/2.);
    TEST_EQUALITY( std::sqrt( std::pow(i/2.,2)*p->getLength() ), p->norm(Belos::TwoNorm) );
  }
}


TEUCHOS_UNIT_TEST( ScalarField, dot ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, true );

  auto p = Pimpact::createScalarField( space );
  auto q = Pimpact::createScalarField( space );

  int Np = p->getLength();
  int Nq = q->getLength();
  double dot;

  TEST_EQUALITY( Np , Nq );

  p->init(1.);
  q->init(1.);
  dot = p->dot(*q);

  TEST_FLOATING_EQUALITY( 1.*Np, dot, errorTolSlack );

}



TEUCHOS_UNIT_TEST( ScalarField, scale ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

  auto p = Pimpact::createScalarField(space);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->scale(2.);
  norm = p->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(4*N), norm)
}


TEUCHOS_UNIT_TEST( ScalarField, random ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

  auto p = Pimpact::createScalarField(space);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->random();
  norm = p->norm(Belos::TwoNorm);
  TEST_INEQUALITY( N, norm)
}


TEUCHOS_UNIT_TEST( ScalarField, add ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

  auto q = Pimpact::createScalarField(space);

  auto r(q);
  auto p(q);

  double norm;
  int N = p->getLength();

  q->init(1.);
  r->init(1./3.);

  p->add( 0., *q, 3., *r);
  norm = p->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm)
}

TEUCHOS_UNIT_TEST( ScalarField, write ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, false );

  auto p = Pimpact::createScalarField(space);

  p->init(1.);
  p->write();

  p->random();
  p->write(1);

  TEST_EQUALITY( 0, 0)
}


TEUCHOS_UNIT_TEST( ScalarField, initField ) {
//  H5open();

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain );

  pl->set( "lx", 2. );
  pl->set( "ly", 2. );
  pl->set( "lz", 1. );


  pl->set("nx", (48*3)+1 );
  pl->set("ny", 49 );
  pl->set("nz", 2 );

//  // processor grid size
//  pl->set("npx", 2 );
//  pl->set("npy", 2 );
//  pl->set("npz", 1 );

  auto space = Pimpact::createSpace( pl, false );

  auto x = Pimpact::createScalarField( space );

  for( int i=0; i<=4; ++i ) {
    x->initField( Pimpact::EScalarField(i) );
    x->write( i );
  }

}


} // namespace

