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


namespace {

bool testMpi = true;
double errorTolSlack = 1e+1;
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
  int rank = Pimpact::init_impact_pre();
  auto bc = Pimpact::createBoudaryConditionsGlobal( Pimpact::EDomainType( domain ) );
  bc->print();
  Pimpact::init_impact_post();  // init impact

  auto space = Pimpact::createSpace();
  space->print();

  auto p = Pimpact::createScalarField(space);

  p->init(rank);

  TEST_EQUALITY( 0, 0 );
}



TEUCHOS_UNIT_TEST( ScalarField, InfNorm_and_init ) {

  auto sVS = Pimpact::createSpace();

  auto p = Pimpact::createScalarField(sVS);

  double norm;

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    p->init(i/2.);
    norm = p->norm(Belos::InfNorm);
    TEST_EQUALITY( i/2., norm );
  }

  // one test with infty-norm
  int rank;
  double init;
  MPI_Comm_rank(sVS->comm(),&rank);
  for( double i = 0.; i<200.1; ++i) {
    init = 3*i-1.;
    init = (init<0)?-init:init;
    p->init(rank*i-1.);
    norm = p->norm(Belos::InfNorm);
    TEST_EQUALITY( init, norm );
  }
}



TEUCHOS_UNIT_TEST( ScalarField, OneNorm_and_init ) {

  auto sVS = Pimpact::createSpace();
  auto p = Pimpact::createScalarField(sVS);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    p->init(i/2.);
    TEST_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm) );
  }
}



TEUCHOS_UNIT_TEST( ScalarField, TwoNorm_and_init ) {

  auto sVS = Pimpact::createSpace();
  auto p = Pimpact::createScalarField(sVS);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    p->init(i/2.);
    TEST_EQUALITY( std::sqrt( std::pow(i/2.,2)*p->getLength() ), p->norm(Belos::TwoNorm) );
  }
}


TEUCHOS_UNIT_TEST( ScalarField, dot ) {

  auto fS = Pimpact::createSpace();

  auto p = Pimpact::createScalarField(fS);
  auto q = Pimpact::createScalarField(fS);

  int Np = p->getLength();
  int Nq = q->getLength();
  double dot;

  TEST_EQUALITY( Np , Nq );

  p->init(1.);
  q->init(2.);
  dot = p->dot(*q);
  TEST_EQUALITY( 2*Np, dot );
}



TEUCHOS_UNIT_TEST( ScalarField, scale ) {

  auto sVS = Pimpact::createSpace();
  auto p = Pimpact::createScalarField(sVS);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->scale(2.);
  norm = p->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(4*N), norm)
}


TEUCHOS_UNIT_TEST( ScalarField, random ) {

  auto sVS = Pimpact::createSpace();
  auto p = Pimpact::createScalarField(sVS);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->random();
  norm = p->norm(Belos::TwoNorm);
  TEST_INEQUALITY( N, norm)
}


TEUCHOS_UNIT_TEST( ScalarField, add ) {

  auto sVS = Pimpact::createSpace();

  auto q = Pimpact::createScalarField(sVS);

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

  auto sVS = Pimpact::createSpace();
  auto p = Pimpact::createScalarField(sVS);

  p->init(1.);
  p->write();

  p->random();
  p->write(1);

  TEST_EQUALITY( 0, 0)
}
} // namespace

