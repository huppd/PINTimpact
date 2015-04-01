#include <iostream>
//#include <vector>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"
//#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"



namespace {

bool testMpi = true;
double eps = 1e-6;
int domain = 1;
int dim = 3;

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

  pl->set( "domain", domain );
  pl->set( "dim", dim );
  pl->set( "lx", 2. );
  pl->set( "ly", 2. );
  pl->set( "lz", 1. );


  pl->set("nx", (48*3)+1 );
  pl->set("ny", 49 );
  pl->set("nz", 17 );

  // processor grid size
  pl->set("npx", 2 );
  pl->set("npy", 2 );
  pl->set("npz", 2 );

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
  TEST_FLOATING_EQUALITY( ( 2==space->dim() )?0.:1., norm, eps );
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
