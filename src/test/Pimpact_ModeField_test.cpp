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
#include "Pimpact_ModeField.hpp"



namespace {

typedef double S;
typedef int O;

typedef typename Pimpact::Space<S,O,3,4> SpaceT;
typedef typename Pimpact::ScalarField<SpaceT> SF;
typedef typename Pimpact::VectorField<SpaceT> VF;
typedef typename Pimpact::ModeField<SF>       MSF;
typedef typename Pimpact::ModeField<VF>       MVF;

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







//TEUCHOS_UNIT_TEST( ScalarModeField, TwoNorm_and_init ) {


  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto pc = Pimpact::createScalarField( space );
  //auto ps = Pimpact::createScalarField( space );

  //auto vel = Pimpact::create<MSF>( space );

  //double norm;
  //int N = vel->getLength();

  //// test different float values, assures that initial and norm work smoothly
  //for( double i=0.; i< 200.1; ++i ) {
    //vel->init(i/2.);
    //norm = vel->norm(Belos::TwoNorm);
    //TEST_EQUALITY( std::sqrt(std::pow(i/2.,2)*N), norm );
  //}
//}



//TEUCHOS_UNIT_TEST( ScalarModeField, dot ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto vel1 = Pimpact::create<MSF>( space );
  //auto vel2 = Pimpact::create<MSF>( space );

  //double dot;

  //TEST_EQUALITY( vel1->getLength(), vel2->getLength() )

  //int N = vel1->getLength();

  //vel1->init(0.);
  //vel2->init(1.);
  //dot = vel1->dot(*vel2);
  //TEST_EQUALITY( 0, dot );

  //vel1->init(1.);
  //vel2->init(1.);
  //dot = vel2->dot(*vel1);
  //TEST_EQUALITY( N, dot );

  //vel1->init(2.);
  //vel2->init(1.);
  //dot = vel1->dot(*vel2);
  //TEST_EQUALITY( 2*N, dot );

  //vel1->init(1.);
  //vel2->init(2.);
  //dot = vel1->dot(*vel2);
  //TEST_EQUALITY( 2*N, dot );

//}



//TEUCHOS_UNIT_TEST( ScalarModeField, scale ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto pc = Pimpact::createScalarField( space );
  //auto ps = Pimpact::createScalarField( space );

  //auto vel = Pimpact::create<MSF>( space );

  //int N = vel->getLength();
  //double norm;

  //vel->init(1.);
  //vel->scale(2.);
  //norm = vel->norm(Belos::TwoNorm);
  //TEST_EQUALITY( std::sqrt(4*N), norm)

//}



//TEUCHOS_UNIT_TEST( ScalarModeField, random ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto pc = Pimpact::createScalarField( space );
  //auto ps = Pimpact::createScalarField( space );

  //auto vel = Pimpact::create<MSF>( space );

  //int N = vel->getLength();
  //double norm;

  //vel->init(1.);
  //vel->random();
  //norm = vel->norm(Belos::TwoNorm);
  //TEST_INEQUALITY( std::sqrt(N), norm)

//}



//TEUCHOS_UNIT_TEST( ScalarModeField, add ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto p1c = Pimpact::createScalarField( space );
  //auto p1s = Pimpact::createScalarField( space );

  //auto p2c = Pimpact::createScalarField( space );
  //auto p2s = Pimpact::createScalarField( space );

  //auto p3c = Pimpact::createScalarField( space );
  //auto p3s = Pimpact::createScalarField( space );

  //auto vel1 = Pimpact::create<MSF>( space );
  //auto vel2 = Pimpact::create<MSF>( space );
  //auto vel3 = Pimpact::create<MSF>( space );

  //TEST_EQUALITY( vel1->getLength(), vel2->getLength() )
  //TEST_EQUALITY( vel2->getLength(), vel3->getLength() )
  //TEST_EQUALITY( vel1->getLength(), vel3->getLength() )

  //double norm;
  //int N = vel1->getLength();

  //vel1->init(0.);
  //vel2->init(1./2.);
  //vel3->init(1./3.);

  //vel1->add( 2., *vel2, 0., *vel3);
  //norm = vel1->norm(Belos::TwoNorm);
  //TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

  //vel1->init(0.);
  //vel2->init(1./2.);
  //vel3->init(1./3.);

  //vel1->add( 0., *vel2, 3., *vel3);
  //norm = vel1->norm(Belos::TwoNorm);
  //TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

  //vel1->init(0.);
  //vel2->init(1.);
  //vel3->init(1.);

  //vel1->add( 0.5, *vel2, 0.5, *vel3);
  //norm = vel1->norm(Belos::TwoNorm);
  //TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );
//}



//TEUCHOS_UNIT_TEST( ScalarModeField, write ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto vel = Pimpact::create<MSF>( space );

  //vel->init( 1. );
  //vel->write();

  //vel->random();
  //vel->write( 2 );

  //TEST_EQUALITY( 0, 0 )

//}



//TEUCHOS_UNIT_TEST( VectorModeField, create_init_print ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto velc = Pimpact::create<Pimpact::VectorField>(space);
  //auto vels = Pimpact::create<Pimpact::VectorField>(space);

  //auto um = Pimpact::create<MVF>( space );

  //um->init();
//}



//TEUCHOS_UNIT_TEST( VectorModeField, InfNorm_and_init ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto velc = Pimpact::create<Pimpact::VectorField>(space);
  //auto vels = Pimpact::create<Pimpact::VectorField>(space);

  //auto um = Pimpact::create<MVF>( space );

  //double norm;
  //// test different float values, assures that initial and norm work smoothly
  //for( double i=0.; i< 200.1; ++i ) {
    //um->init(i/2.);
    //norm = um->norm(Belos::InfNorm);
    //TEST_EQUALITY( i/2., norm );
  //}

  //// one test with infty-norm
  //int rank;
  //double init;
  //MPI_Comm_rank(space->comm(),&rank);
  //for( double i = 0.; i<200.1; ++i) {
    //init = 3*i-1.;
    //init = (init<0)?-init:init;
    //um->init(rank*i-1.);
    //norm = um->norm(Belos::InfNorm);
    //TEST_EQUALITY( init, norm );
  //}
//}



//TEUCHOS_UNIT_TEST( VectorModeField, TwoNorm_and_init ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto velc = Pimpact::create<Pimpact::VectorField>( space );
  //auto vels = Pimpact::create<Pimpact::VectorField>( space );

  //auto vel = Pimpact::create<MVF>( space );

  //double norm;
  //int N = vel->getLength();

  //// test different float values, assures that initial and norm work smoothly
  //for( double i=0.; i< 200.1; ++i ) {
    //vel->init(i/2.);
    //norm = vel->norm(Belos::TwoNorm);
    //TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*N), norm, eps );
  //}
//}



//TEUCHOS_UNIT_TEST( VectorModeField, dot ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );


  //auto vel1 = Pimpact::create<MVF>( space );
  //auto vel2 = Pimpact::create<MVF>( space );

  //double dot;

  //TEST_EQUALITY( vel1->getLength(), vel2->getLength() )

  //int N = vel1->getLength();

  //vel1->init(0.);
  //vel2->init(1.);
  //dot = vel1->dot(*vel2);
  //TEST_EQUALITY( 0, dot );

  //vel1->init(1.);
  //vel2->init(1.);
  //dot = vel2->dot(*vel1);
  //TEST_EQUALITY( N, dot );

  //vel1->init(2.);
  //vel2->init(1.);
  //dot = vel1->dot(*vel2);
  //TEST_EQUALITY( 2*N, dot );

  //vel1->init(1.);
  //vel2->init(2.);
  //dot = vel1->dot(*vel2);
  //TEST_EQUALITY( 2*N, dot );

//}



//TEUCHOS_UNIT_TEST( VectorModeField, scale ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto velc = Pimpact::create<Pimpact::VectorField>( space );
  //auto vels = Pimpact::create<Pimpact::VectorField>( space );

  //auto vel = Pimpact::create<MVF>( space );

  //int N = vel->getLength();
  //double norm;

  //vel->init(1.);
  //vel->scale(2.);
  //norm = vel->norm(Belos::TwoNorm);
  //TEST_FLOATING_EQUALITY( std::sqrt(4*N), norm, eps );

//}


//TEUCHOS_UNIT_TEST( VectorModeField, random ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto velc = Pimpact::create<Pimpact::VectorField>( space );
  //auto vels = Pimpact::create<Pimpact::VectorField>( space );

  //auto vel = Pimpact::create<MVF>( space );

  //int N = vel->getLength();
  //double norm;

  //vel->init(1.);
  //vel->random();
  //norm = vel->norm(Belos::TwoNorm);
  //TEST_INEQUALITY( std::sqrt(N), norm)

//}


//TEUCHOS_UNIT_TEST( VectorModeField, add ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto vel1 = Pimpact::create<MVF>( space );
  //auto vel2 = Pimpact::create<MVF>( space );
  //auto vel3 = Pimpact::create<MVF>( space );


  //TEST_EQUALITY( vel1->getLength(), vel2->getLength() )
  //TEST_EQUALITY( vel2->getLength(), vel3->getLength() )
  //TEST_EQUALITY( vel1->getLength(), vel3->getLength() )

  //double norm;
  //int N = vel1->getLength();

  //vel1->init(0.);
  //vel2->init(1./2.);
  //vel3->init(1./3.);

  //vel1->add( 2., *vel2, 0., *vel3);
  //norm = vel1->norm(Belos::TwoNorm);
  //TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

  //vel1->init(0.);
  //vel2->init(1./2.);
  //vel3->init(1./3.);

  //vel1->add( 0., *vel2, 3., *vel3);
  //norm = vel1->norm(Belos::TwoNorm);
  //TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

  //vel1->init(0.);
  //vel2->init(1.);
  //vel3->init(1.);

  //vel1->add( 0.5, *vel2, 0.5, *vel3);
  //norm = vel1->norm(Belos::TwoNorm);
  //TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

//}



//TEUCHOS_UNIT_TEST( VectorModeField, write ) {

  //auto pl = Teuchos::parameterList();

  //pl->set( "domain", domain);

  //auto space = Pimpact::createSpace( pl, false );

  //auto velc = Pimpact::create<Pimpact::VectorField>( space );
  //auto vels = Pimpact::create<Pimpact::VectorField>( space );

  //auto vel = Pimpact::create<MVF>( space );

  //vel->init( 1. );
  //vel->write();

  //vel->random();
  //vel->write( 2 );

  //TEST_EQUALITY( 0, 0 )

//}


} // end of namespace
