#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"


#include "Pimpact_ScalarField.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_CoarsenStrategy.hpp"



namespace {


bool testMpi = true;
double errorTolSlack = 1e-6;

int domain = 1;
int ftype = 0;

typedef double S;
typedef int O;


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
      "Slack off of machine epsilon used to check test results" );
  clp.setOption(
      "ftype", &ftype,
      "Slack off of machine epsilon used to check test results" );
}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST( MultiGrid, constructor3D ) {

  typedef Pimpact::ScalarField<S,O,3> SF;
  typedef Pimpact::CoarsenStrategy<SF> CS;

  auto pl = Teuchos::parameterList();

  auto space = Pimpact::createSpace( pl );

  auto multiGrid = Pimpact::createMultiGrid<SF,CS>( space, 4 );
  std::cout << "nGridLevels: " << multiGrid->getNGrids() << "\n";
  multiGrid->print();

}


TEUCHOS_UNIT_TEST( MultiGrid, constructor4D ) {

  typedef Pimpact::ScalarField<S,O,4> SF;

  typedef Pimpact::CoarsenStrategy<SF> CS;

  auto pl = Teuchos::parameterList();

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );
  pl->set( "domain", 1 );

  pl->set( "lx", 1. );
  pl->set( "ly", 1. );
  pl->set( "lz", 1. );

  pl->set( "dim", 2 );

  pl->set("nx", 129 );
  pl->set("ny", 65 );
  pl->set("nz", 2 );

  pl->set("nf", 256 );
  //  pl->set("nfs", 0 );
  //  pl->set("nfe", 0 );

  // processor grid size
  pl->set("npx", 2 );
  pl->set("npy", 1 );
  pl->set("npz", 1 );
  pl->set("npf", 2 );


  auto space = Pimpact::createSpace<S,O,4>( pl );

  space->print();

  auto multiGrid = Pimpact::createMultiGrid<SF,CS>( space, 2 );
  std::cout << "nGridLevels: " << multiGrid->getNGrids() << "\n";
  multiGrid->print();


}


TEUCHOS_UNIT_TEST( MultiGrid, Restrictor3D ) {

  typedef Pimpact::ScalarField<S,O,3> SF;

  typedef Pimpact::CoarsenStrategy<SF> CS;

  auto pl = Teuchos::parameterList();

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );
  pl->set( "domain", domain );

  pl->set( "lx", 1. );
  pl->set( "ly", 1. );
  pl->set( "lz", 1. );

  pl->set( "dim", 2 );

  pl->set("nx", 129 );
  pl->set("ny", 129 );
  pl->set("nz", 2 );

  pl->set("nf", 0 );
  //  pl->set("nfs", 0 );
  //  pl->set("nfe", 0 );

  // processor grid size
  pl->set("npx", 2 );
  pl->set("npy", 2 );
  pl->set("npz", 1 );
  pl->set("npf", 1 );

  auto space = Pimpact::createSpace<S,O,3>( pl );

  space->print();

  Pimpact::EField type[] = {Pimpact::EField::S, Pimpact::EField::U, Pimpact::EField::V };

  for( int i=0; i<3; ++i ) {
    //    {
    //      int i = ftype;
    std::cout << "type: " << i << "\n";

    auto multiGrid = Pimpact::createMultiGrid<SF,CS>( space, 2, type[i] );

    auto fieldf = multiGrid->getField( 0 );
    auto fieldc = multiGrid->getField( 1 );

    auto op = multiGrid->getRestrictionOp( 0 );

    // the zero test
    fieldf->init( 0. );
    fieldc->init( 1. );

    op->apply( *fieldf, *fieldc );

    TEST_FLOATING_EQUALITY( 0., fieldf->norm(), errorTolSlack );
    TEST_FLOATING_EQUALITY( 0., fieldc->norm(), errorTolSlack );

    // the random test
    fieldf->random();

    TEST_INEQUALITY( 0., fieldf->norm() );

    op->apply( *fieldf, *fieldc );

    TEST_INEQUALITY( 0., fieldc->norm() );


    //      op->print();
    //      fieldc->print();

    // the init test
    fieldf->initField( Pimpact::ConstField,1. );

    TEST_FLOATING_EQUALITY( 1., fieldf->norm(Belos::InfNorm), errorTolSlack );

    TEST_FLOATING_EQUALITY( (S)fieldf->getLength(), fieldf->norm(Belos::OneNorm), errorTolSlack  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldf->getLength() ), fieldf->norm(Belos::TwoNorm), errorTolSlack  );

    fieldc->init(0.);

    op->apply( *fieldf, *fieldc );

    fieldf->write( 0 );
    fieldc->write( 1 );

    TEST_FLOATING_EQUALITY( 1., fieldc->norm(Belos::InfNorm), errorTolSlack );

    TEST_FLOATING_EQUALITY( (S)fieldc->getLength(), fieldc->norm(Belos::OneNorm), errorTolSlack  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldc->getLength() ), fieldc->norm(Belos::TwoNorm), errorTolSlack  );

  }

}


TEUCHOS_UNIT_TEST( MultiGrid, Interpolator3D ) {

  typedef Pimpact::ScalarField<S,O,3> SF;

  typedef Pimpact::CoarsenStrategy<SF> CS;

  auto pl = Teuchos::parameterList();

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );
  pl->set( "domain", domain );

  pl->set( "lx", 1. );
  pl->set( "ly", 1. );
  pl->set( "lz", 1. );

  pl->set( "dim", 2 );

  pl->set("nx", 33 );
  pl->set("ny", 33 );
  //    pl->set("nx", 17 );
  //    pl->set("ny", 17 );
  pl->set("nz", 2 );

  pl->set("nf", 0 );
  //  pl->set("nfs", 0 );
  //  pl->set("nfe", 0 );

  // processor grid size
  pl->set("npx", 2 );
  pl->set("npy", 2 );
  //    pl->set("npx", 1 );
  //    pl->set("npy", 1 );
  pl->set("npz", 1 );
  pl->set("npf", 1 );

  auto space = Pimpact::createSpace<S,O,3>( pl );

  //    space->print();

  Pimpact::EField type[] = { Pimpact::EField::S, Pimpact::EField::U, Pimpact::EField::V };

  for( int i=1; i<2; ++i ) {
    //    {
    //      int i = ftype;
    std::cout << "type: " << i << "\n";

    auto multiGrid = Pimpact::createMultiGrid<SF,CS>( space, 2, type[i] );

    multiGrid->getSpace(0)->print();
    multiGrid->getSpace(1)->print();

    auto fieldf = multiGrid->getField( 0 );
    auto fieldc = multiGrid->getField( 1 );

    auto op = multiGrid->getInterpolationOp( 0 );
    //      auto op = Pimpact::createInterpolationOp( multiGrid->getSpace(1), multiGrid->getSpace(0) );

    op->print();
    // the zero test

    fieldf->init( 1. );

    fieldc->initField( Pimpact::ConstField, 0. );

    op->apply( *fieldc, *fieldf );

    TEST_FLOATING_EQUALITY( 0., fieldf->norm(), errorTolSlack );
    TEST_FLOATING_EQUALITY( 0., fieldc->norm(), errorTolSlack );


    //      fieldc->print();
    fieldf->print();


    // the random test
    fieldc->random();
    fieldf->init(0.);

    TEST_INEQUALITY( 0., fieldc->norm() );

    op->apply( *fieldc, *fieldf );

    TEST_INEQUALITY( 0., fieldc->norm() );


    // the stronger init test
    fieldc->initField( Pimpact::ConstField, 1. );

    TEST_FLOATING_EQUALITY( 1., fieldc->norm(Belos::InfNorm), errorTolSlack );

    TEST_FLOATING_EQUALITY( (S)fieldc->getLength(), fieldc->norm(Belos::OneNorm), errorTolSlack  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldc->getLength() ), fieldc->norm(Belos::TwoNorm), errorTolSlack  );

    fieldf->init(0.);

    op->apply( *fieldc, *fieldf );

    fieldf->write(0);
    fieldc->write(1);

    TEST_FLOATING_EQUALITY( 1., fieldf->norm(Belos::InfNorm), errorTolSlack );

    TEST_FLOATING_EQUALITY( (S)fieldf->getLength(), fieldf->norm(Belos::OneNorm), errorTolSlack  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldf->getLength() ), fieldf->norm(Belos::TwoNorm), errorTolSlack  );

  }

}

} // end of namespace
