#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_CoarsenStrategy.hpp"




namespace {

typedef double S;
typedef int O;

typedef Pimpact::Space<S,O,3,4> Space3T;
typedef Pimpact::Space<S,O,4,4> Space4T;

bool testMpi = true;
double eps = 1e-6;

int domain = 1;
int ftype = 0;



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
      "Slack off of machine epsilon used to check test results" );
  clp.setOption(
      "ftype", &ftype,
      "Slack off of machine epsilon used to check test results" );
}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST( MultiGrid, constructor3D ) {

  typedef Pimpact::ScalarField<Space3T> SF;
  typedef Pimpact::CoarsenStrategy<Space3T> CS;

  auto pl = Teuchos::parameterList();

  auto space = Pimpact::createSpace( pl );

  auto multiGrid = Pimpact::createMultiGrid<SF,SF,CS>( space, 4 );
  std::cout << "nGridLevels: " << multiGrid->getNGrids() << "\n";
  if( space->rankST()==0 )
    multiGrid->print();

}



TEUCHOS_UNIT_TEST( MultiGrid, constructor4D ) {

  typedef Pimpact::ScalarField<Space4T> SF;
  typedef Pimpact::CoarsenStrategy<Space4T> CS;

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

  auto multiGrid = Pimpact::createMultiGrid<SF,SF,CS>( space, 2 );
  std::cout << "nGridLevels: " << multiGrid->getNGrids() << "\n";
  multiGrid->print();

}



TEUCHOS_UNIT_TEST( MultiGrid, Restrictor3D ) {

  typedef Pimpact::ScalarField<Space3T> SF;
  typedef Pimpact::CoarsenStrategy<Space3T> CS;

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
    std::cout << "type: " << i << "\n";

    auto multiGrid = Pimpact::createMultiGrid<SF,SF,CS>( space, 2, type[i] );

    auto fieldf = multiGrid->getField( 0 );
    auto fieldc = multiGrid->getField( 1 );

    auto op = multiGrid->getRestrictionOp( 0 );

    // the zero test
    fieldf->init( 0. );
    fieldc->init( 1. );

    op->apply( *fieldf, *fieldc );

    TEST_FLOATING_EQUALITY( 0., fieldf->norm(), eps );
    TEST_FLOATING_EQUALITY( 0., fieldc->norm(), eps );

    // the random test
    fieldf->random();

    TEST_INEQUALITY( 0., fieldf->norm() );

    op->apply( *fieldf, *fieldc );

    TEST_INEQUALITY( 0., fieldc->norm() );

    // the strong test
    fieldf->initField( Pimpact::ConstField,1. );

    TEST_FLOATING_EQUALITY( 1., fieldf->norm(Belos::InfNorm), eps );

    TEST_FLOATING_EQUALITY( (S)fieldf->getLength(), fieldf->norm(Belos::OneNorm), eps  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldf->getLength() ), fieldf->norm(Belos::TwoNorm), eps  );

    fieldc->init(0.);

    op->apply( *fieldf, *fieldc );

    fieldf->write( 0 );
    fieldc->write( 1 );

    TEST_FLOATING_EQUALITY( 1., fieldc->norm(Belos::InfNorm), eps );

    TEST_FLOATING_EQUALITY( (S)fieldc->getLength(), fieldc->norm(Belos::OneNorm), eps  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldc->getLength() ), fieldc->norm(Belos::TwoNorm), eps  );

    // the hard test
    fieldf->initField( Pimpact::Grad2D_inX,1. );

    fieldc->init(0.);

    op->apply( *fieldf, *fieldc );

    fieldf->write( 0 );
    fieldc->write( 1 );

//    TEST_FLOATING_EQUALITY( 1., fieldc->norm(Belos::InfNorm), eps );
//
//    TEST_FLOATING_EQUALITY( (S)fieldc->getLength(), fieldc->norm(Belos::OneNorm), eps  );
//
//    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldc->getLength() ), fieldc->norm(Belos::TwoNorm), eps  );


  }

}



TEUCHOS_UNIT_TEST( MultiGrid, Interpolator3D ) {

  typedef Pimpact::ScalarField<Space3T> SF;
  typedef Pimpact::CoarsenStrategy<Space3T> CS;

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
//  pl->set("nx", 17 );
//  pl->set("ny", 17 );
//  pl->set("nx", 9 );
//  pl->set("ny", 9);
  pl->set("nz", 2 );

  pl->set("nf", 0 );
  //  pl->set("nfs", 0 );
  //  pl->set("nfe", 0 );

  // processor grid size
  pl->set("npx", 2 );
  pl->set("npy", 2 );
//  pl->set("npx", 1 );
//  pl->set("npy", 1 );
  pl->set("npz", 1 );
  pl->set("npf", 1 );

  auto space = Pimpact::createSpace<S,O,3>( pl );

  //    space->print();

  Pimpact::EField type[] = { Pimpact::EField::S, Pimpact::EField::U, Pimpact::EField::V };

  for( int i=0; i<1; ++i ) {
    //    {
    //      int i = ftype;
    std::cout << "type: " << i << "\n";

    auto multiGrid = Pimpact::createMultiGrid<SF,SF,CS>( space, 2, type[i] );

//    multiGrid->getSpace(0)->print();
//    multiGrid->getSpace(1)->print();

    auto fieldf = multiGrid->getField( 0 );
    auto fieldc = multiGrid->getField( 1 );

    auto op = multiGrid->getInterpolationOp( 0 );
    //      auto op = Pimpact::createInterpolationOp( multiGrid->getSpace(1), multiGrid->getSpace(0) );

    if( space->rankST()==0 )
      op->print();

    // the zero test

    fieldf->init( 1. );

    fieldc->initField( Pimpact::ConstField, 0. );

    op->apply( *fieldc, *fieldf );

//    fieldf->write(0);
//    fieldc->write(1);

    TEST_FLOATING_EQUALITY( 0., fieldf->norm(), eps );
    TEST_FLOATING_EQUALITY( 0., fieldc->norm(), eps );

    //   fieldc->print();
    //   fieldf->print();


    // the random test
    fieldc->random();
    fieldf->init(0.);

    TEST_INEQUALITY( 0., fieldc->norm() );

    op->apply( *fieldc, *fieldf );

    fieldf->write(0);
    fieldc->write(1);

    TEST_INEQUALITY( 0., fieldc->norm() );


    // the stronger init test
    fieldc->initField( Pimpact::ConstField, 1. );
    fieldf->init(0.);

    op->apply( *fieldc, *fieldf );

    fieldf->write(2);
    fieldc->write(3);

    TEST_FLOATING_EQUALITY( 1., fieldc->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( (S)fieldc->getLength(), fieldc->norm(Belos::OneNorm), eps  );
    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldc->getLength() ), fieldc->norm(Belos::TwoNorm), eps  );

    TEST_FLOATING_EQUALITY( 1., fieldf->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( (S)fieldf->getLength(), fieldf->norm(Belos::OneNorm), eps  );
    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldf->getLength() ), fieldf->norm(Belos::TwoNorm), eps  );

    // hardcore test init test in X
    fieldc->initField( Pimpact::Grad2D_inX );
    fieldf->init(0.);
    auto sol = fieldf->clone();
    sol->initField( Pimpact::Grad2D_inX );

    op->apply( *fieldc, *fieldf );
    fieldf->write(4);
    fieldc->write(5);


    sol->add( 1., *sol, -1., *fieldf );
    sol->write(90);

    std::cout << "error GradX: " << fieldf->norm() << "\n";

    TEST_FLOATING_EQUALITY( fieldf->norm()/std::sqrt( (S)fieldf->getLength() ),
                            fieldc->norm()/std::sqrt( (S)fieldc->getLength() ), eps*1e5  );



    // hardcore test init test in X
    fieldc->initField( Pimpact::Grad2D_inY );
    fieldf->init(0.);
    sol->initField( Pimpact::Grad2D_inY );
    sol->write(81);

    op->apply( *fieldc, *fieldf );

    fieldf->write(6);
    fieldc->write(7);

    sol->add( 1., *sol, -1., *fieldf );
    sol->write(91);

    TEST_FLOATING_EQUALITY( fieldf->norm()/std::sqrt( (S)fieldf->getLength() ),
                            fieldc->norm()/std::sqrt( (S)fieldc->getLength() ), eps*1e5  );

    std::cout << "error GradY: " << fieldf->norm() << "\n";

  }

}

} // end of namespace
