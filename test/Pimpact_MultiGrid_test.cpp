// Pimpact_SalarVectorSpace_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
//#include <Teuchos_Array.hpp>
//#include <Teuchos_Tuple.hpp>
//#include "Teuchos_Range1D.hpp"
//#include <Teuchos_CommHelpers.hpp>

#include "pimpact.hpp"
#include "Pimpact_Space.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
//#include "Pimpact_FieldSpace.hpp"
//#include "Pimpact_IndexSpace.hpp"
//#include "Pimpact_ScalarField.hpp"
//#include "Pimpact_VectorField.hpp"
//#include "Pimpact_ModeField.hpp"
//#include "Pimpact_MultiField.hpp"

#include <iostream>




namespace {


bool testMpi = true;
double errorTolSlack = 1e-6;

bool isImpactInit = false;

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
}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST( MultiGrid, constructor3D ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  typedef Pimpact::ScalarField<S,O,3> SF;

  typedef Pimpact::CoarsenStrategy<SF> CS;

  auto space = Pimpact::createSpace();
  space->print();

  auto asdf = Pimpact::createMultiGrid<SF,CS>( space, 10 );
  std::cout << "nGridLevels: " << asdf->getNGrids() << "\n";

}


TEUCHOS_UNIT_TEST( MultiGrid, constructor4D ) {
  if( !isImpactInit ) {

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

    auto asdf = Pimpact::createMultiGrid<SF,CS>( space, 10 );
    std::cout << "nGridLevels: " << asdf->getNGrids() << "\n";

  }

}


} // end of namespace
