#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Pimpact_Space.hpp"

#include "Pimpact_Domain.hpp"
#include "Pimpact_ProcGrid.hpp"



namespace {

bool testMpi = true;
double eps = 1e+1;

bool isImpactInit = false;
int dim = 2;

typedef int O;
typedef double S;


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
      "dim", &dim,
      "dim" );
}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST( StencilWidths, local_consistency ) {


  auto sW32 = Pimpact::createStencilWidths<3,2>();

  sW32->print();

  auto sW34 = Pimpact::createStencilWidths<3,4>();

  sW34->print();

}



TEUCHOS_UNIT_TEST( IndexSpace, local_consistency ) {

  const int d = 3;

  auto pl = Teuchos::parameterList();
  auto domainSize = Pimpact::createDomainSize(
      pl->get( "dim", dim ),
      pl->get("Re",1.),
      pl->get("alpha2",1.),
      pl->get("lx",2.),
      pl->get("ly",2.),
      pl->get("lz",1.) );

  auto stencilWidths = Pimpact::createStencilWidths<d,4>();

  auto boundaryConditionsGlobal = Pimpact::createBoudaryConditionsGlobal( Pimpact::EDomainType( pl->get("domain",2) ) );

  auto procGridSize = Pimpact::createProcGridSize<O,d>( pl->get("npx",2), pl->get("npy",2), pl->get("npz",1), pl->get("npf",1) );

  auto gridSizeGlobal = Pimpact::createGridSizeGlobal<O,d>( pl->get("nx",33), pl->get("ny",33), pl->get("nz",2), pl->get("nf",32) );

  auto gridSizeLocal = Pimpact::createGridSizeLocal<O,d>( gridSizeGlobal, procGridSize, stencilWidths );

  auto procGrid = Pimpact::createProcGrid<O,d>( gridSizeLocal, boundaryConditionsGlobal, procGridSize );

  auto boundaryConditionsLocal = Pimpact::createBoudaryConditionsLocal( boundaryConditionsGlobal, procGridSize, procGrid );

  auto fieldSpace = Pimpact::createStencilWidths<d,4>();

  auto indexSpace = Pimpact::createIndexSpace<O,d>( fieldSpace, gridSizeLocal, boundaryConditionsLocal, false );

  indexSpace->print();

}



TEUCHOS_UNIT_TEST( Space, GlobalGridCoordinates ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }
  auto space = Pimpact::createSpace();

  auto coord = space->getCoordinatesGlobal();

  coord->print();

}


TEUCHOS_UNIT_TEST( Space, LocalGridCoordinates ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }
  auto space = Pimpact::createSpace();

  auto coord = Pimpact::createGridCoordinatesLocal(
      space->getStencilWidths(),
      space->getDomain()->getDomainSize(),
      space->getGridSizeGlobal(),
      space->getGridSizeLocal(),
      space->getDomain()->getBCGlobal(),
      space->getDomain()->getBCLocal(),
      space->getProcGrid(),
      space->getCoordinatesGlobal() );


  if( space->getProcGrid()->getRank()==3 )
    coord->print();

}

} // end of namespace
