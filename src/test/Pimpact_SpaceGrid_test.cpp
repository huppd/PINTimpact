#include <iostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_Space.hpp"



namespace {

bool testMpi = true;
double eps = 1e+1;

int dim = 3;

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
TEUCHOS_UNIT_TEST( StencilWidths, print ) {

  auto sW32 = Pimpact::createStencilWidths<3,2>(false);

  sW32->print();

  auto sW34 = Pimpact::createStencilWidths<3,4>(false);

  sW34->print();

}



TEUCHOS_UNIT_TEST( IndexSpace, localConsistency ) {

  const int d = 3;

  auto pl = Teuchos::parameterList();
  auto domainSize = Pimpact::createDomainSize(
      pl->get( "dim", dim ),
      pl->get("Re",1.),
      pl->get("alpha2",1.),
      pl->get("lx",2.),
      pl->get("ly",2.),
      pl->get("lz",1.) );

  auto stencilWidths = Pimpact::createStencilWidths<d,4>( false );

  auto boundaryConditionsGlobal =
		Pimpact::createBoudaryConditionsGlobal(
				Pimpact::EDomainType( pl->get("domain",2) ) );

 auto procGridSize =
	 Teuchos::tuple(
			 pl->get("npx",2),
			 pl->get("npy",2),
			 pl->get("npz",1) );

  auto gridSizeGlobal =
		Pimpact::createGridSizeGlobal<O,d>(
				pl->get("nx",33),
				pl->get("ny",33),
				pl->get("nz",33),
				pl->get("nf",32) );

  auto procGrid =
		Pimpact::createProcGrid<O,d>(
				procGridSize,
				boundaryConditionsGlobal );

  auto gridSizeLocal =
		Pimpact::createGridSizeLocal<O,d>(
				gridSizeGlobal,
				procGrid,
				stencilWidths );

  auto boundaryConditionsLocal =
		Pimpact::createBoudaryConditionsLocal(
				boundaryConditionsGlobal,
				procGrid );


  auto indexSpace = Pimpact::createIndexSpace<O,d>(
			stencilWidths,
			gridSizeLocal,
			boundaryConditionsLocal,
			procGrid );

  indexSpace->print();

}


TEUCHOS_UNIT_TEST( ProcGrid, test ) {

	const int d = 4;

	auto pl = Teuchos::parameterList();
	auto domainSize = Pimpact::createDomainSize(
			pl->get( "dim", dim ),
			pl->get("Re",1.),
			pl->get("alpha2",1.),
			pl->get("lx",2.),
			pl->get("ly",2.),
			pl->get("lz",1.) );

	auto stencilWidths = Pimpact::createStencilWidths<d,4>( true );

	auto boundaryConditionsGlobal =
		Pimpact::createBoudaryConditionsGlobal<d>( Pimpact::EDomainType( pl->get("domain",2) ) );

	auto gridSizeGlobal =
		Pimpact::createGridSizeGlobal<O,d>(
				pl->get("nx",17),
				pl->get("ny",17),
				pl->get("nz",17),
				pl->get("nf",32) );

	auto procGridSize =
		Teuchos::tuple(
				pl->get("npx",1),
				pl->get("npy",1),
				pl->get("npz",1),
				pl->get("npf",1) );

	auto procGrid =
		Pimpact::createProcGrid<O,d>(
				procGridSize,
				boundaryConditionsGlobal );


	std::ofstream file;
	std::string fname = "pgr.txt";
	fname.insert( 3, std::to_string( (long long)procGrid->getRank() ) );
	file.open( fname, std::ofstream::out );
	procGrid->print(file);

	auto gridSizeLocal =
		Pimpact::createGridSizeLocal<O,d>(
				gridSizeGlobal,
				procGrid,
				stencilWidths );

//	gridSizeLocal->print();

	auto boundaryConditionsLocal =
		Pimpact::createBoudaryConditionsLocal(
				boundaryConditionsGlobal,
				procGrid );


	auto indexSpace = Pimpact::createIndexSpace<O,d>(
			stencilWidths,
			gridSizeLocal,
			boundaryConditionsLocal,
			procGrid );

	file << "\n\n";
	indexSpace->print( file );
	file.close();

}


TEUCHOS_UNIT_TEST( Space, GlobalGridCoordinates ) {

  auto space = Pimpact::createSpace();

  auto coord = space->getCoordinatesGlobal();

  coord->print();

}


TEUCHOS_UNIT_TEST( Space, LocalGridCoordinates ) {

  auto space = Pimpact::createSpace();

  auto coord = Pimpact::createGridCoordinatesLocal(
      space->getStencilWidths(),
      space->getDomainSize(),
      space->getGridSizeGlobal(),
      space->getGridSizeLocal(),
      space->getBCGlobal(),
      space->getBCLocal(),
      space->getProcGrid(),
      space->getCoordinatesGlobal() );


  if( space->getProcGrid()->getRank()==3 )
    coord->print();

}

} // end of namespace
