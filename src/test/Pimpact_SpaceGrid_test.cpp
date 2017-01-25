#include <iostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_Space.hpp"




namespace {

using ST = double;
using OT = int;

const int sd = 3;
const int d = 3;
const int dNC = 4;

using Space2DT = Pimpact::Space<ST,OT,2,d,dNC>;
using Space3DT = Pimpact::Space<ST,OT,3,d,dNC>;

bool testMpi = true;
double eps = 1e+1;
int rank = 0;

int domain = 0;

ST lx = 1.;
ST ly = 1.;
ST lz = 1.;

OT nx = 33;
OT ny = 33;
OT nz = 33;
OT nf = 2;

int sx = 0;
int sy = 0;
int sz = 0;

int npx = 2;
int npy = 2;
int npz = 1;
int npf = 1;

Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();



TEUCHOS_STATIC_SETUP() {
	Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used." );
	clp.setOption( "error-tol-slack", &eps,
      "Slack off of machine epsilon used to check test results" );
	clp.setOption( "domain", &domain, "domain" );
	clp.setOption( "rank", &rank, "" );

	clp.setOption( "lx", &lx, "" );
	clp.setOption( "ly", &ly, "" );
	clp.setOption( "lz", &lz, "" );

	clp.setOption( "nx", &nx, "" );
	clp.setOption( "ny", &ny, "" );
	clp.setOption( "nz", &nz, "" );
	clp.setOption( "nf", &nf, "" );

	clp.setOption( "sx", &sx, "" );
	clp.setOption( "sy", &sy, "" );
	clp.setOption( "sz", &sz, "" );

	clp.setOption( "npx", &npx, "" );
	clp.setOption( "npy", &npy, "" );
	clp.setOption( "npz", &npz, "" );
	clp.setOption( "npf", &npf, "" );
}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST( StencilWidths, print ) {

  auto sW32 = Pimpact::createStencilWidths<d,2>(false);

  sW32->print();

  auto sW34 = Pimpact::createStencilWidths<d,4>(false);

  sW34->print();

}



TEUCHOS_UNIT_TEST( IndexSpace, localConsistency ) {

	Pimpact::setBoundaryConditions( pl, domain );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set( "nx", nx );
	pl->set( "ny", ny );
	pl->set( "nz", nz );
	pl->set( "nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
		pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
	}
	if( sy!=0 ) {
		pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	}
	if( sz!=0 ) {
		pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	}

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	auto domainSize = Pimpact::createDomainSize<ST,3>(
      pl->get<ST>("Re",1.),
      pl->get<ST>("alpha2",1.),
      pl->get<ST>("lx",2.),
      pl->get<ST>("ly",2.),
      pl->get<ST>("lz",1.),
      pl->get<ST>("origin x",0.),
      pl->get<ST>("origin y",0.),
      pl->get<ST>("origin z",0.) );

  auto stencilWidths = Pimpact::createStencilWidths<d,dNC>( false );

  auto boundaryConditionsGlobal =
		Pimpact::createBoudaryConditionsGlobal( Teuchos::rcpFromRef(pl->sublist("boundary conditions") ) );

 auto procGridSize =
	 Teuchos::tuple(
			 pl->get("npx",2),
			 pl->get("npy",2),
			 pl->get("npz",1) );

  auto gridSizeGlobal =
		Pimpact::createGridSizeGlobal<OT,3>(
				pl->get("nx",33),
				pl->get("ny",33),
				pl->get("nz",33),
				pl->get("nf",32) );

  auto procGrid =
		Pimpact::createProcGrid<OT,d>(
				procGridSize,
				boundaryConditionsGlobal );

  auto gridSizeLocal =
		Pimpact::createGridSizeLocal<OT,sd,d>(
				gridSizeGlobal,
				procGrid,
				stencilWidths );

  auto boundaryConditionsLocal =
		Pimpact::createBoudaryConditionsLocal(
				boundaryConditionsGlobal,
				procGrid );


  auto indexSpace = Pimpact::createIndexSpace<OT,d>(
			stencilWidths,
			gridSizeLocal,
			boundaryConditionsLocal,
			procGrid );

  indexSpace->print();

}



TEUCHOS_UNIT_TEST( ProcGrid, test ) {

	const int d = 4;

	Pimpact::setBoundaryConditions( pl, domain );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set( "nx", nx );
	pl->set( "ny", ny );
	pl->set( "nz", nz );
	pl->set( "nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
		pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
	}
	if( sy!=0 ) {
		pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	}
	if( sz!=0 ) {
		pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	}

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	auto domainSize = Pimpact::createDomainSize<ST,sd>(
			pl->get("Re",1.),
			pl->get("alpha2",1.),
			pl->get("lx",2.),
			pl->get("ly",2.),
			pl->get("lz",1.),
			pl->get("origin x",0.),
			pl->get("origin y",0.),
			pl->get("origin z",0.) );

	auto stencilWidths = Pimpact::createStencilWidths<d,4>( true );

	auto boundaryConditionsGlobal =
		Pimpact::createBoudaryConditionsGlobal<d>( Teuchos::rcpFromRef( pl->sublist("boundary conditions") ) );

	auto gridSizeGlobal =
		Pimpact::createGridSizeGlobal<OT,sd>(
				pl->get("nx",17),
				pl->get("ny",17),
				pl->get("nz",17),
				pl->get("nf",32) );

	auto procGridSize =
		Teuchos::tuple(
				pl->get("npx",2),
				pl->get("npy",1),
				pl->get("npz",1),
				pl->get("npf",2) );

	auto procGrid =
		Pimpact::createProcGrid<OT,d>(
				procGridSize,
				boundaryConditionsGlobal );


	std::ofstream file;
	std::string fname = "pgr.txt";
	fname.insert( 3, std::to_string( (long long)procGrid->getRank() ) );
	file.open( fname, std::ofstream::out );
	procGrid->print(file);

	auto gridSizeLocal =
		Pimpact::createGridSizeLocal<OT,sd,d>(
				gridSizeGlobal,
				procGrid,
				stencilWidths );

//	gridSizeLocal->print();

	auto boundaryConditionsLocal =
		Pimpact::createBoudaryConditionsLocal(
				boundaryConditionsGlobal,
				procGrid );


	auto indexSpace = Pimpact::createIndexSpace<OT,sd,d>(
			stencilWidths,
			gridSizeLocal,
			boundaryConditionsLocal,
			procGrid );

	file << "\n\n";
	indexSpace->print( file );
	file.close();

}


TEUCHOS_UNIT_TEST( Space, CoordinatesGlobal ) {

	Pimpact::setBoundaryConditions( pl, domain );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set( "nx", nx );
	pl->set( "ny", ny );
	pl->set( "nz", nz );
	pl->set( "nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
		pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
	}
	if( sy!=0 ) {
		pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	}
	if( sz!=0 ) {
		pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	}

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	Teuchos::RCP< const Pimpact::Space<ST,OT,sd,4,dNC> > space =
		Pimpact::create< Pimpact::Space<ST,OT,sd,4,dNC> >( pl );

	//space->print();

	auto coord = space->getCoordinatesGlobal();

	if( space->rankST()==0 )
		coord->print();

	//space->getInterpolateV2S()->print();
	auto gsg = space->getGridSizeGlobal();

	auto gridSizeGlobal =
		Pimpact::createGridSizeGlobal<OT,sd>(
				(gsg->get(0)-1)/2+1,
				(gsg->get(1)-1)/2+1,
				gsg->get(2),
				gsg->get(3)/2	);

	auto cordc =
		Pimpact::createCoordinatesGlobal(
				gridSizeGlobal,
				coord );

	//if( space->rankST()==0 )
		//cordc->print();
}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Space, create, SpaceT ) {

	Pimpact::setBoundaryConditions( pl, domain );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set( "nx", nx );
	pl->set( "ny", ny );
	pl->set( "nz", nz );
	pl->set( "nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
		pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
	}
	if( sy!=0 ) {
		pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	}
	if( sz!=0 ) {
		pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	}

	// processor grid size
	pl->set( "npx", npx );
	pl->set( "npy", npy );
	pl->set( "npz", npz );
	pl->set( "npf", npf );

	Teuchos::RCP<const SpaceT > space =
		Pimpact::create<SpaceT>( pl );

	auto coord = Pimpact::createCoordinatesLocal(
			space->getStencilWidths(),
			space->getDomainSize(),
			space->getGridSizeGlobal(),
			space->getGridSizeLocal(),
			space->getBCGlobal(),
			space->getBCLocal(),
			space->getProcGrid(),
			space->getCoordinatesGlobal() );


	if( space->getProcGrid()->getRank()==rank )
		coord->print();
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Space, create, Space2DT )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Space, create, Space3DT )

} // end of namespace
