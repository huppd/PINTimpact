#include <iostream>
#include <cmath>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_DivOp.hpp"
#include "Pimpact_VectorField.hpp"




namespace {


using ST = double;
using OT = int;

const int sd = 3;
const int d = 3;
const int dNC = 4;

using SpaceT = Pimpact::Space<ST,OT,sd,d,dNC>;


bool testMpi = true;
ST eps = 1e-6;

int domain = 0;

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
	clp.setOption( "domain", &domain, "domain" );


	pl->set( "lx", 20. );
	pl->set( "ly", 20. );
	pl->set( "lz", 10. );


	pl->set("nx", 9 );
	pl->set("ny", 17 );
	pl->set("nz", 9 );
}





TEUCHOS_UNIT_TEST( VectorField, InfNorm_and_initvec2d ) {

	Pimpact::setBoundaryConditions( pl, domain );

	// processor grid size
	pl->set("npx", (2==SpaceT::sdim)?4:2 );
	pl->set("npy",            2 );
	pl->set("npz", (2==SpaceT::sdim)?1:2 );

	auto space = Pimpact::create<SpaceT>( pl );

	auto vel = Pimpact::create<Pimpact::VectorField>( space );

	ST norm;


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
	TEST_FLOATING_EQUALITY( ( 2==SpaceT::sdim )?0.:1., norm, eps );

}



TEUCHOS_UNIT_TEST( VectorField, initField ) {

	Pimpact::setBoundaryConditions( pl, domain );

	pl->set<ST>( "Re", 400. );

	// processor grid size
	pl->set("npx", (2==SpaceT::sdim)?4:2 );
	pl->set("npy",            2 );
	pl->set("npz", (2==SpaceT::sdim)?1:2 );

	pl->set( "Re", 400. );

	pl->set( "lx", 30. );
	pl->set( "ly", 20. );
	pl->set( "lz", 20. );

	pl->set("nx", 65 );
	pl->set("ny", 49 );
	pl->set("nz", 97 );

	pl->sublist( "Stretching in X" ).set<std::string>("Stretch Type", "cos");
	pl->sublist( "Stretching in X" ).set<ST>( "N metr L", pl->get<OT>("nx") );
	pl->sublist( "Stretching in X" ).set<ST>( "N metr U", pl->get<OT>("nx") );
	pl->sublist( "Stretching in X" ).set<ST>( "x0 L", 0.05 );
	pl->sublist( "Stretching in X" ).set<ST>( "x0 U", 0. );

	pl->sublist( "Stretching in Y" ).set<std::string>("Stretch Type", "cos");
	pl->sublist( "Stretching in Y" ).set<ST>( "N metr L", 1. );
	//pl->sublist( "Stretching in Y" ).set<ST>( "N metr L", pl->get<OT>("ny") );
	pl->sublist( "Stretching in Y" ).set<ST>( "N metr U", pl->get<OT>("ny") );
	pl->sublist( "Stretching in Y" ).set<ST>( "x0 L", 0. );
	pl->sublist( "Stretching in Y" ).set<ST>( "x0 U", 0. );

	pl->sublist( "Stretching in Z" ).set<std::string>("Stretch Type", "cos");
	pl->sublist( "Stretching in Z" ).set<ST>( "N metr L", 1. );
	pl->sublist( "Stretching in Z" ).set<ST>( "N metr U", pl->get<OT>("nz") );
	pl->sublist( "Stretching in Z" ).set<ST>( "x0 L", 0. );
	pl->sublist( "Stretching in Z" ).set<ST>( "x0 U", 0. );


	auto space = Pimpact::create<SpaceT>( pl );

	space->getInterpolateV2S()->print();
	auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto divVec = Pimpact::create<Pimpact::ScalarField>( space );

	auto divOp = Pimpact::create<Pimpact::DivOp>( space );

	//for( int i=22; i<=23; ++i ) {
	////if( 17==i )
	////vel->initField( Pimpact::EVectorField(i), 1., 1., 0.25 );
	////else
	////vel->initField( Pimpact::EVectorField(i) );
	//vel->write( i );
	//divOp->apply( *vel, *divVec );
	//auto bla = divVec->norm( Belos::InfNorm );
	//if( 0==space->rankST() )
	//std::cout << "EField: " << i << "\tmax div: " << bla << "\n";
	//bla = divVec->norm( Belos::TwoNorm );
	//if( 0==space->rankST() )
	//std::cout << "EField: " << i << "\t||div||: " << bla << "\n";
	//divVec->write( i*2 );
	//}

}


} // end of namespace
