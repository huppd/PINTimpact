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


bool testMpi = true;
double eps = 1e-6;

int dim = 3;
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
	clp.setOption(
			"domain", &domain,
			"domain" );
	clp.setOption(
			"dim", &dim,
			"dim" );

	pl->set( "dim", dim );
	pl->set( "domain", domain );

	pl->set( "lx", 2. );
	pl->set( "ly", 2. );
	pl->set( "lz", 1. );


	pl->set("nx", 25 );
	pl->set("ny", 17 );
	if(  3==dim )
		pl->set("nz", 9 );

	// processor grid size
	pl->set("npx", 2 );
	pl->set("npy", 2 );
	if(  3==dim )
		pl->set("npz", 2 );

}





TEUCHOS_UNIT_TEST( VectorField, InfNorm_and_initvec2d ) {

	pl->set( "dim", dim );
	pl->set( "domain", domain );

	// processor grid size
	pl->set("npx", (2==dim)?4:2 );
	pl->set("npy",            2 );
	pl->set("npz", (2==dim)?1:2 );

	auto space = Pimpact::createSpace( pl );

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

	pl->set( "dim", dim );
	pl->set( "domain", domain );

	pl->set<double>( "Re", 400. );

	// processor grid size
	pl->set("npx", (2==dim)?4:2 );
	pl->set("npy",            2 );
	pl->set("npz", (2==dim)?1:2 );
	pl->set("npx", 1 );
	pl->set("npy", 1 );
	pl->set("npz", 1 );

	pl->set( "Re", 400. );

	pl->set( "lx", 30. );
	pl->set( "ly", 20. );
	pl->set( "lz", 20. );

	pl->set("nx", 65 );
	pl->set("ny", 49 );
	pl->set("nz", 97 );

	pl->set("nx", 1025 );
	pl->set("ny", 7 );
	pl->set("nz", 7 );

	pl->sublist( "Stretching in X" ).set<std::string>("Stretch Type", "Parabola");
	pl->sublist( "Stretching in X" ).set<double>( "N metr L", pl->get<int>("nx") );
	pl->sublist( "Stretching in X" ).set<double>( "N metr U", pl->get<int>("nx") );
	pl->sublist( "Stretching in X" ).set<double>( "x0 L", 0.05 );
	pl->sublist( "Stretching in X" ).set<double>( "x0 U", 0. );

	pl->sublist( "Stretching in Y" ).set<std::string>("Stretch Type", "Parabola");
	pl->sublist( "Stretching in Y" ).set<double>( "N metr L", 1. );
	//pl->sublist( "Stretching in Y" ).set<double>( "N metr L", pl->get<int>("ny") );
	pl->sublist( "Stretching in Y" ).set<double>( "N metr U", pl->get<int>("ny") );
	pl->sublist( "Stretching in Y" ).set<double>( "x0 L", 0. );
	pl->sublist( "Stretching in Y" ).set<double>( "x0 U", 0. );

	pl->sublist( "Stretching in Z" ).set<std::string>("Stretch Type", "Parabola");
	pl->sublist( "Stretching in Z" ).set<double>( "N metr L", 1. );
	pl->sublist( "Stretching in Z" ).set<double>( "N metr U", pl->get<int>("nz") );
	pl->sublist( "Stretching in Z" ).set<double>( "x0 L", 0. );
	pl->sublist( "Stretching in Z" ).set<double>( "x0 U", 0. );


	auto space = Pimpact::createSpace( pl );

	space->getInterpolateV2S()->print();
	auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto divVec = Pimpact::create<Pimpact::ScalarField>( space );

	auto divOp = Pimpact::create<Pimpact::DivOp>( space );

	for( int i=22; i<=23; ++i ) {
		if( 17==i )
			vel->initField( Pimpact::EVectorField(i), 1., 1., 0.25 );
		else
			vel->initField( Pimpact::EVectorField(i) );
		vel->write( i );
		divOp->apply( *vel, *divVec );
		auto bla = divVec->norm( Belos::InfNorm );
		if( 0==space->rankST() )
			std::cout << "EField: " << i << "\tmax div: " << bla << "\n";
		bla = divVec->norm( Belos::TwoNorm );
		if( 0==space->rankST() )
			std::cout << "EField: " << i << "\t||div||: " << bla << "\n";
		divVec->write( i*2 );
	}

}


} // end of namespace
