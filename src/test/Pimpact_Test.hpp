#pragma once
#ifndef PIMPACT_TEST_HPP
#define PIMPACT_TEST_HPP

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Pimpact_Space.hpp"


template<typename ScalarT>
ScalarT order( const std::vector<ScalarT>& x, const std::vector<ScalarT>& y ) {

	TEUCHOS_TEST_FOR_EXCEPT( x.size()!=y.size() );

	const ScalarT n    = x.size();
	if( n<2 ) return( 0. );

	const ScalarT s_x  = std::accumulate(x.begin(), x.end(), 0.0);
	const ScalarT s_y  = std::accumulate(y.begin(), y.end(), 0.0);

	const ScalarT s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
	const ScalarT s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);

	const ScalarT a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);

	return( a );
}


namespace {


using ST = double;
using OT= int;

const int d = 4;
const int dNC = 4;

using D2 = Pimpact::Space<ST,OT,2,d,dNC>;
using D3 = Pimpact::Space<ST,OT,3,d,dNC>;

int print = 0;
int write = 0;
//ST eps = Teuchos::ScalarTraits<ST>::eps()*10;
ST eps = 1.e-8;

int domain = 0;

ST re     = 10.;
ST alpha2 = 10.;

ST lx = 1.;
ST ly = 1.;
ST lz = 1.;

ST x0 = 0.;
ST y0 = 0.;
ST z0 = 0.;

ST omega = 0.8;
ST winds = 1;
int sweeps = 12;
int nIter = 1;
OT ns = 4;

OT nx = 33;
OT ny = 33;
OT nz = 33;
OT nf = 1;

int sx = 0;
int sy = 0;
int sz = 0;

int npx = 1;
int npy = 1;
int npz = 1;
int npf = 1;

// MultiGrid stuff
int fs = 0;
int fe = 4;

int level = -1;


int rankbla = -1;

int maxGrids = 10;

int nMax = 1;

int rank=0; // ???


Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();



TEUCHOS_STATIC_SETUP() {

	Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
	clp.addOutputSetupOptions(true);

	clp.setOption( "eps", &eps,
			"Slack off of machine epsilon used to check test results" );
	clp.setOption( "print", &print, "" );
	clp.setOption( "write", &write, "" );
	clp.setOption( "domain", &domain, "domain" );
	clp.setOption( "omega", &omega, "" );
	clp.setOption( "wind", &winds, "" );
	clp.setOption( "sweeps", &sweeps, "" );
	clp.setOption( "nIter", &nIter, "" );
	clp.setOption( "ns", &ns, "" );

	clp.setOption( "re", &re, "" );
	clp.setOption( "alpha2", &alpha2, "" );

	clp.setOption( "lx", &lx, "" );
	clp.setOption( "ly", &ly, "" );
	clp.setOption( "lz", &lz, "" );

	clp.setOption( "x0", &x0, "" );
	clp.setOption( "y0", &y0, "" );
	clp.setOption( "z0", &z0, "" );

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

	clp.setOption( "fs", &fs, "" );
	clp.setOption( "fe", &fe, "" );
	clp.setOption( "level", &level, "" );

	clp.setOption( "rank", &rankbla, "" );
	clp.setOption( "maxGrids", &maxGrids, "" );
	clp.setOption( "nMax", &nMax, "" );

}




void setStretching() {
	// grid stretching
	
	switch( sx ) {
		case(0):
			pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "none" );
			break;
		case(1):
			pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx"))/2. );
			pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx"))/2. );
			break;
		case(2):
			pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx")) );
			pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx")) );
			break;
		case(3):
			pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(1) );
			pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(1) );
			break;
		case(4):
			pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "para" );
			pl->sublist("Stretching in X").set<ST>( "alpha", static_cast<ST>(0.1) );
			break;
	}
	switch( sy ) {
		case(0):
			pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "none" );
			break;
		case(1):
			pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny"))/2. );
			pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny"))/2. );
			break;
		case(2):
			pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny")) );
			pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny")) );
			break;
		case(3):
			pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(1) );
			pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(1) );
			break;
		case(4):
			pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "para" );
			pl->sublist("Stretching in Y").set<ST>( "alpha", static_cast<ST>(0.1) );
			break;
	}
	switch( sz ) {
		case(0):
			pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "none" );
			break;
		case(1):
			pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz"))/2. );
			pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz"))/2. );
			break;
		case(2):
			pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz")) );
			pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz")) );
			break;
		case(3):
			pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(1) );
			pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(1) );
			break;
		case(4):
			pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "para" );
			pl->sublist("Stretching in Z").set<ST>( "alpha", static_cast<ST>(0.1) );
			break;
	}
}



void setParameter( int dim ) {

	Pimpact::setBoundaryConditions( pl, domain );

	pl->set<ST>( "Re",     re     );
	pl->set<ST>( "alpha2", alpha2 );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	pl->set<ST>( "origin x", x0, "origin in x-direction" );
	pl->set<ST>( "origin y", y0, "origin in y-direction" );
	pl->set<ST>( "origin z", z0, "origin in z-direction" );

	//  grid size
	pl->set( "nx", nx );
	pl->set( "ny", ny );
	pl->set( "nz", nz );
	pl->set( "nf", nf );

	// processor grid size
	pl->set( "npx", (2==dim)?npx*npz:npx );
	pl->set( "npy", npy );
	pl->set( "npz", (2==dim)?1:npz );
	pl->set( "npf", npf );

	setStretching();
}

}

#endif // end of #ifndef PIMPACT_TEST_HPP
