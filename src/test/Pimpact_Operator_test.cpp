#include <functional>
#include <iostream>
#include <cmath>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Fields.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_TransferOp.hpp"
#include "Pimpact_VectorFieldOpWrap.hpp"



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
const int d = 3;
const int dNC = 4;
//const int dNC = 2;

using SpaceT = Pimpact::Space<ST,OT,d,dNC>;

using SF = typename Pimpact::ScalarField<SpaceT>;
using VF = typename Pimpact::VectorField<SpaceT>;
using MSF = typename Pimpact::ModeField<SF>;
using MVF = typename Pimpact::ModeField<VF>;



bool testMpi = true;
int print = 0;
//ST eps = Teuchos::ScalarTraits<ST>::eps()*10;
ST eps = 1.e-8;

int dim = 3;
int domain = 0;

ST lx = 1.;
ST ly = 1.;
ST lz = 1.;

ST omega = 0.8;
ST winds = 1;
int sweeps = 12;
int nIter = 1;
OT ns = 8;

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

Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

//Teuchos::RCP< std::ostream> outp;
//std::ostream& out;

int rank=0;



TEUCHOS_STATIC_SETUP() {

  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used." );
  clp.setOption(
      "eps", &eps,
      "Slack off of machine epsilon used to check test results" );
	clp.setOption( "print", &print, "" );
	clp.setOption( "domain", &domain, "domain" );
	clp.setOption( "dim", &dim, "dim" );
	clp.setOption( "omega", &omega,
      "Slack off of machine epsilon used to check test results" );
	clp.setOption( "wind", &winds,
      "Slack off of machine epsilon used to check test results" );
	clp.setOption( "sweeps", &sweeps, "" );
	clp.setOption( "nIter", &nIter, "" );
	clp.setOption( "ns", &ns, "" );

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

	//int rank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//if( 0==rank )
		//outp =  Teuchos::rcpFromRef( std::cout );
	//else
		//outp = Teuchos::rcp(new Teuchos::oblackholestream() );
	//out = *outp;

}



TEUCHOS_UNIT_TEST( BasicOperator, DivOp ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

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
		//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
		//pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
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

  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	auto p   = Pimpact::create<Pimpact::ScalarField>( space );
	auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto sol = p->clone( Pimpact::ECopyType::Shallow );

	// zero test
	vel->initField();

	auto op = Pimpact::create<Pimpact::DivOp>( space );
	std::cout << "\n\nhello\n\n";

	if( print )
		op->print();

	op->apply( *vel, *p );

	TEST_EQUALITY( p->norm( Belos::InfNorm )<eps, true );


	// random test
	vel->random();

	op->apply(*vel,*p);

	TEST_EQUALITY( p->norm( Belos::InfNorm )>eps, true );

	// circle XY test
	vel->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY, -1. );
	vel->getFieldPtr( Pimpact::V )->initField(  );
	vel->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX, 1. );

	op->apply(*vel,*p);

	TEST_EQUALITY( p->norm( Belos::InfNorm )<eps, true );

	// circle XZ test
	vel->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY, -1. );
	vel->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX, 1. );
	vel->getFieldPtr( Pimpact::W )->initField(  );

	op->apply( *vel, *p );

	TEST_EQUALITY( p->norm( Belos::InfNorm )<eps, true );

	// Grad test
	vel->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
	vel->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY, -2. );
	vel->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );

	p->random();

	op->apply(*vel,*p);

	sol->initField( Pimpact::ConstField, space->dim() );
	sol->add( 1., *sol, -1., *p );

	ST error = sol->norm( Belos::InfNorm );
	if( 0==rank )	
		std::cout << "GradXYZ error: " << error << "\n";
	TEST_EQUALITY( error>eps, true );

	// Grad test
	vel->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inZ );
	vel->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX );
	vel->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inY );

	p->random();

	op->apply(*vel,*p);

	error = p->norm( Belos::InfNorm );
	if( 0==rank )	
		std::cout << "GradZXY error: " << error << "\n";
	TEST_EQUALITY( error<eps, true );

	// Grad test
	vel->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY );
	vel->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inZ );
	vel->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX );

	p->random();

	op->apply(*vel,*p);

	error = p->norm( Belos::InfNorm );
	std::cout << "erro: " << error << "\n";
	TEST_EQUALITY( error<eps, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, InterpolateV2SOp ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  auto p = Pimpact::create<Pimpact::ScalarField>( space );
	auto sol = p->clone();
  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto op = Pimpact::createInterpolateV2S( space );
	
	if( print )
		op->print();

	if( 0==rank )
		std::cout << "\n";
	for( int i=0; i<space->dim(); ++i ) {

		for( int bla=0; bla<=((domain!=1)?6:0); ++bla ) {

			Pimpact::EScalarField type =	static_cast<Pimpact::EScalarField>(bla);

			vel->getFieldPtr( i )->initField( type );
			sol->initField( type );

			p->random();
			op->apply( vel->getConstField(i), *p );
			sol->add( 1., *sol, -1., *p );

			ST error = sol->norm( Belos::InfNorm );
			if( 0==rank )
				std::cout << "field " << Pimpact::toString( static_cast<Pimpact::EField>(i) ) << ", error("<<Pimpact::toString(type)<<"): " << error << "\n";
			TEST_EQUALITY( error<eps, true );
			if( error>= eps ){
				std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
				sol->print(
						*Pimpact::createOstream( "error_int"+Pimpact::toString( static_cast<Pimpact::EField>(i) )+"2S_"+Pimpact::toString(type)+"_r"+r+".txt" ));
			}
		}
		
	}

}



TEUCHOS_UNIT_TEST( BasicOperator, InterpolateS2VOp ) {

	Pimpact::setBoundaryConditions( pl, domain );
	pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

	auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	auto p = Pimpact::create<Pimpact::ScalarField>( space );
	auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto solv = vel->clone();

	auto op = Pimpact::create<Pimpact::InterpolateS2V>( space );

	if( print )
		op->print();

	if( 0==rank )
		std::cout << "\n";

	for( int i=0; i<space->dim(); ++i ) {

		auto sol = solv->getFieldPtr( i );

		for( int bla=0; bla<=(domain!=1?6:0); ++bla ) {

			Pimpact::EScalarField type =	static_cast<Pimpact::EScalarField>(bla);

			vel->getFieldPtr( i )->random();
			sol->initField( type );
			p->initField( type );

			op->apply( *p, vel->getField( i ) );
			sol->add( 1., *sol, -1., vel->getConstField(i) );

			ST error = sol->norm( Belos::InfNorm );
			if( 0==rank )
				std::cout << "field "
					<< Pimpact::toString( static_cast<Pimpact::EField>(i) )
					<< ", error("<<Pimpact::toString(type)<<"): "
					<< error << "\n";

			TEST_EQUALITY( error<eps, true );
			if( error>=eps ) {
				std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
				sol->print(
						*Pimpact::createOstream( "error_int_S2"+Pimpact::toString( static_cast<Pimpact::EField>(i) )+"_"+Pimpact::toString(type)+"_r"+r+".txt" ));
			}

		}

	}

}



TEUCHOS_UNIT_TEST( BasicOperator, TransferOp ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

  using FSpaceT = Pimpact::Space<ST,OT,d,4>;
  using CSpaceT = Pimpact::Space<ST,OT,d,2>;

  auto fSpace = Pimpact::createSpace<ST,OT,d,4>( pl );
  auto cSpace = Pimpact::createSpace<ST,OT,d,2>( pl );

  auto fx = Pimpact::create<Pimpact::ScalarField>( fSpace );
  auto cx = Pimpact::create<Pimpact::ScalarField>( cSpace );

  auto op = Pimpact::create< Pimpact::TransferOp<FSpaceT,CSpaceT> >( fSpace, cSpace );

  // test
  fx->initField( Pimpact::Poiseuille2D_inX );
  cx->random();

  op->apply( *fx, *cx );

  TEST_FLOATING_EQUALITY( fx->norm(Belos::OneNorm), cx->norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::TwoNorm), cx->norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::InfNorm), cx->norm(Belos::InfNorm), eps );

  fx->random();

  op->apply( *cx, *fx );

  TEST_FLOATING_EQUALITY( fx->norm(Belos::OneNorm), cx->norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::TwoNorm), cx->norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::InfNorm), cx->norm(Belos::InfNorm), eps );

}



TEUCHOS_UNIT_TEST( BasicOperator, VectorFieldOpWrap ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

  using FSpaceT = Pimpact::Space<ST,OT,d,4>;
  using CSpaceT = Pimpact::Space<ST,OT,d,2>;

  auto fSpace = Pimpact::createSpace<ST,OT,d,4>( pl );
  auto cSpace = Pimpact::createSpace<ST,OT,d,2>( pl );

  auto fx = Pimpact::create<Pimpact::VectorField>( fSpace );
  auto cx = Pimpact::create<Pimpact::VectorField>( cSpace );

  auto op = Pimpact::create< Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<FSpaceT,CSpaceT> > >( fSpace, cSpace );

  // test
  fx->initField();
  fx->getFieldPtr( Pimpact::U )->initField( Pimpact::Poiseuille2D_inX );
  cx->random();

  op->apply( *fx, *cx );

  TEST_FLOATING_EQUALITY( fx->norm(Belos::OneNorm), cx->norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::TwoNorm), cx->norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::InfNorm), cx->norm(Belos::InfNorm), eps );

  fx->random();

  op->apply( *cx, *fx );

  TEST_FLOATING_EQUALITY( fx->norm(Belos::OneNorm), cx->norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::TwoNorm), cx->norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::InfNorm), cx->norm(Belos::InfNorm), eps );

}



TEUCHOS_UNIT_TEST( BasicOperator, GradOp ) {

	if( domain!=1 ) {
		Pimpact::setBoundaryConditions( pl, domain );
		pl->set( "dim", dim );

		pl->set( "lx", lx );
		pl->set( "ly", ly );
		pl->set( "lz", lz );

		//  grid size
		pl->set("nx", nx );
		pl->set("ny", ny );
		pl->set("nz", nz );
		pl->set("nf", nf );

		// grid stretching
		if( sx!=0 ) {
			pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
			pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

		auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );


		auto p = Pimpact::create<Pimpact::ScalarField>( space );
		auto v = Pimpact::create<Pimpact::VectorField>( space );
		auto sol = v->clone();

		auto op = Pimpact::create<Pimpact::GradOp>( space );

		if( print )
			op->print();

		// op in x test
		p->initField( Pimpact::Grad2D_inX );
		sol->initField();
		sol->getFieldPtr( Pimpact::X )->initField( Pimpact::ConstField, 1. );
		v->random();

		op->apply(*p,*v);

		sol->add( 1., *sol, -1., *v );

		ST error = sol->norm( Belos::InfNorm );
		if( 0==rank )
			std::cout << "error: " << error << "\n";
		TEST_EQUALITY( error<eps, true );

		// op in y test
		p->initField( Pimpact::Grad2D_inY );
		sol->initField();
		sol->getFieldPtr( Pimpact::V )->initField( Pimpact::ConstField, 1. );
		v->random();

		op->apply(*p,*v);

		sol->add( 1., *sol, -1., *v );

		error = sol->norm( Belos::InfNorm );
		if( 0==rank )
			std::cout << "error: " << error << "\n";
		TEST_EQUALITY( error<eps, true );

		// op in z test
		p->initField( Pimpact::Grad2D_inZ );
		sol->initField();
		sol->getFieldPtr( Pimpact::W )->initField( Pimpact::ConstField, 1. );
		v->random();

		op->apply(*p,*v);

		sol->add( 1., *sol, -1., *v );

		error = sol->norm( Belos::InfNorm );
		if( 0==rank )
			std::cout << "error: " << error << "\n";
		TEST_EQUALITY( error<eps, true );
	}

}



TEUCHOS_UNIT_TEST( BasicOperator, HelmholtzOp ) {

	if( domain!= 1 ) {
		Pimpact::setBoundaryConditions( pl, domain );
		pl->set( "dim", dim );

		pl->set( "lx", lx );
		pl->set( "ly", ly );
		pl->set( "lz", lz );

		//  grid size
		pl->set("nx", nx );
		pl->set("ny", ny );
		pl->set("nz", nz );
		pl->set("nf", nf );

		// grid stretching
		if( sx!=0 ) {
			pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
			pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

		ST mulI = 5.;
		ST mulL = 3.;

		pl->set<ST>( "alpha2", mulI );
		pl->set<ST>( "Re", 1./mulL );

		auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

		auto x = Pimpact::create<Pimpact::VectorField>(space);
		auto bs= Pimpact::create<Pimpact::VectorField>(space);
		auto b = Pimpact::create<Pimpact::VectorField>(space);

		auto op= Pimpact::create<Pimpact::HelmholtzOp>( space );

		if( print )
			op->print();

		// test in x direction
		x->initField();
		x->getFieldPtr( Pimpact::U )->initField( Pimpact::Poiseuille2D_inX );
		bs->init( Teuchos::tuple( 8./std::pow(space->getDomainSize()->getSize(Pimpact::Y),2), 0., 0. ) );
		bs->add( mulI, *x, mulL, *bs );

		auto para = Teuchos::parameterList();
		para->set<ST>( "mulI", mulI );
		op->setParameter( para );
		op->apply( *x, *b );

		bs->add( 1., *bs, -1, *b);

		ST error = bs->norm( Belos::InfNorm );
		if( 0==rank )
			std::cout << "error(poiseuilleX): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			bs->print(
					*Pimpact::createOstream( "error_helm_poiX_r"+r+".txt" ));
		}

		// test in y direction
		x->initField();
		x->getFieldPtr( Pimpact::V )->initField( Pimpact::Poiseuille2D_inY );
		bs->init( Teuchos::tuple( 0., 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2), 0. ) );
		bs->add( mulI, *x, mulL, *bs );

		op->apply( *x, *b );

		bs->add( 1., *bs, -1, *b );

		error = bs->norm( Belos::InfNorm );
		if( 0==rank )
			std::cout << "error(poiseuilleY): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			bs->print(
					*Pimpact::createOstream( "error_helm_poiY_r"+r+".txt" ));
		}

		// test in z direction
		x->initField(  );
		x->getFieldPtr( Pimpact::W )->initField( Pimpact::Poiseuille2D_inZ );
		bs->init( Teuchos::tuple( 0., 0., 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2) ) );
		bs->add( mulI, *x, mulL, *bs );

		op->apply( *x, *b );

		bs->add( 1., *bs, -1, *b);

		error = bs->norm( Belos::InfNorm );
		if( 0==rank )
			std::cout << "error(poiseuilleZ): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			bs->print(
					*Pimpact::createOstream( "error_helm_poiZ_r"+r+".txt" ));
		}

		// the circle XY test
		x->initField();
		x->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY, -1. );
		x->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX,  1. );
		bs->initField();
		bs->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY, -1. );
		bs->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX,  1. );
		bs->scale( mulI );

		op->apply( *x, *b );

		bs->add( 1., *bs, -1, *b );

		error = bs->norm( Belos::InfNorm );
		if( 0==rank )
			std::cout << "error(circleXY): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			bs->print(
					*Pimpact::createOstream( "error_helm_circXY_r"+r+".txt" ));
		}

		// the circle XZ test
		x->initField();
		x->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY, -1. );
		x->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX,  1. );
		bs->initField();
		bs->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY, -1. );
		bs->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX,  1. );

		bs->scale( mulI );

		op->apply( *x, *b );

		bs->add( 1., *bs, -1, *b );

		error = bs->norm( Belos::InfNorm );
		if( 0==rank )
			std::cout << "error(circleXZ): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			bs->print(
					*Pimpact::createOstream( "error_helm_circXZ_r"+r+".txt" ));
		}
	}
}


TEUCHOS_UNIT_TEST( BasicOperator, DivGradOp2M ) {

	OT nx = 9;
	OT ny = 7;
	OT nz = 7;

	//OT nx = 9;
	//OT ny = 9;
	//OT nz = 9;

	Pimpact::setBoundaryConditions( pl, domain );
	pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.2 );
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

	auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	auto x = Pimpact::create<Pimpact::ScalarField>( space );
	auto x2= Pimpact::create<Pimpact::ScalarField>( space );
	auto b = Pimpact::create<Pimpact::ScalarField>( space );
	auto ones = Pimpact::create<Pimpact::ScalarField>( space );
	auto diag = Pimpact::create<Pimpact::ScalarField>( space );

	auto op = Pimpact::create<Pimpact::DivGradOp>( space );

	ones->init( 1. );
	op->applyInvDiag( *ones, *diag );

	OT nxS = space->end(Pimpact::S,Pimpact::X) - space->begin(Pimpact::S,Pimpact::X) + 1;
	OT nyS = space->end(Pimpact::S,Pimpact::Y) - space->begin(Pimpact::S,Pimpact::Y) + 1;
	OT nzS = space->end(Pimpact::S,Pimpact::Z) - space->begin(Pimpact::S,Pimpact::Z) + 1;
	OT nS = nxS*nyS*nzS;


	Teuchos::SerialDenseMatrix<OT,ST> DJG( nS, nS );

	// test DJG
	{
		Teuchos::RCP<std::ostream> output = Pimpact::createOstream( "DivGradOp.txt" );
		*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

		for( OT k=space->begin(Pimpact::S,Pimpact::Z); k<=space->end(Pimpact::S,Pimpact::Z); ++k )
			for( OT j=space->begin(Pimpact::S,Pimpact::Y); j<=space->end(Pimpact::S,Pimpact::Y); ++j )
				for( OT i=space->begin(Pimpact::S,Pimpact::X); i<=space->end(Pimpact::S,Pimpact::X); ++i ) {
					OT II = i-space->begin(Pimpact::S,Pimpact::X)
						+ nxS*( j-space->begin(Pimpact::S,Pimpact::Y) )
						+ nxS*nyS*( k-space->begin(Pimpact::S,Pimpact::Z) );
					x->init( 0. );
					x->at(i,j,k) = 1.;
					op->apply( *x, *b );

					for( OT kk=space->begin(Pimpact::S,Pimpact::Z); kk<=space->end(Pimpact::S,Pimpact::Z); ++kk )
						for( OT jj=space->begin(Pimpact::S,Pimpact::Y); jj<=space->end(Pimpact::S,Pimpact::Y); ++jj )
							for( OT ii=space->begin(Pimpact::S,Pimpact::X); ii<=space->end(Pimpact::S,Pimpact::X); ++ii ) {
								OT JJ = ii-space->begin(Pimpact::S,Pimpact::X)
									+ nxS*( jj-space->begin(Pimpact::S,Pimpact::Y) )
									+ nxS*nyS*( kk-space->begin(Pimpact::S,Pimpact::Z) );
								x2->init( 0. );
								x2->at(ii,jj,kk) = 1.;

								DJG(II,JJ) = x2->dot( *b );
								*output << DJG(II,JJ) << "\t";
							}
					*output << "\n";
					ST errorDiag = std::abs( 1./std::abs(DJG(II,II))-diag->at(i,j,k)) /
						std::abs(1./std::abs(DJG(II,II)) );
					if( errorDiag>=eps ) {
						std::cout << "diag("<<i<<", "<<j<< ", " <<k<<")\t" <<
							std::abs( 1./std::abs(DJG(II,II))-diag->at(i,j,k)
									)/std::abs(1./std::abs(DJG(II,II))) << "\n";
					}
					TEST_EQUALITY( errorDiag<eps, true );
				}
	}

	// test DJG^T
	Teuchos::SerialDenseMatrix<OT,ST> DJGT( nx*ny*nz, nx*ny*nz );
	{
		Teuchos::RCP<std::ostream> output = Pimpact::createOstream( "DivGradOpT.txt" );
		*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

		for( OT k=space->begin(Pimpact::S,Pimpact::Z); k<=space->end(Pimpact::S,Pimpact::Z); ++k )
			for( OT j=space->begin(Pimpact::S,Pimpact::Y); j<=space->end(Pimpact::S,Pimpact::Y); ++j )
				for( OT i=space->begin(Pimpact::S,Pimpact::X); i<=space->end(Pimpact::S,Pimpact::X); ++i ) {
					OT II = i-space->begin(Pimpact::S,Pimpact::X)
						+ nx*( j-space->begin(Pimpact::S,Pimpact::Y) )
						+ nx*ny*( k-space->begin(Pimpact::S,Pimpact::Z) );
					x->init( 0. );
					x->at(i,j,k) = 1.;
					op->apply( *x, *b, Belos::TRANS );

					for( OT kk=space->begin(Pimpact::S,Pimpact::Z); kk<=space->end(Pimpact::S,Pimpact::Z); ++kk )
						for( OT jj=space->begin(Pimpact::S,Pimpact::Y); jj<=space->end(Pimpact::S,Pimpact::Y); ++jj )
							for( OT ii=space->begin(Pimpact::S,Pimpact::X); ii<=space->end(Pimpact::S,Pimpact::X); ++ii ) {
								OT JJ = ii-space->begin(Pimpact::S,Pimpact::X)
									+ nx*( jj-space->begin(Pimpact::S,Pimpact::Y) )
									+ nx*ny*( kk-space->begin(Pimpact::S,Pimpact::Z) );
								x2->init( 0. );
								x2->at(ii,jj,kk) = 1.;
								DJGT(JJ,II) = x2->dot( *b );
								*output << DJGT(JJ,II) << "\t";
							}
					*output << "\n";
				}
	}
	DJG -= DJGT;

	ST errOne = DJG.normOne();
	std::cout << "\n\n||DJG-DJG^TT||_1 = " << errOne << "\n";
	TEST_EQUALITY( errOne<eps, true );

	ST errInf = DJG.normInf();
	std::cout << "||DJG-DJG^TT||_infty = " << errInf << "\n";
	TEST_EQUALITY( errInf<eps, true );
}



TEUCHOS_UNIT_TEST( BasicOperator, DivOp2M ) {

	OT nx = 7;
	OT ny = 7;
	OT nz = 7;

	//OT nx = 9;
	//OT ny = 9;
	//OT nz = 9;

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.2 );
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

  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	auto x = Pimpact::create<Pimpact::VectorField>( space );
	auto x2= Pimpact::create<Pimpact::ScalarField>( space );
	auto b = Pimpact::create<Pimpact::ScalarField>( space );
	auto ones = Pimpact::create<Pimpact::ScalarField>( space );

  auto op = Pimpact::create<Pimpact::DivOp>( space );


	OT nxS = space->end(Pimpact::S,Pimpact::X,true) - space->begin(Pimpact::S,Pimpact::X,true) + 1;
	OT nyS = space->end(Pimpact::S,Pimpact::Y,true) - space->begin(Pimpact::S,Pimpact::Y,true) + 1;
	OT nzS = space->end(Pimpact::S,Pimpact::Z,true) - space->begin(Pimpact::S,Pimpact::Z,true) + 1;
	OT nS = nxS*nyS*nzS;

	OT nxU = space->end(Pimpact::U,Pimpact::X,true) - space->begin(Pimpact::U,Pimpact::X,true) + 1;
	OT nyU = space->end(Pimpact::U,Pimpact::Y,true) - space->begin(Pimpact::U,Pimpact::Y,true) + 1;
	OT nzU = space->end(Pimpact::U,Pimpact::Z,true) - space->begin(Pimpact::U,Pimpact::Z,true) + 1;
	OT nU = nxU*nyU*nzU;

	OT nxV = space->end(Pimpact::V,Pimpact::X,true) - space->begin(Pimpact::V,Pimpact::X,true) + 1;
	OT nyV = space->end(Pimpact::V,Pimpact::Y,true) - space->begin(Pimpact::V,Pimpact::Y,true) + 1;
	OT nzV = space->end(Pimpact::V,Pimpact::Z,true) - space->begin(Pimpact::V,Pimpact::Z,true) + 1;
	OT nV = nxV*nyV*nzV;

	OT nxW = space->end(Pimpact::W,Pimpact::X,true) - space->begin(Pimpact::W,Pimpact::X,true) + 1;
	OT nyW = space->end(Pimpact::W,Pimpact::Y,true) - space->begin(Pimpact::W,Pimpact::Y,true) + 1;
	OT nzW = space->end(Pimpact::W,Pimpact::Z,true) - space->begin(Pimpact::W,Pimpact::Z,true) + 1;
	OT nW = nxW*nyW*nzW;

	Teuchos::Tuple<OT,3> NX = Teuchos::tuple( nxU,nxV,nxW );
	Teuchos::Tuple<OT,3> NY = Teuchos::tuple( nyU,nyV,nyW );
	Teuchos::Tuple<OT,3> NZ = Teuchos::tuple( nzU,nzV,nzW );
	Teuchos::Tuple<OT,3> offset = Teuchos::tuple( 0, nU, nU+nV );

	Teuchos::SerialDenseMatrix<OT,ST> D  ( nS      , nU+nV+nW );

	// test D

	Teuchos::RCP<std::ostream> output = Pimpact::createOstream( "DivOp.txt" );
	*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

	for( int dir=0; dir<3; ++dir ) {
		for( OT k=space->begin(dir,Pimpact::Z,true); k<=space->end(dir,Pimpact::Z,true); ++k )
			for( OT j=space->begin(dir,Pimpact::Y,true); j<=space->end(dir,Pimpact::Y,true); ++j )
				for( OT i=space->begin(dir,Pimpact::X,true); i<=space->end(dir,Pimpact::X,true); ++i ) {
					OT II =               i-space->begin(dir,Pimpact::X,true)
						+ NX[dir]*(         j-space->begin(dir,Pimpact::Y,true) )
						+ NX[dir]*NY[dir]*( k-space->begin(dir,Pimpact::Z,true) ) + offset[dir];
					x->initField();
					x->getFieldPtr(dir)->at(i,j,k) = 1.;
					op->apply( *x, *b );

					for( OT kk=space->begin(Pimpact::S,Pimpact::Z,true); kk<=space->end(Pimpact::S,Pimpact::Z,true); ++kk )
						for( OT jj=space->begin(Pimpact::S,Pimpact::Y,true); jj<=space->end(Pimpact::S,Pimpact::Y,true); ++jj )
							for( OT ii=space->begin(Pimpact::S,Pimpact::X,true); ii<=space->end(Pimpact::S,Pimpact::X,true); ++ii ) {
								OT JJ =       ii-space->begin(Pimpact::S,Pimpact::X,true)
									+ nxS*(     jj-space->begin(Pimpact::S,Pimpact::Y,true) )
									+ nxS*nyS*( kk-space->begin(Pimpact::S,Pimpact::Z,true) );
								x2->initField();
								x2->at(ii,jj,kk) = 1.;

								D(JJ,II) = x2->dot( *b );
								*output << D(JJ,II) << "\t";
							}
					*output << "\n";
				}
	}
}



TEUCHOS_UNIT_TEST( BasicOperator, GradOp2M ) {

	OT nx = 7;
	OT ny = 7;
	OT nz = 7;

	//OT nx = 9;
	//OT ny = 9;
	//OT nz = 9;

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.2 );
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

  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	auto x = Pimpact::create<Pimpact::VectorField>( space );
	auto x2= Pimpact::create<Pimpact::ScalarField>( space );
	auto b = Pimpact::create<Pimpact::VectorField>( space );

  auto op = Pimpact::create<Pimpact::GradOp>( space );


	OT nxS = space->end(Pimpact::S,Pimpact::X,true) - space->begin(Pimpact::S,Pimpact::X,true) + 1;
	OT nyS = space->end(Pimpact::S,Pimpact::Y,true) - space->begin(Pimpact::S,Pimpact::Y,true) + 1;
	OT nzS = space->end(Pimpact::S,Pimpact::Z,true) - space->begin(Pimpact::S,Pimpact::Z,true) + 1;
	OT nS = nxS*nyS*nzS;

	OT nxU = space->end(Pimpact::U,Pimpact::X,true) - space->begin(Pimpact::U,Pimpact::X,true) + 1;
	OT nyU = space->end(Pimpact::U,Pimpact::Y,true) - space->begin(Pimpact::U,Pimpact::Y,true) + 1;
	OT nzU = space->end(Pimpact::U,Pimpact::Z,true) - space->begin(Pimpact::U,Pimpact::Z,true) + 1;
	OT nU = nxU*nyU*nzU;

	OT nxV = space->end(Pimpact::V,Pimpact::X,true) - space->begin(Pimpact::V,Pimpact::X,true) + 1;
	OT nyV = space->end(Pimpact::V,Pimpact::Y,true) - space->begin(Pimpact::V,Pimpact::Y,true) + 1;
	OT nzV = space->end(Pimpact::V,Pimpact::Z,true) - space->begin(Pimpact::V,Pimpact::Z,true) + 1;
	OT nV = nxV*nyV*nzV;

	OT nxW = space->end(Pimpact::W,Pimpact::X,true) - space->begin(Pimpact::W,Pimpact::X,true) + 1;
	OT nyW = space->end(Pimpact::W,Pimpact::Y,true) - space->begin(Pimpact::W,Pimpact::Y,true) + 1;
	OT nzW = space->end(Pimpact::W,Pimpact::Z,true) - space->begin(Pimpact::W,Pimpact::Z,true) + 1;
	OT nW = nxW*nyW*nzW;

	Teuchos::Tuple<OT,3> NX = Teuchos::tuple( nxU,nxV,nxW );
	Teuchos::Tuple<OT,3> NY = Teuchos::tuple( nyU,nyV,nyW );
	Teuchos::Tuple<OT,3> NZ = Teuchos::tuple( nzU,nzV,nzW );
	Teuchos::Tuple<OT,3> offset = Teuchos::tuple( 0, nU, nU+nV );

	Teuchos::SerialDenseMatrix<OT,ST> G  ( nU+nV+nW, nS);

	// test G

	Teuchos::RCP<std::ostream> output = Pimpact::createOstream( "GradOp.txt" );
	*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

	for( OT kk=space->begin(Pimpact::S,Pimpact::Z,true); kk<=space->end(Pimpact::S,Pimpact::Z,true); ++kk )
		for( OT jj=space->begin(Pimpact::S,Pimpact::Y,true); jj<=space->end(Pimpact::S,Pimpact::Y,true); ++jj )
			for( OT ii=space->begin(Pimpact::S,Pimpact::X,true); ii<=space->end(Pimpact::S,Pimpact::X,true); ++ii ) {
				OT JJ =       ii-space->begin(Pimpact::S,Pimpact::X,true)
					+ nxS*(     jj-space->begin(Pimpact::S,Pimpact::Y,true) )
					+ nxS*nyS*( kk-space->begin(Pimpact::S,Pimpact::Z,true) );
				x2->initField();
				x2->at(ii,jj,kk) = 1.;
				op->applyG( *x2, *b );

				for( int dir=0; dir<3; ++dir ) {
					for( OT k=space->begin(dir,Pimpact::Z,true); k<=space->end(dir,Pimpact::Z,true); ++k )
						for( OT j=space->begin(dir,Pimpact::Y,true); j<=space->end(dir,Pimpact::Y,true); ++j )
							for( OT i=space->begin(dir,Pimpact::X,true); i<=space->end(dir,Pimpact::X,true); ++i ) {
								OT II =               i-space->begin(dir,Pimpact::X,true)
									+ NX[dir]*(         j-space->begin(dir,Pimpact::Y,true) )
									+ NX[dir]*NY[dir]*( k-space->begin(dir,Pimpact::Z,true) ) + offset[dir];
								x->initField();
								x->getFieldPtr(dir)->at(i,j,k) = 1.;


								G(II,JJ) = x->dot( *b, true );
								//*output << G(II,JJ) << "\t";
							}
					//*output << "\n";
				}
	}
	*output << G;

	// test J
	Teuchos::RCP<std::ostream> outputJ = Pimpact::createOstream( "JOp.txt" );
	*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

	Teuchos::SerialDenseMatrix<OT,ST> J  ( nU+nV+nW, nU+nV+nW );
	for( int dir=0; dir<3; ++dir ) {
		for( OT kk=space->begin(dir,Pimpact::Z,true); kk<=space->end(dir,Pimpact::Z,true); ++kk )
			for( OT jj=space->begin(dir,Pimpact::Y,true); jj<=space->end(dir,Pimpact::Y,true); ++jj )
				for( OT ii=space->begin(dir,Pimpact::X,true); ii<=space->end(dir,Pimpact::X,true); ++ii ) {
					OT JJ =           ii-space->begin(dir,Pimpact::X,true)
						+ NX[dir]*(     jj-space->begin(dir,Pimpact::Y,true) )
						+ NX[dir]*NY[dir]*( kk-space->begin(dir,Pimpact::Z,true) )+offset[dir];
					b->initField();
					b->getFieldPtr(dir)->at(ii,jj,kk) = 1.;
					op->applyJ( *b );

					for( OT k=space->begin(dir,Pimpact::Z,true); k<=space->end(dir,Pimpact::Z,true); ++k )
						for( OT j=space->begin(dir,Pimpact::Y,true); j<=space->end(dir,Pimpact::Y,true); ++j )
							for( OT i=space->begin(dir,Pimpact::X,true); i<=space->end(dir,Pimpact::X,true); ++i ) {
								OT II =               i-space->begin(dir,Pimpact::X,true)
									+ NX[dir]*(         j-space->begin(dir,Pimpact::Y,true) )
									+ NX[dir]*NY[dir]*( k-space->begin(dir,Pimpact::Z,true) ) + offset[dir];
								x->initField();
								x->getFieldPtr(dir)->at(i,j,k) = 1.;

								J(II,JJ) = x->dot( *b, true );
								//*outputJ << J(II,JJ) << "\t";
							}
					//*outputJ << "\n";
				}
	}
	*outputJ << J ;
}



TEUCHOS_UNIT_TEST( BasicOperator, DivGradO2Op ) {

	Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

	const int dNC = 2;

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.2 );
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


  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  auto x = Pimpact::create<Pimpact::ScalarField>( space );
  auto x2= Pimpact::create<Pimpact::ScalarField>( space );
  auto b = Pimpact::create<Pimpact::ScalarField>( space );


  auto op = Pimpact::create<Pimpact::DivGradOp>( space );
  auto op2 = Pimpact::create<Pimpact::DivGradO2Op>( space );

	if( print ) {
		op->print();
		op2->print();
	}

  // zero test
  b->initField( Pimpact::ConstField, 0. );
  x->random();

  op->apply( *b, *x );

	ST error = x->norm( Belos::InfNorm );
	if( 0==rank )
		std::cout << "error(0): " << error << "\n";
  TEST_EQUALITY( error<eps, true );

  b->initField( Pimpact::ConstField, 0. );
  x->random();

  op2->apply( *b, *x );

	error = x->norm( Belos::InfNorm );
	if( 0==rank )
		std::cout << "error(0) O2: " << error << "\n";
  TEST_EQUALITY( error<eps, true );

	// one test
  b->initField( Pimpact::ConstField, 1. );
  x->random();

  op->apply( *b, *x );

	error = x->norm( Belos::InfNorm );
	if( 0==rank )
		std::cout << "error(1): " << error << "\n";
	TEST_EQUALITY( error<eps, true );

  b->initField( Pimpact::ConstField, 1. );
  x->random();

  op2->apply( *b, *x );

	error = x->norm( Belos::InfNorm );
	if( 0==rank )
		std::cout << "error(1) O2: " << error << "\n";
	if( error>= eps ) {
		x->print();
	}
	TEST_EQUALITY( error<eps, true );

	// consistency test
	b->random();
	b->random();
	op->apply( *b, *x );
	op2->apply( *b, *x2 );

	x2->add( 1., *x2, -1., *x );
	ST errInf = x2->norm( Belos::InfNorm );
	ST err2 = x2->norm( Belos::TwoNorm )/std::sqrt( x2->getLength() );
	if( 0==rank )
		std::cout << "consistency error random: inf: " << errInf << ", two: " << err2 << "\n";
	TEST_EQUALITY( errInf<eps, true );
	TEST_EQUALITY( err2<eps, true );
	if( errInf>=eps ) {
		x2->write();
	}

	for( int dir=0; dir<=6; ++dir ) {

		Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>( dir );

		b->initField( type );
		op->apply( *b, *x );
		op2->apply( *b, *x2 );

		x2->add( 1., *x2, -1., *x );
		errInf = x2->norm( Belos::InfNorm );
		err2 = x2->norm( Belos::TwoNorm )/std::sqrt( x2->getLength() );
		if( 0==rank )
			std::cout << "consistency error("+Pimpact::toString(type)+"): inf: " << errInf << ", two: " << err2 << "\n";
		TEST_EQUALITY( errInf<eps, true );
		TEST_EQUALITY( err2<eps, true );
		if( errInf>=eps && err2>=eps ) {
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			x2->print( *Pimpact::createOstream(
						"error_dgo2_"+Pimpact::toString(type)+"_r"+r+".txt" ) );
			x2->write(dir+1);
			x->write( (dir+1)*10 );
			x2->add( 1., *x2, 1., *x );
			x2->write( (dir+1)*100 );
		}
	}

	// InvDiag consistency test
	//x->init( 1. );
	//op->applyInvDiag( *x, *b );
	//op2->applyInvDiag( *x, *x2 );
	//b->write();
	//x2->write(1);

	//x2->add( 1., *x2, -1., *b );
	//x2->write(2);
	//ST diff = x2->norm( Belos::InfNorm );
	//if( 0==space->rankST() ) 
		//std::cout << "diff InvDiag: " << diff << "\n";
	//TEST_EQUALITY( diff<eps, true );
}



TEUCHOS_UNIT_TEST( BasicOperator, DivGradTransposeOp ) {

	//const int dNC = 2;

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.2 );
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

  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  auto xp = Pimpact::create<Pimpact::ScalarField>( space );
  auto xv = Pimpact::create<Pimpact::VectorField>( space );
  auto bp = Pimpact::create<Pimpact::ScalarField>( space );
  auto bp2 = Pimpact::create<Pimpact::ScalarField>( space );
  auto bv = Pimpact::create<Pimpact::VectorField>( space );
  auto bv2 = Pimpact::create<Pimpact::VectorField>( space );


  auto div = Pimpact::create<Pimpact::DivOp>( space );
  auto grad = Pimpact::create<Pimpact::GradOp>( space );

	auto divGrad = Pimpact::createDivGradOp( div, grad );

	if( print ) {
		div->print();
		grad->print();
	}

	xp->init(1.);
	xv->init(1.);
	////xp->random();

	div->apply(  *xp, *bv  );
	std::cout << "||div^T(ones)||" << bv->norm() << "\n";
	grad->apply( *xp, *bv2 );
	std::cout << "||grad(ones)||" << bv2->norm() << "\n";

	//bv->write();
	//bv2->write(1);

	grad->apply( *xv, *bp);
	std::cout << "||grad^T(ones)||" << bp->norm() << "\n";
	div->apply(  *xv, *bp2  );
	std::cout << "||div(ones)||" << bp2->norm() << "\n";

	//bp->write();
	//bp2->write(1);

	divGrad->apply( *xp, *bp );
	divGrad->apply( *xp, *bp2, Belos::TRANS );
	std::cout << "||divGrad(ones): " << bp->norm() << "\n";
	std::cout << "||divGrad^T(ones): " << bp2->norm() << "\n";

	//bp->write(2);
	bp2->write(3);

	ST pi2 = std::atan(1)*8;
	xp->initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
	//xp->random();
	divGrad->apply( *xp, *bp );
	divGrad->apply( *xp, *bp2, Belos::TRANS );

	//bp->write(1);
	//bp2->write(2);
	bp->add( 1., *bp, -1., *bp2 );
	std::cout << "difference(divgrad, divgrad^T): " << bp->norm() << "\n";

	bp->write(3);

}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, DivGradO2Smoother, SType ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05);
	}
	if( sy!=0 ) {
		pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	}
	if( sz!=0 ) {
		if( 1==sz ) {
			pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
			pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
		}
		else {
			pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz) );
			pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz) );
		}
	}

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	Teuchos::RCP< const Pimpact::Space<ST,OT,d,dNC> > space =
		Pimpact::createSpace<ST,OT,d,dNC>( pl );

	auto coord = space->getCoordinatesGlobal();
	if( 0==rank ) {
		Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "coord.txt" );
		coord->print( *fstream );
	}

	Teuchos::RCP< Pimpact::ScalarField< Pimpact::Space<ST,OT,d,dNC> > > x =
		Pimpact::create<Pimpact::ScalarField>( space );

  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );
	//op->print2Mat();

	Teuchos::RCP<Teuchos::ParameterList> ppl = Teuchos::parameterList();

	// compute EV
	ST evMax;
	ST evMin;

	Teuchos::RCP<Pimpact::TeuchosEigenvalues<Pimpact::DivGradO2Op<SpaceT> > > ev = 
		Teuchos::rcp( new Pimpact::TeuchosEigenvalues<Pimpact::DivGradO2Op<SpaceT> >( op ) );

	if( 0==rank )	
		std::cout << "\n";
	for( int i=0; i<3; ++i ) {
		ev->computeEV( static_cast<Pimpact::ECoord>(i), evMax, evMin );
		if( 0==rank )	
			std::cout << Pimpact::toString( static_cast<Pimpact::ECoord>(i) ) << ": " << evMax << "\t" <<evMin << "\n";
	}

	ev->computeEV( evMax, evMin );
	if( 0==rank )	
		std::cout << "glob : " << evMax << "\t" <<evMin << "\n";

	if( nx<=20 && ny<=20 && nz<=20 ) {
		ev->computeFullEV( evMax, evMin );
		if( 0==rank )	
			std::cout << "glob2: " << evMax << "\t" <<evMin << "\n";
	}

	ppl->set<int>( "numIters", sweeps );

	// good result with no stretch + own
	ppl->set<ST>( "max EV", evMax );
	ppl->set<ST>( "min EV", evMin*1.1 );
	//ppl->set<ST>( "max EV", evMax*1.  );
	//ppl->set<ST>( "min EV", evMin*2.1 );

	ppl->set<bool>( "with output", true );
	
	// JSmoother
	ppl->set<int>( "BC smoothing", 0 );
	ppl->set<OT>( "depth", 2 );

  auto smoother = Pimpact::create<SType>( op, ppl );

	// --- zero rhs test --- 
	x->random();
	b->init( 0. );

	//Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "conv_"+smoother->getLabel()+".txt", 0 );
	if( 0==rank ) {
		if( 0==rank ) {
			std::cout << "\nstep\terror\t\trate\n";
			std::cout << 0 << "\t" << 1. << "\n";
		}
	}
	ST err0 = x->norm( Belos::InfNorm );
	ST errP = err0;

	x->write( 0 );
	for( int i=1; i<=nIter; ++i ) {
		smoother->apply( *b, *x );
		x->write(i);
		ST err = x->norm( Belos::InfNorm );
		if( 0==rank ) {
				std::cout << i << "\t" << err/err0 << "\t" << err/errP << "\n";
				//fstream << i << "\t" << err/err0 << "\t" << err/errP << "\n";
		}
		errP = err;
	}

	// --- grad test ---
	auto sol = x->clone();
	auto e = x->clone();
	auto res = x->clone();
	for( int dir=1; dir<=3; ++dir ) {

		Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

		x->initField( type );
		sol->initField( type );
		sol->level();

		op->apply( *x, *b );

		//x->init( 0 );
		x->random();

		// residual
		op->apply( *x, *res );
		res->add( -1, *b, 1., *res );
		ST residual = res->norm();
		ST res0 = res->norm();
		ST resP = res0;
		//error
		res->add( 1., *sol, -1., *x );
		ST err0 = res->norm();
		ST errP = err0;
		//ST err0 = 1.;
		//ST errP =1.;

		if( 0==rank ) {
			std::cout << "\n\n\t\t\t--- " << Pimpact::toString(type) << " test ---\n";
			std::cout << "\tresidual:\trate:\t\t\terror:\t\trate:\n";
			std::cout << "\t"  << 1.<< "\t\t\t\t" << 1.  << "\n";
		}
		for( int i=0; i<nIter; ++i ) {
			smoother->apply( *b, *x );
			x->level();

			op->apply( *x, *res );
			res->add( -1, *b, 1., *res );
			residual = res->norm();
			res->add( 1., *sol, -1., *x );
			ST err = res->norm();

			if( 0==rank )
				std::cout << "\t" << residual/res0 << "\t" <<  residual/resP << "\t\t" << err/err0 << "\t" <<  err/errP  << "\n";
			resP = residual;
			errP = err;
		}
		TEST_EQUALITY( residual/res0<1.e-1, true );
	}

	// --- consistency test ---
	for( int dir=1; dir<=6; ++dir ) {

		x->initField( static_cast<Pimpact::EScalarField>(dir) );
		x->level();
		auto xp = x->clone( Pimpact::ECopyType::Deep );

		op->apply( *x, *b );

		smoother->apply( *b, *x );
		x->level();

		xp->add( 1., *xp, -1., *x );
		ST err2 = xp->norm( Belos::InfNorm )/std::sqrt( static_cast<ST>(xp->getLength()) );
		ST errInf = xp->norm( Belos::InfNorm );

		if( 0==rank ) {
			std::cout << "consistency for " << dir << ": ||" << err2 << "||_2, ||" << errInf << "||_inf\n";
		}

		TEST_EQUALITY( err2<eps, true );
		TEST_EQUALITY( errInf<eps, true );
		if( err2>eps && errInf>eps ) 
			x->write( nIter + dir + 1 );

	}

}

using JT = Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op<SpaceT> >;
//using SORT = Pimpact::DivGradO2SORSmoother< Pimpact::DivGradO2Op<SpaceT> >;
using CheT = Pimpact::Chebyshev< Pimpact::DivGradO2Op<SpaceT> >;
using LT = Pimpact::DivGradO2LSmoother< Pimpact::DivGradO2Op<SpaceT> >;

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, JT )
//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, SORT )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, CheT )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, LT )



TEUCHOS_UNIT_TEST( BasicOperator, DivGradO2Inv ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", 33 );
	pl->set("ny", 17 );
	pl->set("nz", 9 );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", 33./2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", 33./2. );
	}
	if( sy!=0 ) {
		pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Y").set<ST>( "N metr L", 17./2. );
		pl->sublist("Stretching in Y").set<ST>( "N metr U", 17./2. );
	}
	if( sz!=0 ) {
		pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in Z").set<ST>( "N metr L", 9./2. );
		pl->sublist("Stretching in Z").set<ST>( "N metr U", 9./2. );
	}

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	const int dNC=2;

  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	Teuchos::RCP< Pimpact::ScalarField< Pimpact::Space<ST,OT,d,dNC> > > x =
		Pimpact::create<Pimpact::ScalarField>( space );

	Teuchos::RCP< Pimpact::ScalarField< Pimpact::Space<ST,OT,d,dNC> > > xp = x->clone();
	Teuchos::RCP< Pimpact::ScalarField< Pimpact::Space<ST,OT,d,dNC> > > b = x->clone();

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );

	if( print )
		op->print2Mat();

  auto ppl = Teuchos::parameterList();

  auto solver = Pimpact::create<Pimpact::DivGradO2Inv>( op, ppl );

  // --- consistency test ---
	for( int field=0; field<=(domain!=1?6:0); ++field ) {
		x->initField( static_cast<Pimpact::EScalarField>(field) );
		x->level();
		xp->random();

		op->apply( *x, *b );
		solver->apply( *b, *xp );

		xp->add( 1., *xp, -1., *x );
		xp->level();

		ST err2 = xp->norm( Belos::InfNorm )/std::sqrt( static_cast<ST>(xp->getLength()) );
		ST errInf = xp->norm( Belos::InfNorm );
		if( errInf>=eps ) {
			xp->write();
			//std::cout << "rank: " << rank << "\te: " << xp->norm( Belos::InfNorm, false ) << "\n";
			//if( 0==rank )
				//xp->print();
		}


		if( 0==rank )
			std::cout << "consistency for " << Pimpact::toString(static_cast<Pimpact::EScalarField>(field)) << ": ||" << err2 << "||_2, ||" << errInf << "||_inf\n";

		if( npx*npy*npz==1 ) {
			TEST_EQUALITY( err2<eps, true );
			TEST_EQUALITY( errInf<eps, true );
		}
	}

}



TEUCHOS_UNIT_TEST( BasicOperator, ForcingOp ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

	auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto force = Pimpact::create<Pimpact::VectorField>( space );
	auto res = Pimpact::create<Pimpact::VectorField>( space );

	vel->init(1.);
	//force->initField( Pimpact::EVectorField(11) );
	force->random();
	res->random();

	auto op = Pimpact::createForcingOp( force, 1. );

	op->apply( *vel, *res );
	vel->add( 1., *res, -1., *force );

	TEST_EQUALITY( vel->norm( Belos::InfNorm )<eps, true );
}





TEUCHOS_UNIT_TEST( MultiOperator, Add2Op ) {

	Pimpact::setBoundaryConditions( pl, domain );
	pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

	auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing

	auto velx = Pimpact::create<Pimpact::VectorField>(space);


	auto x = Pimpact::createMultiField(*velx->clone(),10);
	auto y = Pimpact::createMultiField(*velx->clone(),10);

	auto op = Pimpact::createOperatorBase(
			Pimpact::createMultiOpWrap(
				Pimpact::createAdd2Op(
					Pimpact::create<Pimpact::HelmholtzOp>( space ),
					Pimpact::create<Pimpact::HelmholtzOp>( space )
					)
				)
			);

	x->random();
	op->apply( *x, *y);

}


template<class T> using WUP = Pimpact::MultiOpUnWrap<Pimpact::MultiOpWrap<T> >;


TEUCHOS_UNIT_TEST( MultiOperator, MulitOpUnWrap ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  auto x = Pimpact::create<Pimpact::VectorField>(space);
  auto y = Pimpact::create<Pimpact::VectorField>(space);

  auto op = Pimpact::createOperatorBase(
			Pimpact::createMultiOpWrap(
				Pimpact::create<Pimpact::HelmholtzOp>( space)  ) );

  auto opUW = Pimpact::create<Pimpact::MultiOpUnWrap>( op );

  auto opUW2 =
      Pimpact::create<WUP>(
          Pimpact::create<Pimpact::HelmholtzOp>(space)
      );


  //x->initField(Pimpact::RankineVortex2D );
	x->random();

  opUW2->apply( *x, *y );

}



TEUCHOS_UNIT_TEST( MultiModeOperator, HelmholtzOp ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  auto velc = Pimpact::create<Pimpact::VectorField>(space);
  auto vels = Pimpact::create<Pimpact::VectorField>(space);

  auto vel = Pimpact::create<MVF>( space );

  auto mv = Pimpact::createMultiField(*vel,10);

  auto mv2 = mv->clone(10);

  auto op = Pimpact::createMultiOpWrap(
      Pimpact::createModeOpWrap(
          Pimpact::create<Pimpact::HelmholtzOp>(space) ) );

  op->apply(*mv,*mv2);

}




TEUCHOS_UNIT_TEST( MultiModeOperator, DtModeOp ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

	pl->set<ST>("Re", 1., "Reynolds number");
	pl->set<ST>("alpha2", 1.,
			"Womersley square \alpha^2");
  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  using MF = Pimpact::MultiField< Pimpact::ModeField< Pimpact::VectorField<SpaceT> > >;
  using BOp = Pimpact::OperatorBase<MF>;


  auto velc = Pimpact::create<Pimpact::VectorField>(space);
  auto vels = Pimpact::create<Pimpact::VectorField>(space);

  auto vel = Pimpact::create<MVF>( space );

  auto mv = Pimpact::createMultiField(*vel,1);

  mv->getFieldPtr(0)->getCFieldPtr()->init(1.);
  mv->getFieldPtr(0)->getSFieldPtr()->init(0.);

  auto mv2 = mv->clone(1);
  mv2->init(0.);

  auto op = Pimpact::createMultiOperatorBase(
      Pimpact::create<Pimpact::DtModeOp>(space) );

  op->apply(*mv,*mv2);

  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm( Belos::InfNorm ), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm( Belos::InfNorm ) );
  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm( Belos::InfNorm ), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm( Belos::InfNorm ) );
  Belos::OperatorTraits<ST,MF,BOp>::Apply(*op,*mv,*mv2);

}




//TEUCHOS_UNIT_TEST( MultiModeOperator, DtLapOp ) {

	//Pimpact::setBoundaryConditions( pl, domain );
	//pl->set( "dim", dim );

	//pl->set( "lx", lx );
	//pl->set( "ly", ly );
	//pl->set( "lz", lz );

	////  grid size
	//pl->set("nx", nx );
	//pl->set("ny", ny );
	//pl->set("nz", nz );
	//pl->set("nf", nf );

	//// grid stretching
	//if( sx!=0 ) {
		//pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
	//}
	//if( sy!=0 ) {
		//pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		//pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	//}
	//if( sz!=0 ) {
		//pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		//pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	//}

	//// processor grid size
	//pl->set( "npx", npx );
	//pl->set( "npy", npy );
	//pl->set( "npz", npz );
	//pl->set( "npf", npf );

	//pl->set( "Re", 1. );
	//auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	//auto velc = Pimpact::create<Pimpact::VectorField>(space);
	//auto vels = Pimpact::create<Pimpact::VectorField>(space);

	//auto vel = Pimpact::create<MVF>( space );

	//auto mv = Pimpact::createMultiField(*vel,1);

	//mv->random();

	//auto mv2 = mv->clone(1);
	//auto mv3 = mv->clone(1);
	//mv2->init(0.);
	//mv3->init(0.);

	//auto A = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp( space, 1., 0. ) );

	//mv->random();

	//A->apply( *mv, *mv2 );

	//auto diff = mv->getConstFieldPtr(0)->getConstCFieldPtr()->clone(Pimpact::ECopyType::Shallow);

	//diff->add( 1., mv->getConstFieldPtr(0)->getConstCField(), -1.,  mv2->getConstFieldPtr(0)->getConstSField() );

	//TEST_EQUALITY( diff->norm(Belos::InfNorm)<eps, true );

	//diff->add( 1., mv->getConstFieldPtr(0)->getConstSField(), -1.,  mv2->getConstFieldPtr(0)->getConstCField() );
	//TEST_EQUALITY( diff->norm(Belos::InfNorm)<eps, true );

	//mv2->init(0.);

	//auto A2 = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp( space, 0., 1. ) );
	//auto A3 = Pimpact::createMultiOpWrap( Pimpact::createModeOpWrap( Pimpact::create<Pimpact::HelmholtzOp>( space ) ) );

	//A2->apply(*mv,*mv2);
	//A3->apply(*mv,*mv3);


	//auto diff2 = mv2->clone();
	//diff2->add( 1., *mv2, -1., *mv3 );

	//TEST_EQUALITY( diff2->norm(Belos::InfNorm)<eps, true );

//}




//TEUCHOS_UNIT_TEST( MultiModeOperator, TripleCompostion ) {

  //Pimpact::setBoundaryConditions( pl, domain );
  //pl->set( "dim", dim );

	//pl->set( "lx", lx );
	//pl->set( "ly", ly );
	//pl->set( "lz", lz );

	////  grid size
	//pl->set("nx", nx );
	//pl->set("ny", ny );
	//pl->set("nz", nz );
	//pl->set("nf", nf );

	//// grid stretching
	//if( sx!=0 ) {
		//pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
	//}
	//if( sy!=0 ) {
		//pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		//pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	//}
	//if( sz!=0 ) {
		//pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		//pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	//}

  //// processor grid size
  //pl->set( "npx", npx );
  //pl->set( "npy", npy );
  //pl->set( "npz", npz );
  //pl->set( "npf", npf );

	//pl->set( "Re", 1. );
	//pl->set( "alpha2", 1. );
  //auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  //using MVF = Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > >;
  //using MSF = Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<SpaceT> > >;

  //auto X = Pimpact::create<MSF>( space );
  //auto B = Pimpact::create<MSF>( space );

  //auto temp = Pimpact::create<MVF>( space );

  //X->init(0.);
////  B->random();
  //B->init(1.);


  //auto H =
      //Pimpact::createMultiOperatorBase(
          //Pimpact::createDtLapOp( space, 0., 10. ) );

  //H->apply( *temp, *temp );

  //// Make an empty new parameter list.
  ////auto solverParams = Teuchos::parameterList();
  //auto solverParams = Pimpact::createLinSolverParameter( "GMRES", 9.e-1 );

  //// Create the Pimpact::LinearSolver solver.
  //auto Hprob = Pimpact::createLinearProblem<MVF>( H, temp, temp, solverParams,"GMRES" );
  //auto Hinv  = Pimpact::createInverseOperator( Hprob );

  //auto schur = Pimpact::createTripleCompositionOp(
      //Pimpact::createMultiModeOpWrap( Pimpact::create<Pimpact::DivOp>( space ) ),
      //Hinv,
      //Pimpact::createMultiModeOpWrap( Pimpact::create<Pimpact::GradOp>( space ) )
  //);

  //schur->apply( *B, *X );

//}




//TEUCHOS_UNIT_TEST( MultiModeOperator, InverseOperator ) {

	//Pimpact::setBoundaryConditions( pl, domain );
	//pl->set( "dim", dim );

	//pl->set( "lx", lx );
	//pl->set( "ly", ly );
	//pl->set( "lz", lz );

	////  grid size
	//pl->set("nx", nx );
	//pl->set("ny", ny );
	//pl->set("nz", nz );
	//pl->set("nf", nf );

	//// grid stretching
	//if( sx!=0 ) {
		//pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
	//}
	//if( sy!=0 ) {
		//pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		//pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	//}
	//if( sz!=0 ) {
		//pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		//pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	//}

	//// processor grid size
	//pl->set( "npx", npx );
	//pl->set( "npy", npy );
	//pl->set( "npz", npz );
	//pl->set( "npf", npf );

	//auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	//using Teuchos::ParameterList;
	//using Teuchos::parameterList;
	//using Teuchos::RCP;
	//using Teuchos::rcp; // Save some typing


	//using MVF = Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > >;

	//auto X = Pimpact::create<MVF>( space );
	//auto B = Pimpact::create<MVF>( space );

	//X->init(0.);
	//B->random();

	//auto op = Pimpact::createMultiModeOperatorBase<MVF>(
			//Pimpact::create<Pimpact::HelmholtzOp>(space) );

	//// Make an empty new parameter list.
	//auto solverParams = Pimpact::createLinSolverParameter( "CG", 9.e-1 );

	//// Create the Pimpact::LinearSolver solver.
	//auto prob =
		//Pimpact::createLinearProblem<MVF>(
				//op,
				//X,X,
				////          Pimpact::createMultiField(X->getFieldPtr(0)->getCFieldPtr()),
				////          Pimpact::createMultiField(B->getFieldPtr(0)->getCFieldPtr()),
				//solverParams );

	//auto opinv = Pimpact::createInverseOperator<MVF>( prob );
	//auto opp = Pimpact::createOperatorBase( opinv );

//}





//TEUCHOS_UNIT_TEST( MultiModeOperator, EddyPrec ) {

	//Pimpact::setBoundaryConditions( pl, domain );
	//pl->set( "dim", dim );

	//pl->set( "lx", lx );
	//pl->set( "ly", ly );
	//pl->set( "lz", lz );

	////  grid size
	//pl->set("nx", nx );
	//pl->set("ny", ny );
	//pl->set("nz", nz );
	//pl->set("nf", nf );

	//// grid stretching
	//if( sx!=0 ) {
		//pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
	//}
	//if( sy!=0 ) {
		//pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		//pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	//}
	//if( sz!=0 ) {
		//pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		//pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	//}

  //// processor grid size
  //pl->set( "npx", npx );
  //pl->set( "npy", npy );
  //pl->set( "npz", npz );
  //pl->set( "npf", npf );

	////pl->set("alpha2",1.);
	////pl->set("Re",1./14.);
	//auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	//auto temp = Pimpact::createMultiModeVectorField( space );

	//auto X = Pimpact::createMultiModeVectorField( space );
	//auto B = Pimpact::createMultiModeVectorField( space );

	//X->init(0.);
	//B->init(1.);

	//// Make an empty new parameter list.
	//auto solverParams = Pimpact::createLinSolverParameter("CG",9.e-1);

	//// Create the Pimpact::LinearSolver solver.
	//auto A =
		//Pimpact::createMultiOperatorBase(
				//Pimpact::create<Pimpact::HelmholtzOp>( space )
				//);

	//A->apply( *Pimpact::createMultiField( B->getFieldPtr(0)->getCFieldPtr() ), 
			//*Pimpact::createMultiField( X->getFieldPtr(0)->getCFieldPtr() ) );

	//auto prob =
		//Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::VectorField<SpaceT> > >(
				//A,
				//Teuchos::null,
				//Teuchos::null,
				//solverParams,
				//"CG" );

	//prob->solve( Pimpact::createMultiField( B->getFieldPtr(0)->getCFieldPtr() ), 
			//Pimpact::createMultiField( X->getFieldPtr(0)->getCFieldPtr() ) );
	
	//auto invOp = Pimpact::createInverseOperator( prob );

	//invOp->apply( Pimpact::createMultiField( B->getFieldPtr(0)->getCFieldPtr() ), 
			//Pimpact::createMultiField( X->getFieldPtr(0)->getCFieldPtr() ) );

	//auto op =
		//Pimpact::create<Pimpact::EddyPrec>(
				//Pimpact::createMultiOpUnWrap(
					//Pimpact::createInverseOperatorBase(prob) ) );
////	auto op =
////		Pimpact::create<Pimpact::EddyPrec>( 
////				Pimpact::createMultiOpUnWrap( A )	);
	
	//op->apply( X->getField(0), B->getField(0) );
	
	//auto schur = Pimpact::createMultiOperatorBase( op );

	//schur->apply( *B, *X );

//}




//TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicConvectionOp ) {

	//Pimpact::setBoundaryConditions( pl, domain );
	//pl->set( "dim", dim );

	//pl->set( "lx", lx );
	//pl->set( "ly", ly );
	//pl->set( "lz", lz );

	//pl->set<bool>( "spectral in time", true );
  //Pimpact::setBoundaryConditions( pl, domain );
  //pl->set( "dim", dim );

	////  grid size
	//pl->set("nx", nx );
	//pl->set("ny", ny );
	//pl->set("nz", nz );
	//pl->set("nf", nf );

	//// grid stretching
	//if( sx!=0 ) {
		//pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
	//}
	//if( sy!=0 ) {
		//pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		//pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	//}
	//if( sz!=0 ) {
		//pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		//pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	//}

  //// processor grid size
  //pl->set( "npx", npx );
  //pl->set( "npy", npy );
  //pl->set( "npz", npz );
  //pl->set( "npf", npf );

  //auto space = Pimpact::createSpace<ST,OT,4,dNC>( pl );

  //auto vel = Pimpact::create<Pimpact::VectorField>( space );

  //auto mv1 = Pimpact::createMultiHarmonicVectorField( space );
  //auto mv2 = Pimpact::createMultiHarmonicVectorField( space );

  //auto op = Pimpact::createMultiHarmonicConvectionOp( space );

  //op->assignField( *mv1 );
  //op->apply( *mv1, *mv2 );

//}



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicOpWrap ) {

	using SpaceT = Pimpact::Space<ST,OT,4,dNC>;

	Pimpact::setBoundaryConditions( pl, domain );
	pl->set( "dim", dim );

	pl->set<bool>( "spectral in time", true );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	if( sx!=0 ) {
		pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
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

  auto space = Pimpact::createSpace<ST,OT,4,dNC>( pl );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto mv1 = Pimpact::createMultiHarmonicVectorField( space );
  auto mv2 = Pimpact::createMultiHarmonicVectorField( space );

  auto op = Pimpact::createMultiHarmonicOpWrap< Pimpact::HelmholtzOp<SpaceT> >(
      Pimpact::create<Pimpact::HelmholtzOp>(space) );

  op->apply( *mv1, *mv2 );

}



//TEUCHOS_UNIT_TEST( CompoundOperator, CompoundOpWrap ) {

	//using SpaceT = Pimpact::Space<ST,OT,4,dNC>;

	//using VF = Pimpact::VectorField<SpaceT>;
	//using SF = Pimpact::ScalarField<SpaceT>;

	//Pimpact::setBoundaryConditions( pl, domain );
	//pl->set( "dim", dim );

	//pl->set( "lx", lx );
	//pl->set( "ly", ly );
	//pl->set( "lz", lz );

	//pl->set<bool>( "spectral in time", true );

	////  grid size
	//pl->set("nx", nx );
	//pl->set("ny", ny );
	//pl->set("nz", nz );
	//pl->set("nf", nf );

	//// grid stretching
	//if( sx!=0 ) {
		//pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(nx)/2. );
		//pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(nx)/2. );
	//}
	//if( sy!=0 ) {
		//pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(ny)/2. );
		//pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(ny)/2. );
	//}
	//if( sz!=0 ) {
		//pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
		//pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(nz)/2. );
		//pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(nz)/2. );
	//}

  //// processor grid size
  //pl->set( "npx", npx );
  //pl->set( "npy", npy );
  //pl->set( "npz", npz );
  //pl->set( "npf", npf );

	//auto space = Pimpact::createSpace<ST,OT,4,dNC>( pl );

	//auto x =
		//Pimpact::createMultiField(
				//Pimpact::createCompoundField(
					//Pimpact::createMultiHarmonic<VF>( space ),
					//Pimpact::createMultiHarmonic<SF>( space )
					//)
				//);

	//auto fu   = x->clone();
	//x->init(1.);
	//x->random();

	//auto opV2V =
		//Pimpact::createAdd2Op(
				//Pimpact::createMultiDtHelmholtz( space, 1., 1. ),
				//Pimpact::createMultiHarmonicConvectionOp(space)	);

	//auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
	//auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );

	//auto op =
		//Pimpact::createMultiOperatorBase(
				//Pimpact::createCompoundOpWrap(
					//opV2V,
					//opS2V,
					//opV2S )
				//);

	//// vector to vector operator
	//fu->init(0.);
	//TEST_EQUALITY( 0., fu->norm( Belos::InfNorm ) );

	//opV2V->apply( x->getConstFieldPtr(0)->getConstVField(), fu->getFieldPtr(0)->getVField() );

	//TEST_INEQUALITY( 0., fu->norm( Belos::InfNorm ) );


	//// scalar to vector operator
	//fu->init(0.);
	//TEST_EQUALITY( 0., fu->norm( Belos::InfNorm ) );

	//opS2V->apply( x->getConstFieldPtr(0)->getConstSField(), fu->getFieldPtr(0)->getVField() );

	//TEST_INEQUALITY( 0., fu->norm( Belos::InfNorm ) );

	//// vector to scalar operator
	//fu->init(0.);
	//TEST_EQUALITY( 0., fu->norm( Belos::InfNorm ) );

	//opV2S->apply( x->getConstFieldPtr(0)->getConstVField(), fu->getFieldPtr(0)->getSField() );

	//TEST_INEQUALITY( 0., fu->norm( Belos::InfNorm ) );

	//// compound operator
	//fu->init(0.);
	//TEST_EQUALITY( 0., fu->norm( Belos::InfNorm ) );
	//op->apply(*x,*fu);
	//TEST_INEQUALITY( 0., fu->norm( Belos::InfNorm ) );

//}



TEUCHOS_UNIT_TEST( Convergence, DivOp ) {

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );


  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	//  grid size
	for( int dir=0; dir<3; ++dir ) {

		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );

		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> dofs( ns );

		ST pi2 = std::atan(1)*8;

		if( 0==rank )	
			std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

		for( OT n=0; n<ns; ++n ) {

			if( 0==dir ) 
				pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
			else if( 1==dir )
				pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
			else if( 2==dir )
				pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

			// grid stretching
			if( sx!=0 ) {
				pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx"))/2. );
				pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx"))/2. );
				//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
				//pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
			}
			if( sy!=0 ) {
				pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny"))/2. );
				pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny"))/2. );
			}
			if( sz!=0 ) {
				pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz"))/2. );
				pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz"))/2. );
			}

			auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

			auto vel = Pimpact::create<Pimpact::VectorField>( space );
			auto p   = Pimpact::create<Pimpact::ScalarField>( space );
			auto sol = p->clone( Pimpact::ECopyType::Shallow );

			// init 
			if( 0==dir ) {
				vel->getFieldPtr(Pimpact::U)->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
				sol->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST { return( -std::sin(x*pi2)*pi2 ); } );
			}
			else if( 1==dir ) {
				vel->getFieldPtr(Pimpact::V)->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST { return( -std::sin(y*pi2)*pi2 ); } );
			}
			else if( 2==dir ) {
				vel->getFieldPtr(Pimpact::Z)->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
				sol->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST { return( -std::sin(z*pi2)*pi2 ); } );
			}


			auto op = Pimpact::create<Pimpact::DivOp>( space );

			op->apply( *vel, *p );

			// compute error
			p->add( 1., *sol, -1., *p );
			p->write(n);
			error2[n] = std::log10( p->norm( Belos::TwoNorm ) / sol->norm( Belos::TwoNorm ) );
			errorInf[n] = std::log10( p->norm( Belos::InfNorm ) / sol->norm( Belos::InfNorm ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
			if( 0==rank )	
				std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";

		}
		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << "DivOp: order two norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << "DivOp: order inf norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << orderInf << "\n";
		// test
		TEST_EQUALITY( -order2>3., true );
		TEST_EQUALITY( -orderInf>3., true );
	}
	
}



TEUCHOS_UNIT_TEST( Convergence, InterpolateV2SOp ) { 

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	//  grid size
	for( int dir=0; dir<3; ++dir ) {

		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );

		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> dofs( ns );

		ST pi2 = std::atan(1)*8;

		if( 0==rank )	
			std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

		for( OT n=0; n<ns; ++n ) {

			if( 0==dir ) 
				pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
			else if( 1==dir )
				pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
			else if( 2==dir )
				pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

			// grid stretching
			if( sx!=0 ) {
				pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx"))/2. );
				pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx"))/2. );
				//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
				//pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
			}
			if( sy!=0 ) {
				pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny"))/2. );
				pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny"))/2. );
			}
			if( sz!=0 ) {
				pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz"))/2. );
				pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz"))/2. );
			}

			auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

			auto vel = Pimpact::create<Pimpact::VectorField>( space );
			auto p   = Pimpact::create<Pimpact::ScalarField>( space );
			auto sol = p->clone( Pimpact::ECopyType::Shallow );

			// init 
			if( 0==dir ) {
				vel->getFieldPtr(Pimpact::U)->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
				sol->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
			}
			else if( 1==dir ) {
				vel->getFieldPtr(Pimpact::V)->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
			}
			else if( 2==dir ) {
				vel->getFieldPtr(Pimpact::W)->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
				sol->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
			}

			auto op = Pimpact::createInterpolateV2S( space );

			op->apply( vel->getConstField(dir), *p );

			// compute error
			p->add( 1., *sol, -1., *p );
			error2[n] = std::log10( p->norm( Belos::TwoNorm ) / sol->norm( Belos::TwoNorm ) );
			errorInf[n] = std::log10( p->norm( Belos::InfNorm ) / sol->norm( Belos::InfNorm ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
			if( 0==rank )	
				std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";

		}
		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << "InterpolateV2SOp: order two norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << "InterpolateV2SOp: order inf norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << orderInf << "\n";
		// test
		TEST_EQUALITY( -order2  >4., true );
		TEST_EQUALITY( -orderInf>4., true );
	}
	
}



TEUCHOS_UNIT_TEST( Convergence, InterpolateS2VOp ) { 

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	//  grid size
	for( int dir=0; dir<3; ++dir ) {

		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );


		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> dofs( ns );

		ST pi2 = std::atan(1)*8;

		if( 0==rank )	
			std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

		for( OT n=0; n<ns; ++n ) {

			if( 0==dir ) 
				pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
			else if( 1==dir )
				pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
			else if( 2==dir )
				pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

			// grid stretching
			if( sx!=0 ) {
				pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx"))/2. );
				pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx"))/2. );
				//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
				//pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
			}
			if( sy!=0 ) {
				pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny"))/2. );
				pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny"))/2. );
			}
			if( sz!=0 ) {
				pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz"))/2. );
				pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz"))/2. );
			}

			auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

			auto vel = Pimpact::create<Pimpact::VectorField>( space );
			auto p   = Pimpact::create<Pimpact::ScalarField>( space );
			auto sol = vel->clone( Pimpact::ECopyType::Shallow );

			// init 
			if( 0==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
				sol->getFieldPtr(Pimpact::U)->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
			}
			else if( 1==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol->getFieldPtr(Pimpact::V)->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
			}
			else if( 2==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
				sol->getFieldPtr(Pimpact::W)->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
			}

			auto op = Pimpact::create<Pimpact::InterpolateS2V>( space );

			op->apply(  *p, vel->getField(dir) );

			// compute error
			vel->getFieldPtr(dir)->add( 1., sol->getConstField(dir), -1., vel->getConstField(dir), true );
			vel->write(n);

			error2[n] = std::log10( vel->getConstFieldPtr(dir)->norm( Belos::TwoNorm, true ) / sol->getConstFieldPtr(dir)->norm( Belos::TwoNorm, true ) );
			errorInf[n] = std::log10( vel->getConstFieldPtr(dir)->norm( Belos::InfNorm, true ) / sol->getConstFieldPtr(dir)->norm( Belos::InfNorm, true ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
			if( 0==rank )	
				std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";

		}
		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << "InterpolateS2V: order two norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << "InterpolateS2V: order inf norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << orderInf << "\n";
		// test
		TEST_EQUALITY( -order2  >4., true );
		TEST_EQUALITY( -orderInf>3.9, true );
	}
	
}



TEUCHOS_UNIT_TEST( Convergence, extrapolateBC ) { 

	if( domain!=1 ) {
		//const int dNC = 2;
		Pimpact::setBoundaryConditions( pl, domain );
		pl->set( "dim", dim );

		pl->set( "lx", lx );
		pl->set( "ly", ly );
		pl->set( "lz", lz );

		// processor grid size
		pl->set( "npx", npx );
		pl->set( "npy", npy );
		pl->set( "npz", npz );
		pl->set( "npf", npf );

		//  grid size
		for( int dir=0; dir<3; ++dir ) {

			pl->set<OT>( "nx", 9 );
			pl->set<OT>( "ny", 9 );
			pl->set<OT>( "nz", 9 );
			pl->set<OT>( "nf", 1 );

			std::vector<ST> error2( ns );
			std::vector<ST> errorInf( ns );
			std::vector<ST> dofs( ns );

			ST pi2 = std::atan(1)*8;

			if( 0==rank )	
				std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

			for( OT n=0; n<ns; ++n ) {

				if( 0==dir ) 
					pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
				else if( 1==dir )
					pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
				else if( 2==dir )
					pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

				// grid stretching
				if( sx!=0 ) {
					pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
					pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx"))/2. );
					pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx"))/2. );
					//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
					//pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
				}
				if( sy!=0 ) {
					pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
					pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny"))/2. );
					pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny"))/2. );
				}
				if( sz!=0 ) {
					pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
					pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz"))/2. );
					pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz"))/2. );
				}
				auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

				auto vel = Pimpact::create<Pimpact::VectorField>( space );
				auto sol = vel->clone( Pimpact::ECopyType::Shallow );

				// init 
				if( 0==dir ) {
					vel->getFieldPtr(Pimpact::U)->initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::sin(x*pi2) ); } );
					sol->getFieldPtr(Pimpact::U)->initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::sin(x*pi2) ); } );
				}
				else if( 1==dir ) {
					vel->getFieldPtr(Pimpact::V)->initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return( std::sin(y*pi2) ); } );
					sol->getFieldPtr(Pimpact::V)->initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return( std::sin(y*pi2) ); } );
				}
				else if( 2==dir ) {
					vel->getFieldPtr(Pimpact::W)->initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::sin(z*pi2) ); } );
					sol->getFieldPtr(Pimpact::W)->initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::sin(z*pi2) ); } );
				}
				vel->extrapolateBC();


				// compute error
				vel->getFieldPtr(dir)->add( 1., sol->getConstField(dir), -1., vel->getConstField(dir), true );
				vel->write( n );

				error2[n]   = std::log10( vel->getConstFieldPtr(dir)->norm( Belos::TwoNorm, true ) / sol->getConstFieldPtr(dir)->norm( Belos::TwoNorm, true ) );
				errorInf[n] = std::log10( vel->getConstFieldPtr(dir)->norm( Belos::InfNorm, true ) / sol->getConstFieldPtr(dir)->norm( Belos::InfNorm, true ) );

				dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
				if( 0==rank )	
					std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";
			}
			// compute order
			ST order2 = order<ST>( dofs, error2 );
			if( 0==rank )	
				std::cout << "extrapolateBC: order two norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << order2 << "\n";

			ST orderInf = order<ST>( dofs, errorInf );
			if( 0==rank )	
				std::cout << "extrapolateBC: order inf norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << orderInf << "\n";
			// test
			TEST_EQUALITY( -order2  >3., true );
			TEST_EQUALITY( -orderInf>3., true );
		}
	}
	
}
TEUCHOS_UNIT_TEST( Convergence, GradOp ) { 

	//const int dNC = 2;
  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	//  grid size
	for( int dir=0; dir<3; ++dir ) {

		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );

		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> dofs( ns );

		ST pi2 = std::atan(1)*8;

		if( 0==rank )	
			std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

		for( OT n=0; n<ns; ++n ) {

			if( 0==dir ) 
				pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
			else if( 1==dir )
				pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
			else if( 2==dir )
				pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

			// grid stretching
			if( sx!=0 ) {
				pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx"))/2. );
				pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx"))/2. );
				//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
				//pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
			}
			if( sy!=0 ) {
				pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny"))/2. );
				pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny"))/2. );
			}
			if( sz!=0 ) {
				pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz"))/2. );
				pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz"))/2. );
			}
			auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

			auto vel = Pimpact::create<Pimpact::VectorField>( space );
			auto p   = Pimpact::create<Pimpact::ScalarField>( space );
			auto sol = vel->clone( Pimpact::ECopyType::Shallow );

			// init 
			if( 0==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::cos(x*pi2) ); } );
				sol->getFieldPtr(Pimpact::U)->initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {
							if( (std::abs( y )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
									(std::abs( y-ly )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
							    (std::abs( z )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
									(std::abs( z-lz )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 )  )
								return( -std::sin(x*pi2)*pi2*0.1 ); 
							else
								return( -std::sin(x*pi2)*pi2 );
						} );
			}
			else if( 1==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol->getFieldPtr(Pimpact::V)->initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {
							if( (std::abs( x )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
									(std::abs( x-lx )<Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
							    (std::abs( z )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
									(std::abs( z-lz )<Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 )  )
								return( -std::sin(y*pi2)*pi2*0.1 ); 
							else
								return( -std::sin(y*pi2)*pi2 );
						} );
			}
			else if( 2==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return(  std::cos(z*pi2) ); } );
				sol->getFieldPtr(Pimpact::W)->initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST { 
							if( (std::abs( x )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
									(std::abs( x-lx )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
							    (std::abs( y )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
									(std::abs( y-ly )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 )  )
								return( -std::sin(z*pi2)*pi2*0.1 ); 
							else
								return( -std::sin(z*pi2)*pi2 );
						}  );
			}

			auto op = Pimpact::create<Pimpact::GradOp>( space );

			op->apply(  *p, *vel );

			// compute error
			vel->getFieldPtr(dir)->add( 1., sol->getConstField(dir), -1., vel->getConstField(dir), true );
			vel->write( n );
			
			error2[n] = std::log10( vel->getConstFieldPtr(dir)->norm( Belos::TwoNorm, true ) / sol->getConstFieldPtr(dir)->norm( Belos::TwoNorm, true ) );
			errorInf[n] = std::log10( vel->getConstFieldPtr(dir)->norm( Belos::InfNorm, true ) / sol->getConstFieldPtr(dir)->norm( Belos::InfNorm, true ) );
			//errorInf[n] = std::log10( p->norm( Belos::InfNorm ) / sol->norm( Belos::InfNorm ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
			if( 0==rank )	
				std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";

		}
		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << "GradOp: order two norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << "GradOp: order inf norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << orderInf << "\n";
		// test
		TEST_EQUALITY( -order2  >3., true );
		TEST_EQUALITY( -orderInf>3., true );
	}
	
}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, DivGradOp, OperatorT ) { 

	ST pi2 = std::atan(1)*8;
	std::string label;

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	//  grid size
	for( int dir=0; dir<3; ++dir ) {

		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );

		// grid stretching
		if( sx!=0 ) {
			pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx"))/2. );
			pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx"))/2. );
			//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
			//pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
		}
		if( sy!=0 ) {
			pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny"))/2. );
			pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny"))/2. );
		}
		if( sz!=0 ) {
			pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
			pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz"))/2. );
			pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz"))/2. );
		}


		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> error2Ref( ns );
		std::vector<ST> errorInfRef( ns );
		std::vector<ST> dofs( ns );

		// reference
		Teuchos::RCP<Teuchos::ParameterList> plRef = Teuchos::parameterList();
		*plRef = *pl;

		if( 0==dir ) 
			plRef->set<OT>( "nx", 8*std::pow(2,ns)+1 );
		else if( 1==dir )
			plRef->set<OT>( "ny", 8*std::pow(2,ns)+1 );
		else if( 2==dir )
			plRef->set<OT>( "nz", 8*std::pow(2,ns)+1 );

		// grid stretching
		if( sx!=0 ) {
			plRef->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
			plRef->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(plRef->get<OT>("nx"))/2. );
			plRef->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(plRef->get<OT>("nx"))/2. );
			//plRef->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
			//plRef->sublist("Stretching in X").set<ST>( "x0 U", 0. );
		}
		if( sy!=0 ) {
			plRef->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
			plRef->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(plRef->get<OT>("ny"))/2. );
			plRef->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(plRef->get<OT>("ny"))/2. );
		}
		if( sz!=0 ) {
			plRef->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
			plRef->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(plRef->get<OT>("nz"))/2. );
			plRef->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(plRef->get<OT>("nz"))/2. );
		}

		auto spaceRef = Pimpact::createSpace<ST,OT,d,dNC>( plRef );

		auto yRef = Pimpact::create<Pimpact::ScalarField>( spaceRef );
		auto xRef   = Pimpact::create<Pimpact::ScalarField>( spaceRef );

		auto opRef = Pimpact::create<OperatorT>( spaceRef );
		//opRef->print();
		//auto opRef = Pimpact::create<Pimpact::DivGradOp>( spaceRef );

		// init 
		if( 0==dir ) {
			xRef->initFromFunction(
					[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
		}
		else if( 1==dir ) {
			xRef->initFromFunction(
					[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
		}
		else if( 2==dir ) {
			xRef->initFromFunction(
					[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
		}


		yRef->random();
		opRef->apply( *xRef, *yRef );
		if( 0==rank )	
			std::cout << "\n\n"<<std::setw(5)<< "N\t||e||_2\t\t||e||_inf\n";

		// loop of discretization
		for( OT n=0; n<ns; ++n ) {

			pl->set<OT>( "nx", 9 );
			pl->set<OT>( "ny", 9 );
			pl->set<OT>( "nz", 9 );

			if( 0==dir ) 
				pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
			else if( 1==dir )
				pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
			else if( 2==dir )
				pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

			// grid stretching
			if( sx!=0 ) {
				pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx"))/2. );
				pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx"))/2. );
				//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
				//pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
			}
			if( sy!=0 ) {
				pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny"))/2. );
				pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny"))/2. );
			}
			if( sz!=0 ) {
				pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
				pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz"))/2. );
				pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz"))/2. );
			}

			auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

			auto y = Pimpact::create<Pimpact::ScalarField>( space );
			auto p   = Pimpact::create<Pimpact::ScalarField>( space );
			auto e   = Pimpact::create<Pimpact::ScalarField>( space );
			auto eRef= Pimpact::create<Pimpact::ScalarField>( space );
			auto sol = y->clone( Pimpact::ECopyType::Shallow );

			// init 
			if( 0==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
				sol->initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {  
							if( (std::abs( y )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
									(std::abs( y-ly )<Teuchos::ScalarTraits<ST>::eps()  && space->bcu(1)>0 ) ||
							    (std::abs( z )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
									(std::abs( z-lz )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 )  )
								return( - pi2*pi2*0.1*std::cos(x*pi2) );
							else
								return( - pi2*pi2*std::cos(x*pi2) );
						}
					);
			}
			else if( 1==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol->initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {
							if( (std::abs( x )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
									(std::abs( x-lx )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
							    (std::abs( z )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
									(std::abs( z-lz )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 )  )
								return( - pi2*pi2*0.1*std::cos(y*pi2) );
							else
								return( - pi2*pi2*std::cos(y*pi2) );
						}
					);
			}
			else if( 2==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
				sol->initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {  
							if( (std::abs( x )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
									(std::abs( x-lx )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
							    (std::abs( y )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
									(std::abs( y-ly )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 )  )
								return( - pi2*pi2*0.1*std::cos(z*pi2) );
							else
								return( - pi2*pi2*std::cos(z*pi2) );
						}
					);
			}

			auto op = Pimpact::create<OperatorT>( space );
			label = op->getLabel();

			y->random();
			op->apply( *p, *y );

			// compute error
			e->add( 1., *sol, -1., *y );
			error2[n] = std::log10( e->norm( Belos::TwoNorm ) / sol->norm( Belos::TwoNorm ) );
			errorInf[n] = std::log10( e->norm( Belos::InfNorm ) / sol->norm( Belos::InfNorm ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
			eRef->add( 1., *yRef, -1., *y );
			error2Ref[n] = std::log10( eRef->norm( Belos::TwoNorm ) / std::sqrt( dofs[n] ) );
			errorInfRef[n] = std::log10( eRef->norm( Belos::InfNorm ) );
			if( 0==rank )	
				std::cout  << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\t";
				std::cout                                  << std::pow(10.,error2Ref[n]) << "\t" << std::pow(10.,errorInfRef[n]) << "\n";
			//y->write(0);
			//sol->write(1);
			e->write(2);
			//eRef->write(3);
			//yRef->write(4);
		}

		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << label << ": order two norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << label << ": order inf norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << orderInf << "\n";

		ST order2Ref = order<ST>( dofs, error2Ref );
		if( 0==rank )	
			std::cout << label << ": ref order two norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << order2Ref << "\n";

		ST orderInfRef = order<ST>( dofs, errorInfRef );
		if( 0==rank )	
			std::cout << label << ": ref order inf norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << orderInfRef << "\n";
		// test
		TEST_EQUALITY( -order2  >2., true );
		TEST_EQUALITY( -orderInf>2., true );
	}
	
}

using DivGradOpT = Pimpact::DivGradOp<SpaceT>;
using DivGradO2OpT = Pimpact::DivGradO2Op<SpaceT>;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, DivGradOp, DivGradOpT )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, DivGradOp, DivGradO2OpT )


TEUCHOS_UNIT_TEST( Convergence, HelmholtzOp ) { 

  Pimpact::setBoundaryConditions( pl, domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

	for( int field=0; field<3; ++field ) {
		//  grid size
		for( int dir=0; dir<3; ++dir ) {

			pl->set<OT>( "nx", 9 );
			pl->set<OT>( "ny", 9 );
			pl->set<OT>( "nz", 9 );
			pl->set<OT>( "nf", 1 );

			std::vector<ST> error2( ns );
			std::vector<ST> errorInf( ns );
			std::vector<ST> dofs( ns );

			if( 0==rank )	
				std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

			for( OT n=0; n<ns; ++n ) {

				if( 0==dir ) 
					pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
				else if( 1==dir )
					pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
				else if( 2==dir )
					pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

				// grid stretching
				if( sx!=0 ) {
					pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
					pl->sublist("Stretching in X").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nx")) );
					pl->sublist("Stretching in X").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nx")) );
					//pl->sublist("Stretching in X").set<ST>( "x0 L", 0.05 );
					//pl->sublist("Stretching in X").set<ST>( "x0 U", 0. );
				}
				if( sy!=0 ) {
					pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
					pl->sublist("Stretching in Y").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("ny"))/2. );
					pl->sublist("Stretching in Y").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("ny"))/2. );
				}
				if( sz!=0 ) {
					pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
					pl->sublist("Stretching in Z").set<ST>( "N metr L", static_cast<ST>(pl->get<OT>("nz"))/2. );
					pl->sublist("Stretching in Z").set<ST>( "N metr U", static_cast<ST>(pl->get<OT>("nz"))/2. );
				}

				auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

				auto x   = Pimpact::create<Pimpact::VectorField>( space );
				auto y   = Pimpact::create<Pimpact::VectorField>( space );
				auto sol = x->clone( Pimpact::ECopyType::Shallow );

				// init 
				ST pi2 = std::atan(1)*8;
				ST ire = 1./space->getDomainSize()->getRe();

				x->initField();
				sol->initField();
				if( 0==dir ) {
					x->getFieldPtr(field)->initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST {  return(  std::cos(x*pi2) ); } );
					sol->getFieldPtr(field)->initFromFunction(
							[&pi2,&ire]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2)*pi2*pi2*ire ); } );
				}
				else if( 1==dir ) {
					x->getFieldPtr(field)->initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST {  return(  std::cos(y*pi2) ); } );
					sol->getFieldPtr(field)->initFromFunction(
							[&pi2,&ire]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2)*pi2*pi2*ire ); } );
				}
				else if( 2==dir ) {
					x->getFieldPtr(field)->initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST {  return(  std::cos(z*pi2) ); } );
					sol->getFieldPtr(field)->initFromFunction(
							[&pi2,&ire]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2)*pi2*pi2*ire ); } );
				}

				auto op = Pimpact::create<Pimpact::HelmholtzOp>( space );

				op->apply(  *x, *y );

				// compute error
				y->getFieldPtr(field)->add( 1., sol->getConstField(field), -1., y->getConstField(field) );
				error2[n]   = std::log10( y->getConstFieldPtr(field)->norm( Belos::TwoNorm ) / sol->getConstFieldPtr(field)->norm( Belos::TwoNorm ) );
				errorInf[n] = std::log10( y->getConstFieldPtr(field)->norm( Belos::InfNorm ) / sol->getConstFieldPtr(field)->norm( Belos::InfNorm ) );
				dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
				if( 0==rank )	
					std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";

			}
			// compute order
			ST order2 = order<ST>( dofs, error2 );
			if( 0==rank )	
				std::cout << "Helmholtz(" << Pimpact::toString(static_cast<Pimpact::EField>(field)) << "): order two norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << order2 << "\n";

			ST orderInf = order<ST>( dofs, errorInf );
			if( 0==rank )	
				std::cout << "Helmholtz(" << Pimpact::toString(static_cast<Pimpact::EField>(field)) << "): order inf norm in "<< Pimpact::toString(static_cast<Pimpact::ECoord>(dir)) << "-dir: " << orderInf << "\n";
			// test
			TEST_EQUALITY( -order2  >3., true );
			TEST_EQUALITY( -orderInf>3., true );
		}

	}
}



} // namespace
