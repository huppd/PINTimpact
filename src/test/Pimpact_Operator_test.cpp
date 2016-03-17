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
ST eps = 1e-8;

int dim = 3;
int domain = 0;

ST lx = 1.;
ST ly = 1.;
ST lz = 1.;

ST omega = 0.8;
ST winds = 1;
int sweeps = 1;
int nIter = 1;

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
	clp.setOption( "domain", &domain, "domain" );
	clp.setOption( "dim", &dim, "dim" );
	clp.setOption( "omega", &omega,
      "Slack off of machine epsilon used to check test results" );
	clp.setOption( "wind", &winds,
      "Slack off of machine epsilon used to check test results" );
	clp.setOption( "sweeps", &sweeps, "" );
	clp.setOption( "nIter", &nIter, "" );

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



TEUCHOS_UNIT_TEST( BasicOperator, DivOp ) {

  pl->set( "domain", domain );
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

  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  auto p   = Pimpact::create<Pimpact::ScalarField>( space );
  auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto sol = p->clone( Pimpact::ShallowCopy );

  // zero test
  vel->initField();

  auto op = Pimpact::create<Pimpact::DivOp>( space );

  op->print();

  op->apply( *vel, *p );

  TEST_EQUALITY( p->norm( Belos::InfNorm )<eps, true );


  // random test
  vel->random();

  op->apply(*vel,*p);


  TEST_EQUALITY( p->norm( Belos::InfNorm )>eps, true );

  // circle XY test
  vel->initField( Pimpact::Circle2D );

  op->apply(*vel,*p);

  TEST_EQUALITY( p->norm( Belos::InfNorm )<eps, true );

  // circle XZ test
  vel->initField( Pimpact::Circle2D_inXZ );

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
	if( space->rankST()==0 )	
		std::cout << "GradXYZ error: " << error << "\n";
	TEST_EQUALITY( error>eps, true );

	// Grad test
	vel->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inZ );
	vel->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX );
	vel->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inY );

  p->random();

  op->apply(*vel,*p);

	error = p->norm( Belos::InfNorm );
	if( space->rankST()==0 )	
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

  pl->set( "domain", domain );
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
  //  auto op = Pimpact::create<Pimpact::InterpolateV2S>( space );
	
	if( 0==space->rankST() )
		std::cout << "\n";
	for( int i=0; i<space->dim(); ++i ) {

		for( int bla=1; bla<=6; ++bla ) {

			Pimpact::EScalarField type =	static_cast<Pimpact::EScalarField>(bla);

			vel->getFieldPtr( i )->initField( type );
			sol->initField( type );

			p->random();
			op->apply( vel->getConstField(i), *p );
			sol->add( 1., *sol, -1., *p );

			ST error = sol->norm( Belos::InfNorm );
			if( 0==space->rankST() )
				std::cout << "field " << Pimpact::toString( static_cast<Pimpact::EField>(i) ) << ", error("<<Pimpact::toString(type)<<"): " << error << "\n";
			TEST_EQUALITY( error<eps, true );
			if( error>= eps ){
				std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
				sol->print(
						*Pimpact::createOstream( "error_int"+Pimpact::toString( static_cast<Pimpact::EField>(i) )+"2S_"+Pimpact::toString(type)+"_r"+r+".txt" ));
			}
		}
		
	}

}



TEUCHOS_UNIT_TEST( BasicOperator, InterpolateS2VOp ) {

	pl->set( "domain", domain );
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
	op->print();

	if( 0==space->rankST() )
		std::cout << "\n";

	for( int i=0; i<space->dim(); ++i ) {

		auto sol = solv->getFieldPtr( i );

		for( int bla=1; bla<=6; ++bla ) {

			Pimpact::EScalarField type =	static_cast<Pimpact::EScalarField>(bla);

			vel->getFieldPtr( i )->random();
			sol->initField( type );
			p->initField( type );

			op->apply( *p, vel->getField( i ) );
			sol->add( 1., *sol, -1., vel->getConstField(i) );

			ST error = sol->norm( Belos::InfNorm );
			if( 0==space->rankST() )
				std::cout << "field "
					<< Pimpact::toString( static_cast<Pimpact::EField>(i) )
					<< ", error("<<Pimpact::toString(type)<<"): "
					<< error << "\n";

			TEST_EQUALITY( error<eps, true );
			if( error>=eps ) {
				std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
				sol->print(
						*Pimpact::createOstream( "error_int_S2"+Pimpact::toString( static_cast<Pimpact::EField>(i) )+"_"+Pimpact::toString(type)+"_r"+r+".txt" ));
			}

		}

	}

}



TEUCHOS_UNIT_TEST( BasicOperator, TransferOp ) {

  pl->set( "domain", domain );
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

  pl->set( "domain", domain );
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
  fx->initField( Pimpact::PoiseuilleFlow2D_inX );
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

  pl->set( "domain", domain );
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
  op->print();

  // op in x test
  p->initField( Pimpact::Grad2D_inX );
	sol->initField( Pimpact::ConstFlow, 1., 0., 0. );
  v->random();

  op->apply(*p,*v);

	sol->add( 1., *sol, -1., *v );

	ST error = sol->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error: " << error << "\n";
  TEST_EQUALITY( error<eps, true );

  // op in y test
  p->initField( Pimpact::Grad2D_inY );
	sol->initField( Pimpact::ConstFlow, 0., 1., 0. );
  v->random();

  op->apply(*p,*v);

	sol->add( 1., *sol, -1., *v );

	error = sol->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error: " << error << "\n";
  TEST_EQUALITY( error<eps, true );

  // op in z test
  p->initField( Pimpact::Grad2D_inZ );
	sol->initField( Pimpact::ConstFlow, 0., 0., 1. );
  v->random();

  op->apply(*p,*v);

	sol->add( 1., *sol, -1., *v );

	error = sol->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error: " << error << "\n";
  TEST_EQUALITY( error<eps, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, HelmholtzOp ) {

  pl->set( "domain", domain );
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

  op->print();

  // test in x direction
  x->initField( Pimpact::PoiseuilleFlow2D_inX );
  bs->init( Teuchos::tuple( 8./std::pow(space->getDomainSize()->getSize(Pimpact::Y),2), 0., 0. ) );
  bs->add( mulI, *x, mulL, *bs );

	auto para = Teuchos::parameterList();
	para->set<ST>( "mulI", mulI );
	op->setParameter( para );
  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);

	ST error = bs->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error(poiseuilleX): " << error << "\n";
  TEST_EQUALITY( error<eps, true );
	if( error>= eps ){
		std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
		bs->print(
				*Pimpact::createOstream( "error_helm_poiX_r"+r+".txt" ));
	}

  // test in y direction
  x->initField( Pimpact::PoiseuilleFlow2D_inY );
  bs->init( Teuchos::tuple( 0., 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2), 0. ) );
  bs->add( mulI, *x, mulL, *bs );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b );

	error = bs->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error(poiseuilleY): " << error << "\n";
  TEST_EQUALITY( error<eps, true );
	if( error>= eps ){
		std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
		bs->print(
				*Pimpact::createOstream( "error_helm_poiY_r"+r+".txt" ));
	}

  // test in z direction
  x->initField( Pimpact::PoiseuilleFlow2D_inZ );
  bs->init( Teuchos::tuple( 0., 0., 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2) ) );
  bs->add( mulI, *x, mulL, *bs );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);

	error = bs->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error(poiseuilleZ): " << error << "\n";
  TEST_EQUALITY( error<eps, true );
	if( error>= eps ){
		std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
		bs->print(
				*Pimpact::createOstream( "error_helm_poiZ_r"+r+".txt" ));
	}

  // the circle XY test
  x->initField( Pimpact::Circle2D );
  bs->initField( Pimpact::Circle2D );
  bs->scale( mulI );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b );

	error = bs->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error(circleXY): " << error << "\n";
  TEST_EQUALITY( error<eps, true );
	if( error>= eps ){
		std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
		bs->print(
				*Pimpact::createOstream( "error_helm_circXY_r"+r+".txt" ));
	}

  // the circle XZ test
  x->initField( Pimpact::Circle2D_inXZ );
  bs->initField( Pimpact::Circle2D_inXZ );
  bs->scale( mulI );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b );

	error = bs->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error(circleXZ): " << error << "\n";
  TEST_EQUALITY( error<eps, true );
	if( error>= eps ){
		std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
		bs->print(
				*Pimpact::createOstream( "error_helm_circXZ_r"+r+".txt" ));
	}

}



TEUCHOS_UNIT_TEST( BasicOperator, DivGradO2Op ) {

	const int dNC = 2;

  pl->set( "domain", domain );
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

  op->print();

  // zero test
  b->initField( Pimpact::ConstField, 0. );
  x->random();

  op->apply( *b, *x );

	ST error = x->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error(0): " << error << "\n";
  TEST_EQUALITY( error<eps, true );

  b->initField( Pimpact::ConstField, 0. );
  x->random();

  op2->apply( *b, *x );

	error = x->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error(0) O2: " << error << "\n";
  TEST_EQUALITY( error<eps, true );

	// one test
  b->initField( Pimpact::ConstField, 1. );
  x->random();

  op->apply( *b, *x );

	error = x->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error(1): " << error << "\n";
	TEST_EQUALITY( error<eps, true );

  b->initField( Pimpact::ConstField, 1. );
  x->random();

  op2->apply( *b, *x );

	error = x->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error(1) O2: " << error << "\n";
	if( error>= eps ) {
		x->print();
	}
	TEST_EQUALITY( error<eps, true );

	// consistency test
	b->random();
	op->apply( *b, *x );
	op2->apply( *b, *x2 );

	x2->add( 1., *x2, -1., *x );
	ST errInf = x2->norm( Belos::InfNorm );
	ST err2 = x2->norm( Belos::TwoNorm )/std::sqrt( x2->getLength() );
	if( 0==space->rankST() )
		std::cout << "consistency error random: inf: " << errInf << ", two: " << err2 << "\n";
	TEST_EQUALITY( errInf<eps, true );
	TEST_EQUALITY( err2<eps, true );

	for( int dir=0; dir<=6; ++dir ) {

		Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>( dir );

		b->initField( type );
		op->apply( *b, *x );
		op2->apply( *b, *x2 );

		x2->add( 1., *x2, -1., *x );
		errInf = x2->norm( Belos::InfNorm );
		err2 = x2->norm( Belos::TwoNorm )/std::sqrt( x2->getLength() );
		if( 0==space->rankST() )
			std::cout << "consistency error("+Pimpact::toString(type)+"): inf: " << errInf << ", two: " << err2 << "\n";
		TEST_EQUALITY( errInf<eps, true );
		TEST_EQUALITY( err2<eps, true );
		if( errInf>=eps && err2>=eps ) {
			std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
			x2->print(
					*Pimpact::createOstream( "error_dgo2_"+Pimpact::toString(type)+"_r"+r+".txt" ));
		}

	}

}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, DivGradO2Smoother, SType ) {

  pl->set( "domain", domain );
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
	if( 0==space->rankST() ) {
		Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "coord.txt" );
		coord->print( *fstream );
	}

	Teuchos::RCP< Pimpact::ScalarField< Pimpact::Space<ST,OT,d,dNC> > > x =
		Pimpact::create<Pimpact::ScalarField>( space );

  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );
	op->print2Mat();

	// compute EV
	ST evMax;
	ST evMin;


	Teuchos::RCP<Pimpact::TeuchosEigenvalues<Pimpact::DivGradO2Op<SpaceT> > > ev = 
		Teuchos::rcp( new Pimpact::TeuchosEigenvalues<Pimpact::DivGradO2Op<SpaceT> >( op ) );

	std::cout << "\n";
	for( int i=0; i<3; ++i ) {
		ev->computeEV( static_cast<Pimpact::ECoord>(i), evMax, evMin );
		std::cout << Pimpact::toString( static_cast<Pimpact::ECoord>(i) ) << ": " << evMax << "\t" <<evMin << "\n";
	}

	ev->computeEV( evMax, evMin );
	std::cout << "glob : " << evMax << "\t" <<evMin << "\n";

	if( nx<=20 && ny<=20 && nz<=20 ) {
		ev->computeFullEV( evMax, evMin );
		std::cout << "glob2: " << evMax << "\t" <<evMin << "\n";
	}

	Teuchos::RCP<Teuchos::ParameterList> ppl = Teuchos::parameterList();
	ppl->set<int>( "numIters", sweeps );

	// good result with no stretch + own
	ppl->set<ST>( "max EV", evMin*1.1 );
	ppl->set<ST>( "min EV", evMax );
	//ppl->set<ST>( "min EV", evMin*1.1/30. );
	//ppl->set<ST>( "max EV", evMin*1.1 );
	//ppl->set<ST>( "min EV", evMin*1.1 );

	ppl->set<bool>( "with output", true );
	
	// JSmoother
	ppl->set<int>( "BC smoothing", 1 );
	ppl->set<OT>( "depth", 2 );

  auto smoother = Pimpact::create<SType>( op, ppl );

	// --- zero rhs test --- 
	x->random();
	b->init( 0. );

	//Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "conv_"+smoother->getLabel()+".txt", 0 );
	if( 0==space->rankST() ) {
		std::cout << "\nstep\terror\t\trate\n";
		std::cout << 0 << "\t" << 1. << "\n";
	}
	ST err0 = x->norm( Belos::InfNorm );
	ST errP = err0;

	//x->write( 0 );
	for( int i=1; i<=nIter; ++i ) {
		smoother->apply( *b, *x );
		//x->write(i);
		ST err = x->norm( Belos::InfNorm );
		if( 0==space->rankST() ) {
				std::cout << i << "\t" << err/err0 << "\t" << err/errP << "\n";
				//*fstream << i << "\t" << err/err0 << "\t" << err/errP << "\n";
		}
		errP = err;
	}


	// --- consistency test ---
	for( int dir=1; dir<=1; ++dir ) {

		x->initField( static_cast<Pimpact::EScalarField>(dir) );
		x->level();
		auto xp = x->clone( Pimpact::DeepCopy );

		op->apply( *x, *b );

		smoother->apply( *b, *x );
		x->level();

		xp->add( 1., *xp, -1., *x );
		ST err2 = xp->norm( Belos::InfNorm )/std::sqrt( static_cast<ST>(xp->getLength()) );
		ST errInf = xp->norm( Belos::InfNorm );

		if( 0==space->rankST() ) {
			std::cout << "consistency for " << dir << ": ||" << err2 << "||_2, ||" << errInf << "||_inf\n";
		}

		TEST_EQUALITY( err2<eps, true );
		TEST_EQUALITY( errInf<eps, true );
		//if( err2>eps && errInf>eps ) 
			//xp->write( dir );

	}

}

using JT = Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op<SpaceT> >;
using SORT = Pimpact::DivGradO2SORSmoother< Pimpact::DivGradO2Op<SpaceT> >;
using CheT = Pimpact::Chebyshev< Pimpact::DivGradO2Op<SpaceT> >;
using LT = Pimpact::DivGradO2LSmoother< Pimpact::DivGradO2Op<SpaceT> >;

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, JT )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, SORT )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, CheT )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, LT )



TEUCHOS_UNIT_TEST( BasicOperator, DivGradO2Inv ) {

  pl->set( "domain", domain );
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

	op->print2Mat();

  auto ppl = Teuchos::parameterList();

  auto solver = Pimpact::create<Pimpact::DivGradO2Inv>( op, ppl );

  // --- consistency test ---
	//for( int dir=1; dir<=6; ++dir ) {
	for( int dir=1; dir<=6; ++dir ) {
		x->initField( static_cast<Pimpact::EScalarField>(dir) );
		x->level();
		xp->random();

		op->apply( *x, *b );
		solver->apply( *b, *xp );

		xp->add( 1., *xp, -1., *x );
		//xp->level();

		ST err2 = xp->norm( Belos::InfNorm )/std::sqrt( static_cast<ST>(xp->getLength()) );
		ST errInf = xp->norm( Belos::InfNorm );
		//if( xp->norm( Belos::InfNorm, false )>=eps ) {
			//std::cout << "rank: " << space->rankST() << "\te: " << xp->norm( Belos::InfNorm, false ) << "\n";
			//if( 0==space->rankST() )
				//xp->print();
		//}


		if( 0==space->rankST() )
			std::cout << "consistency for " << dir << ": ||" << err2 << "||_2, ||" << errInf << "||_inf\n";

		TEST_EQUALITY( err2<eps, true );
		TEST_EQUALITY( errInf<eps, true );
	}

}



TEUCHOS_UNIT_TEST( BasicOperator, ForcingOp ) {

  pl->set( "domain", domain );
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
	force->initField( Pimpact::EVectorField(11) );
	res->random();

	auto op = Pimpact::createForcingOp( force, 1. );

	op->apply( *vel, *res );
	vel->add( 1., *res, -1., *force );

	TEST_EQUALITY( vel->norm( Belos::InfNorm )<eps, true );
}





template<class T> using MSFT = Pimpact::MultiField< Pimpact::ScalarField<T> >;
template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;

TEUCHOS_UNIT_TEST( MultiOperator, InverseOp ) {

  pl->set( "domain", domain );
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

  auto x = Pimpact::create<Pimpact::ScalarField>( space );
  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );

  auto opinv = Pimpact::create<MOP>( op );

  // Grad in x
  x->initField( Pimpact::Grad2D_inX );

  op->apply(*x,*b);

  x->init();

  opinv->apply( *b, *x );

  b->initField( Pimpact::Grad2D_inX );
  x->add( -1, *x, 1., *b );

	ST error = x->norm( Belos::InfNorm )/b->norm( Belos::InfNorm );
	//if( 0==space->rankST() )
		//std::cout << "\nerror: " << error << "\n";
	//TEST_EQUALITY( error<0.5, true );

  // Grad in y
  x->initField( Pimpact::Grad2D_inY );

  op->apply(*x,*b);

  x->init();

  opinv->apply( *b, *x );

  b->initField( Pimpact::Grad2D_inY );
  x->add( -1, *x, 1., *b );

	error = x->norm( Belos::InfNorm )/b->norm( Belos::InfNorm );
	if( 0==space->rankST() )
		std::cout << "error: " << error << "\n";
	TEST_EQUALITY( error<0.5, true );

}



TEUCHOS_UNIT_TEST( MultiOperator, Add2Op ) {

  pl->set( "domain", domain );
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


 for( OT i=0; i<10; ++i ) {
	 x->getFieldPtr(i)->initField(Pimpact::RankineVortex2D );
 }

 op->apply( *x, *y);

}


template<class T> using WUP = Pimpact::MultiOpUnWrap<Pimpact::MultiOpWrap<T> >;


TEUCHOS_UNIT_TEST( MultiOperator, MulitOpUnWrap ) {

  pl->set( "domain", domain );
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


  x->initField(Pimpact::RankineVortex2D );


  opUW2->apply( *x, *y);

}



TEUCHOS_UNIT_TEST( MultiModeOperator, HelmholtzOp ) {

  pl->set( "domain", domain );
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

  pl->set( "domain", domain );
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

	//pl->set( "domain", domain );
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

	//auto diff = mv->getConstFieldPtr(0)->getConstCFieldPtr()->clone(Pimpact::ShallowCopy);

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




TEUCHOS_UNIT_TEST( MultiModeOperator, TripleCompostion ) {

  pl->set( "domain", domain );
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

	pl->set( "Re", 1. );
	pl->set( "alpha2", 1. );
  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  using MVF = Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > >;
  using MSF = Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<SpaceT> > >;

  auto X = Pimpact::create<MSF>( space );
  auto B = Pimpact::create<MSF>( space );

  auto temp = Pimpact::create<MVF>( space );

  X->init(0.);
//  B->random();
  B->init(1.);


  auto H =
      Pimpact::createMultiOperatorBase(
          Pimpact::createDtLapOp( space, 0., 10. ) );

  H->apply( *temp, *temp );

  // Make an empty new parameter list.
  //auto solverParams = Teuchos::parameterList();
  auto solverParams = Pimpact::createLinSolverParameter( "GMRES", 1.e-1 );

  // Create the Pimpact::LinearSolver solver.
  auto Hprob = Pimpact::createLinearProblem<MVF>( H, temp, temp, solverParams,"GMRES" );
  auto Hinv  = Pimpact::createInverseOperator( Hprob );

  auto schur = Pimpact::createTripleCompositionOp(
      Pimpact::createMultiModeOpWrap( Pimpact::create<Pimpact::DivOp>( space ) ),
      Hinv,
      Pimpact::createMultiModeOpWrap( Pimpact::create<Pimpact::GradOp>( space ) )
  );

  schur->apply( *B, *X );

}




TEUCHOS_UNIT_TEST( MultiModeOperator, InverseOperator ) {

	pl->set( "domain", domain );
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


	using MVF = Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > >;

	auto X = Pimpact::create<MVF>( space );
	auto B = Pimpact::create<MVF>( space );

	X->init(0.);
	B->random();

	auto op = Pimpact::createMultiModeOperatorBase<MVF>(
			Pimpact::create<Pimpact::HelmholtzOp>(space) );

	// Make an empty new parameter list.
	auto solverParams = Pimpact::createLinSolverParameter( "CG", 1.e-1 );

	// Create the Pimpact::LinearSolver solver.
	auto prob =
		Pimpact::createLinearProblem<MVF>(
				op,
				X,X,
				//          Pimpact::createMultiField(X->getFieldPtr(0)->getCFieldPtr()),
				//          Pimpact::createMultiField(B->getFieldPtr(0)->getCFieldPtr()),
				solverParams );

	auto opinv = Pimpact::createInverseOperator<MVF>( prob );
	auto opp = Pimpact::createOperatorBase( opinv );

}





TEUCHOS_UNIT_TEST( MultiModeOperator, EddyPrec ) {

	pl->set( "domain", domain );
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

	//pl->set("alpha2",1.);
	//pl->set("Re",1./14.);
	auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

	auto temp = Pimpact::createMultiModeVectorField( space );

	auto X = Pimpact::createMultiModeVectorField( space );
	auto B = Pimpact::createMultiModeVectorField( space );

	X->init(0.);
	B->init(1.);

	// Make an empty new parameter list.
	auto solverParams = Pimpact::createLinSolverParameter("CG",1.e-1);

	// Create the Pimpact::LinearSolver solver.
	auto A =
		Pimpact::createMultiOperatorBase(
				Pimpact::create<Pimpact::HelmholtzOp>( space )
				);

	A->apply( *Pimpact::createMultiField( B->getFieldPtr(0)->getCFieldPtr() ), 
			*Pimpact::createMultiField( X->getFieldPtr(0)->getCFieldPtr() ) );

	auto prob =
		Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::VectorField<SpaceT> > >(
				A,
				Teuchos::null,
				Teuchos::null,
				solverParams,
				"CG" );

	prob->solve( Pimpact::createMultiField( B->getFieldPtr(0)->getCFieldPtr() ), 
			Pimpact::createMultiField( X->getFieldPtr(0)->getCFieldPtr() ) );
	
	auto invOp = Pimpact::createInverseOperator( prob );

	invOp->apply( Pimpact::createMultiField( B->getFieldPtr(0)->getCFieldPtr() ), 
			Pimpact::createMultiField( X->getFieldPtr(0)->getCFieldPtr() ) );

	auto op =
		Pimpact::create<Pimpact::EddyPrec>(
				Pimpact::createMultiOpUnWrap(
					Pimpact::createInverseOperatorBase(prob) ) );
//	auto op =
//		Pimpact::create<Pimpact::EddyPrec>( 
//				Pimpact::createMultiOpUnWrap( A )	);
	
	op->apply( X->getField(0), B->getField(0) );
	
	auto schur = Pimpact::createMultiOperatorBase( op );

	schur->apply( *B, *X );

}




TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicConvectionOp ) {

	pl->set( "domain", domain );
	pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	pl->set<bool>( "spectral in time", true );
  pl->set( "domain", domain );
  pl->set( "dim", dim );

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

  auto op = Pimpact::createMultiHarmonicConvectionOp( space );

  op->assignField( *mv1 );
  op->apply( *mv1, *mv2 );

}



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicOpWrap ) {

	using SpaceT = Pimpact::Space<ST,OT,4,dNC>;

	pl->set( "domain", domain );
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



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiDtHelmholtz) {

	pl->set( "domain", domain );
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

	auto mv1 = Pimpact::createMultiHarmonicVectorField( space );
	auto mv2 = Pimpact::createMultiHarmonicVectorField( space );

	auto op = Pimpact::createMultiDtHelmholtz( space );

	op->apply( *mv1, *mv2 );

}



TEUCHOS_UNIT_TEST( CompoundOperator, CompoundOpWrap ) {

	using SpaceT = Pimpact::Space<ST,OT,4,dNC>;

	using VF = Pimpact::VectorField<SpaceT>;
	using SF = Pimpact::ScalarField<SpaceT>;

	pl->set( "domain", domain );
	pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	pl->set<bool>( "spectral in time", true );

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

	auto x =
		Pimpact::createMultiField(
				Pimpact::createCompoundField(
					Pimpact::createMultiHarmonic<VF>( space ),
					Pimpact::createMultiHarmonic<SF>( space )
					)
				);

	auto fu   = x->clone();
	x->init(1.);
	x->random();

	auto opV2V =
		Pimpact::createAdd2Op(
				Pimpact::createMultiDtHelmholtz( space, 1., 1. ),
				Pimpact::createMultiHarmonicConvectionOp(space)	);

	auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
	auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );

	auto op =
		Pimpact::createMultiOperatorBase(
				Pimpact::createCompoundOpWrap(
					opV2V,
					opS2V,
					opV2S )
				);

	// vector to vector operator
	fu->init(0.);
	TEST_EQUALITY( 0., fu->norm( Belos::InfNorm ) );

	opV2V->apply( x->getConstFieldPtr(0)->getConstVField(), fu->getFieldPtr(0)->getVField() );

	TEST_INEQUALITY( 0., fu->norm( Belos::InfNorm ) );


	// scalar to vector operator
	fu->init(0.);
	TEST_EQUALITY( 0., fu->norm( Belos::InfNorm ) );

	opS2V->apply( x->getConstFieldPtr(0)->getConstSField(), fu->getFieldPtr(0)->getVField() );

	TEST_INEQUALITY( 0., fu->norm( Belos::InfNorm ) );

	// vector to scalar operator
	fu->init(0.);
	TEST_EQUALITY( 0., fu->norm( Belos::InfNorm ) );

	opV2S->apply( x->getConstFieldPtr(0)->getConstVField(), fu->getFieldPtr(0)->getSField() );

	TEST_INEQUALITY( 0., fu->norm( Belos::InfNorm ) );

	// compound operator
	fu->init(0.);
	TEST_EQUALITY( 0., fu->norm( Belos::InfNorm ) );
	op->apply(*x,*fu);
	TEST_INEQUALITY( 0., fu->norm( Belos::InfNorm ) );

}


} // namespace
