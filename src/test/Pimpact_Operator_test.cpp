#include <iostream>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Fields.hpp"

#include "Pimpact_TransferOp.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_VectorFieldOpWrap.hpp"

#include "Pimpact_LinSolverParameter.hpp"




namespace {


using S = double;
using O = int;
const int d = 3;
const int dNC = 4;
//const int dNC = 2;

using SpaceT = Pimpact::Space<S,O,d,dNC>;

using SF = typename Pimpact::ScalarField<SpaceT>;
using VF = typename Pimpact::VectorField<SpaceT>;
using MSF = typename Pimpact::ModeField<SF>;
using MVF = typename Pimpact::ModeField<VF>;



bool testMpi = true;
S eps = 1e-8;

int dim = 3;
int domain = 0;

S omega = 0.8;
S winds = 1;
int nIter = 1000;

int nx = 97;
int ny = 25;
int nz = 49;
int nf = 32;


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
	clp.setOption( "nIter", &nIter,
      "Slack off of machine epsilon used to check test results" );


	clp.setOption( "nx", &nx, "" );
	clp.setOption( "ny", &ny, "" );
	clp.setOption( "nz", &nz, "" );
	clp.setOption( "nf", &nf, "" );

	pl->set( "lx", 30. );
	pl->set( "ly", 20. );
	pl->set( "lz", 20. );

	//pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
	//pl->sublist("Stretching in X").set<S>( "N metr L", nx );
	////pl->sublist("Stretching in X").set<S>( "N metr L", nx/2 );
	//pl->sublist("Stretching in X").set<S>( "N metr U", nx  );
	//pl->sublist("Stretching in X").set<S>( "x0 L", 0.05 );
	////pl->sublist("Stretching in X").set<S>( "x0 L", 0. );
	////pl->sublist("Stretching in X").set<S>( "x0 U", 0. );

	//pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
	//pl->sublist("Stretching in Y").set<S>( "N metr L", 1 );
	//pl->sublist("Stretching in Y").set<S>( "N metr U", ny  );
	//pl->sublist("Stretching in Y").set<S>( "x0 L", 0. );
	//pl->sublist("Stretching in Y").set<S>( "x0 U", 0. );

	//pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
	//pl->sublist("Stretching in Z").set<S>( "N metr L", 1 );
	//pl->sublist("Stretching in Z").set<S>( "N metr U", nz  );
	//pl->sublist("Stretching in Z").set<S>( "x0 L", 0. );
	//pl->sublist("Stretching in Z").set<S>( "x0 U", 0. );
}



TEUCHOS_UNIT_TEST( BasicOperator, DivOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set( "nx", nx );
	pl->set( "ny", ny );
	pl->set( "nz", nz );
	pl->set( "nf", nf );

  // processor grid size
  pl->set( "npx", (2==dim)?4:2 );
  pl->set( "npy",            2 );
  pl->set( "npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p   = Pimpact::create<Pimpact::ScalarField>( space );
  auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto sol = p->clone( Pimpact::ShallowCopy );

  // zero test
  vel->initField();

  auto op = Pimpact::create<Pimpact::DivOp>( space );

  op->print();

  op->apply(*vel,*p);

  TEST_EQUALITY( p->norm()<eps, true );


  // random test
  vel->random();

  op->apply(*vel,*p);


  TEST_EQUALITY( p->norm()>eps, true );

  // circle XY test
  vel->initField( Pimpact::Circle2D );

  op->apply(*vel,*p);

  TEST_EQUALITY( p->norm()<eps, true );

  // circle XZ test
  vel->initField( Pimpact::Circle2D_inXZ );

  op->apply(*vel,*p);

  TEST_EQUALITY( p->norm()<eps, true );

	// Grad test
	vel->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
	vel->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY, -2. );
	vel->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );

  p->random();

  op->apply(*vel,*p);

	sol->initField( Pimpact::ConstField, space->dim() );
	sol->add( 1., *sol, -1., *p );

	S error = sol->norm();
	if( space->rankST()==0 )	
		std::cout << "GradXYZ error: " << error << "\n";
	TEST_EQUALITY( error>eps, true );

	// Grad test
	vel->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inZ );
	vel->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX );
	vel->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inY );

  p->random();

  op->apply(*vel,*p);

	error = p->norm();
	if( space->rankST()==0 )	
		std::cout << "GradZXY error: " << error << "\n";
	TEST_EQUALITY( error<eps, true );

	// Grad test
	vel->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY );
	vel->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inZ );
	vel->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX );

  p->random();

  op->apply(*vel,*p);

	std::cout << "erro: " << p->norm() << "\n";
	TEST_EQUALITY( p->norm()<eps, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, InterpolateV2SOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<Pimpact::ScalarField>( space );
	auto sol = p->clone();
  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto op = Pimpact::createInterpolateV2S( space );
  //  auto op = Pimpact::create<Pimpact::InterpolateV2S>( space );

	for( int i=0; i<space->dim(); ++i ) {

		// X test
		vel->getFieldPtr( i )->initField( Pimpact::Grad2D_inX );
		sol->initField( Pimpact::Grad2D_inX );

		p->random();
		op->apply( vel->getConstField( i ), *p );
		sol->add( 1., *sol, -1., *p );

		std::cout << "dir: " << i << " GradX error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		vel->getFieldPtr( i )->initField( Pimpact::Poiseuille2D_inX );
		sol->initField( Pimpact::Poiseuille2D_inX );

		p->random();
		op->apply( vel->getConstField( i ), *p );
		sol->add( 1., *sol, -1., *p );

		std::cout << "dir: " << i << " PoiseuilleX error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		// Y test
		vel->getFieldPtr( i )->initField( Pimpact::Grad2D_inY );
		sol->initField( Pimpact::Grad2D_inY );

		p->random();
		op->apply( vel->getConstField( i ), *p );
		sol->add( 1., *sol, -1., *p );

		std::cout << "dir: " << i << " GradY error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		vel->getFieldPtr( i )->initField( Pimpact::Poiseuille2D_inY );
		sol->initField( Pimpact::Poiseuille2D_inY );

		p->random();
		op->apply( vel->getConstField( i ), *p );
		sol->add( 1., *sol, -1., *p );

		std::cout << "dir: " << i << " PoiseuilleY error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		// Z test
		vel->getFieldPtr( i )->initField( Pimpact::Grad2D_inZ );
		sol->initField( Pimpact::Grad2D_inZ );

		p->random();
		op->apply( vel->getConstField( i ), *p );
		sol->add( 1., *sol, -1., *p );

		std::cout << "dir: " << i << " GradZ error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		vel->getFieldPtr( i )->initField( Pimpact::Poiseuille2D_inZ );
		sol->initField( Pimpact::Poiseuille2D_inZ );

		p->random();
		op->apply( vel->getConstField( i ), *p );
		sol->add( 1., *sol, -1., *p );

		std::cout << "dir: " << i << " PoiseuilleZ error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

	}

}



TEUCHOS_UNIT_TEST( BasicOperator, InterpolateS2VOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<Pimpact::ScalarField>( space );
  auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto solv = vel->clone();

  auto op = Pimpact::create<Pimpact::InterpolateS2V>( space );
  op->print();


	for( int i=0; i<space->dim(); ++i ) {

		auto sol = solv->getFieldPtr( i );

		// X test
		vel->getFieldPtr( i )->random();
		sol->initField( Pimpact::Grad2D_inX );
		p->initField( Pimpact::Grad2D_inX );

		op->apply( *p, vel->getField( i ) );
		sol->add( 1., *sol, -1., vel->getConstField(i) );

		std::cout << "dir: " << i << " GradX error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		vel->getFieldPtr( i )->random();
		sol->initField( Pimpact::Poiseuille2D_inX );
		p->initField( Pimpact::Poiseuille2D_inX );

		op->apply( *p, vel->getField( i ) );
		sol->add( 1., *sol, -1., vel->getConstField(i) );

		std::cout << "dir: " << i << " PoiseuilleX error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		// Y test
		vel->getFieldPtr( i )->random();
		sol->initField( Pimpact::Grad2D_inY );
		p->initField( Pimpact::Grad2D_inY );

		op->apply( *p, vel->getField( i ) );
		sol->add( 1., *sol, -1., vel->getConstField(i) );

		std::cout << "dir: " << i << " GradY error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		vel->getFieldPtr( i )->random();
		sol->initField( Pimpact::Poiseuille2D_inY );
		p->initField( Pimpact::Poiseuille2D_inY );

		op->apply( *p, vel->getField( i ) );
		sol->add( 1., *sol, -1., vel->getConstField(i) );

		std::cout << "dir: " << i << " PoiseuilleY error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		// Z test
		vel->getFieldPtr( i )->random();
		sol->initField( Pimpact::Grad2D_inZ );
		p->initField( Pimpact::Grad2D_inZ );

		op->apply( *p, vel->getField( i ) );
		sol->add( 1., *sol, -1., vel->getConstField(i) );

		std::cout << "dir: " << i << " GradZ error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

		vel->getFieldPtr( i )->random();
		sol->initField( Pimpact::Poiseuille2D_inZ );
		p->initField( Pimpact::Poiseuille2D_inZ );

		op->apply( *p, vel->getField( i ) );
		sol->add( 1., *sol, -1., vel->getConstField(i) );

		std::cout << "dir: " << i << " PoiseuilleZ error: " << sol->norm() << "\n";
		TEST_EQUALITY( sol->norm()<eps, true );

	}

}



TEUCHOS_UNIT_TEST( BasicOperator, TransferOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  typedef Pimpact::Space<S,O,d,4> FSpaceT;
  typedef Pimpact::Space<S,O,d,2> CSpaceT;

  auto fSpace = Pimpact::createSpace<S,O,d,4>( pl );
  auto cSpace = Pimpact::createSpace<S,O,d,2>( pl );

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

//  cx->write(0);

  fx->random();

  op->apply( *cx, *fx );

  TEST_FLOATING_EQUALITY( fx->norm(Belos::OneNorm), cx->norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::TwoNorm), cx->norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::InfNorm), cx->norm(Belos::InfNorm), eps );

//  fx->write(1);

}



TEUCHOS_UNIT_TEST( BasicOperator, VectorFieldOpWrap ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );
  typedef Pimpact::Space<S,O,d,4> FSpaceT;
  typedef Pimpact::Space<S,O,d,2> CSpaceT;

  auto fSpace = Pimpact::createSpace<S,O,d,4>( pl );
  auto cSpace = Pimpact::createSpace<S,O,d,2>( pl );

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

//  cx->write(0);

  fx->random();

  op->apply( *cx, *fx );

  TEST_FLOATING_EQUALITY( fx->norm(Belos::OneNorm), cx->norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::TwoNorm), cx->norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::InfNorm), cx->norm(Belos::InfNorm), eps );

//  fx->write(1);

}



TEUCHOS_UNIT_TEST( BasicOperator, GradOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );


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

	std::cout << "error: " << sol->norm() << "\n";
  TEST_EQUALITY( sol->norm()<eps, true );

  // op in y test
  p->initField( Pimpact::Grad2D_inY );
	sol->initField( Pimpact::ConstFlow, 0., 1., 0. );
  v->random();

  op->apply(*p,*v);

	sol->add( 1., *sol, -1., *v );

	std::cout << "error: " << sol->norm() << "\n";
  TEST_EQUALITY( sol->norm()<eps, true );

  // op in z test
  p->initField( Pimpact::Grad2D_inZ );
	sol->initField( Pimpact::ConstFlow, 0., 0., 1. );
  v->random();

  op->apply(*p,*v);

	sol->add( 1., *sol, -1., *v );

	std::cout << "error: " << sol->norm() << "\n";
  TEST_EQUALITY( sol->norm()<eps, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, HelmholtzOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );


  S mulI = 5.;
  S mulL = 3.;

  pl->set<S>( "alpha2", mulI );
  pl->set<S>( "Re", 1./mulL );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto x = Pimpact::create<Pimpact::VectorField>(space);
  auto bs= Pimpact::create<Pimpact::VectorField>(space);
  auto b = Pimpact::create<Pimpact::VectorField>(space);

	S n = std::sqrt( (S)x->getLength() );

  auto op= Pimpact::create<Pimpact::HelmholtzOp>( space );

  op->print();

  // test in x direction
  x->initField( Pimpact::PoiseuilleFlow2D_inX );
  bs->init( Teuchos::tuple( 8./std::pow(space->getDomainSize()->getSize(Pimpact::Y),2), 0., 0. ) );
  bs->add( mulI, *x, mulL, *bs );

	auto para = Teuchos::parameterList();
	para->set<S>( "mulI", mulI );
	op->setParameter( para );
  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);
//	b->write(); bs->write(1);

  std::cout << "error: " << bs->norm()/n << "\n";
  TEST_EQUALITY( bs->norm()/n<eps, true );

  // test in y direction
  x->initField( Pimpact::PoiseuilleFlow2D_inY );
  bs->init( Teuchos::tuple( 0., 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2), 0. ) );
  bs->add( mulI, *x, mulL, *bs );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);
//	b->write(2); bs->write(3);

  std::cout << "error: " << bs->norm()/n << "\n";
  TEST_EQUALITY( bs->norm()/n<eps, true );

  // test in z direction
  x->initField( Pimpact::PoiseuilleFlow2D_inZ );
  bs->init( Teuchos::tuple( 0., 0., 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2) ) );
  bs->add( mulI, *x, mulL, *bs );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);
//  b->write(4); bs->write(5);

  std::cout << "error: " << bs->norm()/n << "\n";
  TEST_EQUALITY( bs->norm()/n<eps, true );

  // the circle XY test
  x->initField( Pimpact::Circle2D );
  bs->initField( Pimpact::Circle2D );
  bs->scale( mulI );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b );
//	x->write(4); b->write(5); bs->write(6);

  std::cout << "error: " << bs->norm()/n << "\n";
  TEST_EQUALITY( bs->norm()/n<eps, true );

  // the circle XZ test
  x->initField( Pimpact::Circle2D_inXZ );
  bs->initField( Pimpact::Circle2D_inXZ );
  bs->scale( mulI );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b );
//	x->write(4); b->write(5); bs->write(6);

  std::cout << "error: " << bs->norm()/n << "\n";
  TEST_EQUALITY( bs->norm()/n<eps, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, DivGradO2Op ) {

	const int dNC = 2;

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto x = Pimpact::create<Pimpact::ScalarField>( space );
  auto x2= Pimpact::create<Pimpact::ScalarField>( space );
  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  // zero test
  b->initField( Pimpact::ConstField, 0. );
  x->init(2.);

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );

  op->print();

  op->apply( *b, *x );

	S n = std::sqrt( x->getLength() );
	std::cout << "error: " << x->norm()/n << "\n";
  TEST_EQUALITY( x->norm()/n<eps, true );

  // random test
  b->random();
  x->init(2.);

  op->apply(*b,*x);

  TEST_EQUALITY( x->norm()/n>eps, true );

  //b->initField( Pimpact::Grad2D_inX );
//  b->initField( Pimpact::Grad2D_inX );
  b->initField( Pimpact::ConstField, 1. );
  x->random();

  op->apply( *b, *x );

	std::cout << "error: " << x->norm()/n << "\n";
	TEST_EQUALITY( x->norm()/n<eps, true );

	// consistency test
  auto op2 = Pimpact::create<Pimpact::DivGradOp>( space );
	for( int dir=0; dir<=6; ++dir ) {
		b->initField( static_cast<Pimpact::EScalarField>(dir) );
		op->apply( *b, *x );
		op2->apply( *b, *x2 );

		x2->add( 1., *x2, -1., *x );
		S errInf = x2->norm( Belos::InfNorm );
		S err2 = x2->norm( Belos::TwoNorm )/std::sqrt( x2->getLength() );
		if( 0==space->rankST() )
			std::cout << "consistency error: inf: " << errInf << ", two: " << err2 << "\n";
		TEST_EQUALITY( errInf<eps, true );
		TEST_EQUALITY( err2<eps, true );

	}

}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, DivGradO2Smoother, SType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	Teuchos::RCP< Pimpact::ScalarField< Pimpact::Space<S,O,d,dNC> > > x =
		Pimpact::create<Pimpact::ScalarField>( space );

  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );

  auto ppl = Teuchos::parameterList();

  auto smoother = Pimpact::create<SType>( op, ppl );


	// --- zero rhs test --- 
	x->random();
	b->init( 0. );

	if( 0==space->rankST() ) {
		std::cout << "\nstep\terror\t\trate\n";
		std::cout << 0 << "\t" << 1. << "\n";
	}
	S err0 = x->norm();
	S errP = err0;

	for( int i=1; i<=20; ++i ) {
		smoother->apply( *b, *x );
		S err = x->norm();
		if( 0==space->rankST() ) {
				std::cout << i << "\t" << err/err0 << "\t" << err/errP << "\n";
		}
		errP = err;
	}


  // --- consistency test ---
	for( int dir=1; dir<=6; ++dir ) {
		x->initField( static_cast<Pimpact::EScalarField>(dir) );
		x->level();
		auto xp = x->clone( Pimpact::DeepCopy );

		op->apply( *x, *b );

		smoother->apply( *b, *x );
		x->level();

		xp->add( 1., *xp, -1., *x );
		S err2 = xp->norm()/std::sqrt( static_cast<S>(xp->getLength()) );
		S errInf = xp->norm( Belos::InfNorm );

		if( 0==space->rankST() )
			std::cout << "consistency for " << dir << ": ||" << err2 << "||_2, ||" << errInf << "||_inf\n";

		TEST_EQUALITY( err2<eps, true );
		TEST_EQUALITY( errInf<eps, true );
	}

}

using JType = Pimpact::DivGradO2JSmoother<Pimpact::DivGradO2Op<SpaceT> >;
using SORType = Pimpact::DivGradO2SORSmoother<Pimpact::DivGradO2Op<SpaceT> >;

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, JType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, SORType )



TEUCHOS_UNIT_TEST( BasicOperator, ForcingOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

	auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto force = Pimpact::create<Pimpact::VectorField>( space );
	auto res = Pimpact::create<Pimpact::VectorField>( space );

	vel->init(1.);
	force->initField( Pimpact::EVectorField(11) );
	res->random();

	auto op = Pimpact::createForcingOp( force, 1. );

	op->apply( *vel, *res );
	vel->add( 1., *res, -1., *force );

	TEST_EQUALITY( vel->norm()<eps, true );
}





template<class T> using MSFT = Pimpact::MultiField< Pimpact::ScalarField<T> >;
template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;

TEUCHOS_UNIT_TEST( MultiOperator, InverseOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto x = Pimpact::create<Pimpact::ScalarField>( space );
  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );

  auto opinv = Pimpact::create<MOP>( op );

  // Grad in x
  x->initField( Pimpact::Grad2D_inX );
//  x->write(0);

  op->apply(*x,*b);
//  b->write(1);

  x->init();

  opinv->apply( *b, *x );

//  x->write(2);

  b->initField( Pimpact::Grad2D_inX );
  x->add( -1, *x, 1., *b );
  std::cout << "error: " << x->norm()/b->norm() << "\n";
//  TEST_EQUALITY( x->norm()/b->norm()<0.5, true );

  // Grad in y
  x->initField( Pimpact::Grad2D_inY );
//  x->write(3);

  op->apply(*x,*b);
//  b->write(4);

  x->init();

  opinv->apply( *b, *x );

//  x->write(5);

  b->initField( Pimpact::Grad2D_inY );
  x->add( -1, *x, 1., *b );
  std::cout << "error: " << x->norm()/b->norm() << "\n";
//  TEST_EQUALITY( x->norm()/b->norm()<0.5, true );

}



TEUCHOS_UNIT_TEST( MultiOperator, Add2Op ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );


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


 for( O i=0; i<10; ++i ) {
	 x->getFieldPtr(i)->initField(Pimpact::RankineVortex2D );
 }

// x->getFieldPtr(0)->write();

 op->apply( *x, *y);

// y->getFieldPtr(0)->write(99);

}


template<class T> using WUP = Pimpact::MultiOpUnWrap<Pimpact::MultiOpWrap<T> >;


TEUCHOS_UNIT_TEST( MultiOperator, MulitOpUnWrap ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

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

//  x->write();

  opUW2->apply( *x, *y);

//  y->write(99);

}



TEUCHOS_UNIT_TEST( MultiModeOperator, HelmholtzOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );


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

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

	pl->set<S>("Re", 1., "Reynolds number");
	pl->set<S>("alpha2", 1.,
			"Womersley square \alpha^2");
  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  typedef Pimpact::MultiField< Pimpact::ModeField< Pimpact::VectorField<SpaceT> > > MF;
  typedef Pimpact::OperatorBase<MF> BOp;


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

  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm() );
  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm() );
  Belos::OperatorTraits<S,MF,BOp>::Apply(*op,*mv,*mv2);

}




TEUCHOS_UNIT_TEST( MultiModeOperator, DtLapOp ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

//	pl->set( "Re", 1. );
//  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );
//
//  auto velc = Pimpact::create<Pimpact::VectorField>(space);
//  auto vels = Pimpact::create<Pimpact::VectorField>(space);
//
//  auto vel = Pimpact::create<MVF>( space );
//
//  auto mv = Pimpact::createMultiField(*vel,1);
//
//  mv->random();
//
//  auto mv2 = mv->clone(1);
//  auto mv3 = mv->clone(1);
//  mv2->init(0.);
//  mv3->init(0.);
//
//  auto A = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp( space, 1., 0. ) );
//
//	mv->random();
//
//  A->apply( *mv, *mv2 );
//
//	auto diff = mv->getConstFieldPtr(0)->getConstCFieldPtr()->clone(Pimpact::ShallowCopy);
//
//	diff->add( 1., mv->getConstFieldPtr(0)->getConstCField(), -1.,  mv2->getConstFieldPtr(0)->getConstSField() );
//
//  TEST_EQUALITY( diff->norm(Belos::InfNorm)<eps, true );
//
//	diff->add( 1., mv->getConstFieldPtr(0)->getConstSField(), -1.,  mv2->getConstFieldPtr(0)->getConstCField() );
//  TEST_EQUALITY( diff->norm(Belos::InfNorm)<eps, true );
//
//  mv2->init(0.);
//
//  auto A2 = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp( space, 0., 1. ) );
//  auto A3 = Pimpact::createMultiOpWrap( Pimpact::createModeOpWrap( Pimpact::create<Pimpact::HelmholtzOp>( space ) ) );
//
//  A2->apply(*mv,*mv2);
//  A3->apply(*mv,*mv3);
//
//
//	auto diff2 = mv2->clone();
//	diff2->add( 1., *mv2, -1., *mv3 );
//
//  TEST_EQUALITY( diff2->norm(Belos::InfNorm)<eps, true );

}




TEUCHOS_UNIT_TEST( MultiModeOperator, TripleCompostion ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

	pl->set( "Re", 1. );
	pl->set( "alpha2", 1. );
  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > > MVF;
  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<SpaceT> > > MSF;

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

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

 // processor grid size
 pl->set("npx", (2==dim)?4:2 );
 pl->set("npy",            2 );
 pl->set("npz", (2==dim)?1:2 );

 auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

 using Teuchos::ParameterList;
 using Teuchos::parameterList;
 using Teuchos::RCP;
 using Teuchos::rcp; // Save some typing


 typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > > MVF;

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

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", (2==dim)?4:2 );
	pl->set("npy",            2 );
	pl->set("npz", (2==dim)?1:2 );

	//pl->set("alpha2",1.);
	//pl->set("Re",1./14.);
	auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

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

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto mv1 = Pimpact::createMultiHarmonicVectorField( space );
  auto mv2 = Pimpact::createMultiHarmonicVectorField( space );

  auto op = Pimpact::createMultiHarmonicConvectionOp( space );

  op->assignField( *mv1 );
  op->apply( *mv1, *mv2 );

}



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicOpWrap ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

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

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

	auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto mv1 = Pimpact::createMultiHarmonicVectorField( space );
	auto mv2 = Pimpact::createMultiHarmonicVectorField( space );

	auto op = Pimpact::createMultiDtHelmholtz( space );

	op->apply( *mv1, *mv2 );

}



TEUCHOS_UNIT_TEST( CompoundOperator, CompoundOpWrap ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

	auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

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
	TEST_EQUALITY( 0., fu->norm() );

	opV2V->apply( x->getConstFieldPtr(0)->getConstVField(), fu->getFieldPtr(0)->getVField() );

	TEST_INEQUALITY( 0., fu->norm() );


	// scalar to vector operator
	fu->init(0.);
	TEST_EQUALITY( 0., fu->norm() );

	opS2V->apply( x->getConstFieldPtr(0)->getConstSField(), fu->getFieldPtr(0)->getVField() );

	TEST_INEQUALITY( 0., fu->norm() );

	// vector to scalar operator
	fu->init(0.);
	TEST_EQUALITY( 0., fu->norm() );

	opV2S->apply( x->getConstFieldPtr(0)->getConstVField(), fu->getFieldPtr(0)->getSField() );

	TEST_INEQUALITY( 0., fu->norm() );

	// compound operator
	fu->init(0.);
	TEST_EQUALITY( 0., fu->norm() );
	op->apply(*x,*fu);
	TEST_INEQUALITY( 0., fu->norm() );

}


} // namespace
