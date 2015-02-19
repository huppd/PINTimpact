#include <iostream>

#include "mpi.h"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosMVOPTester.hpp"

#include "BelosPimpactAdapter.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

#include "Pimpact_TimeField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_TimeOpWrap.hpp"
#include "Pimpact_DtTimeOp.hpp"
#include "Pimpact_TimeNonlinearJacobianOp.hpp"




namespace {


typedef double S;
typedef int O;
const int d = 4;
const int dNC = 4;

bool testMpi = true;
double eps = 3e-1;
int domain = 1;
int dim = 2;
//bool isImpactInit=false;


typedef Pimpact::Space<S,O,d,dNC> SpaceT;

typedef Pimpact::ScalarField<SpaceT> SF;
typedef Pimpact::VectorField<SpaceT> VF;

typedef Pimpact::TimeField<SF> TSF;
typedef Pimpact::TimeField<VF> TVF;

typedef Pimpact::MultiField<TSF> MTSF;
typedef Pimpact::MultiField<TVF> MTVF;

typedef Pimpact::OperatorBase<MTSF> SOpBase;
typedef Pimpact::OperatorBase<MTVF> VOpBase;

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

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );
  pl->set( "domain", domain );

  pl->set( "lx", 1. );
  pl->set( "ly", 1. );
  pl->set( "lz", 1. );

  pl->set( "dim", dim );

  pl->set("nx", 33 );
  pl->set("ny", 17 );
  pl->set("nz", 2 );

  pl->set("nf", 8 );

  // processor grid size
  pl->set("npx", 2 );
  pl->set("npy", 2 );
  pl->set("npz", 1 );
  pl->set("npf", 2 );

}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, create_init_print, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, false );

  space->print();

  auto p = Pimpact::create<FType>(space);

  p->init( space->rankST() );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, InfNorm_and_init, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, false );

  auto p = Pimpact::create<FType>(space);

  double norm;

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    p->init(i/2.);
    norm = p->norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( i/2., norm, eps );

  }

  // one test with infty-norm
  int rank;
  int size;
  double init;
  MPI_Comm_rank(space->comm(),&rank);
  MPI_Comm_size(space->comm(),&size);
  for( double i = 0.; i<200.1; ++i) {
    init = (size-1)*i-1.;
    init = (init<0)?-init:init;
    p->init(rank*i-1.);
    norm = p->norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( init, norm, eps );

  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, TSF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, TVF ) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, OneNorm_and_init, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, false );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    p->init(i/2.);
//    TEST_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm) );
    TEST_FLOATING_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm), eps );

  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, TSF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, TVF ) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, TwoNorm_and_init, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, false );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 200.1; ++i ) {
    p->init(i/2.);
    TEST_FLOATING_EQUALITY( std::sqrt( std::pow(i/2.,2)*p->getLength() ), p->norm(Belos::TwoNorm), eps );
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, dot, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, false );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);

  int Np = vel1->getLength();
  int Nq = vel2->getLength();
  double dot;

  TEST_EQUALITY( Np , Nq );
  int N = Np;

  vel1->init(0.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_EQUALITY( dot<eps, true );

  vel1->init(1.);
  vel2->init(1.);
  dot = vel2->dot(*vel1);
  TEST_FLOATING_EQUALITY( (S)N, dot, eps );

  vel1->init(2.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY( 2.*N, dot, eps );

  vel1->init(1.);
  vel2->init(2.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY( 2.*N, dot, eps );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, scale, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, false );

  auto p = Pimpact::create<FType>(space);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->scale(2.);
  norm = p->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(4*N), norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, random, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, false );

  auto p = Pimpact::create<FType>(space);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->random();
  norm = p->norm(Belos::TwoNorm);
  TEST_INEQUALITY( N, norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TemplateField, add, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, false );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);
  auto vel3 = Pimpact::create<FType>(space);

  TEST_EQUALITY( vel1->getLength(), vel2->getLength() )
  TEST_EQUALITY( vel2->getLength(), vel3->getLength() )
  TEST_EQUALITY( vel1->getLength(), vel3->getLength() )

  double norm;
  int N = vel1->getLength();

  vel1->init(0.);
  vel2->init(1./2.);
  vel3->init(1./3.);

  vel1->add( 2., *vel2, 0., *vel3);
  norm = vel1->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm )

  vel1->init(0.);
  vel2->init(1./2.);
  vel3->init(1./3.);

  vel1->add( 0., *vel2, 3., *vel3);
  norm = vel1->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm )

  vel1->init(0.);
  vel2->init(1.);
  vel3->init(1.);

  vel1->add( 0.5, *vel2, 0.5, *vel3);
  norm = vel1->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm )

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, write, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, false );

  auto p = Pimpact::create<FType>( space );

  p->init(1.);
  p->write();

  p->random();
  p->write(1);

  TEST_EQUALITY( 0, 0)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, TVF )



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TimeField, all, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	space->print();
  //---------------------------------------------------
  auto field1 = Pimpact::createTimeField<FType>( space );

  auto field2 = field1->clone();

  std::cout << "field1: length: " << field1->getLength() << "\n";
  std::cout << "field2: length: " << field2->getLength() << "\n";

  for( S i=0.; i< 200.1; ++i ) {
    field1->init(i/2.);
    S norm_ = field1->norm(Belos::InfNorm);
    TEST_EQUALITY( i/2., norm_ );
    field2->init(i/2.);
    norm_ = field2->norm(Belos::InfNorm);
    TEST_EQUALITY( i/2., norm_ );
  }
  field2->init(1.);
  for( double i=0.; i< 200.1; ++i ) {
    field1->init(i/2.);
    TEST_EQUALITY( (i/2.)*field1->getLength(), field1->norm(Belos::OneNorm) );
    TEST_EQUALITY( (i/2.)*field1->getLength(), field1->dot(*field2) );
  }

	// some test
	auto fields = Pimpact::createTimeField<FType>( space );

	auto msca = Pimpact::createMultiField(*fields,5);

	auto field2s = fields->clone();
	field2s->random();
	MPI_Barrier( MPI_COMM_WORLD );
	field2s->norm();
	field2s->init( space->procCoordinate()[3] );
	field2s->exchange();
	field2s->write();



//	// op test
// auto dt = Pimpact::create<Pimpact::DtTimeOp>( space );
//
////  field1->random();
// Pimpact::initVectorTimeField( field1, Pimpact::OscilatingDisc2D );
// field2->init(0);
//
// dt->apply( *field1, *field2 );
//
// field2->write( 50 );


}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TimeField, all, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TimeField, all, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TimeField, BelosMVTest, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto vel = Pimpact::create<FType>( space );

  auto mvec = Pimpact::createMultiField( *vel, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM =
		Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<S,Pimpact::MultiField<FType> >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TimeField, BelosMVTest, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TimeField, BelosMVTest, TVF )



TEUCHOS_UNIT_TEST( TimeOpearotr, TimeOpWrap ) {

	auto space = Pimpact::createSpace<S,O,d,dNC>( pl );


	auto field = Pimpact::create<TVF>( space );
	auto field1 = field->clone();
	auto field2 = field->clone();

	auto mv = Pimpact::createMultiField( *field, 10 );

	// op test
	auto op = Pimpact::createTimeOpWrap<Pimpact::HelmholtzOp<SpaceT>,true>(
			Pimpact::create<Pimpact::HelmholtzOp>(space ) );

	field2->random();
	Pimpact::initVectorTimeField( field1, Pimpact::Poiseuille_inX );
	op->apply( *field1, *field2 );
	field2->write();
	Pimpact::initVectorTimeField( field2, Pimpact::Zero2DFlow );
	Pimpact::initVectorTimeField( field2, Pimpact::Poiseuille_inX );
	field2->write(10);
	Pimpact::initVectorTimeField( field2, Pimpact::Poiseuille_inY );
	field2->write(20);
	Pimpact::initVectorTimeField( field2, Pimpact::Streaming2DFlow );
	field2->write(30);
	Pimpact::initVectorTimeField( field2, Pimpact::OscilatingDisc2D );
	field2->write(40);

	auto bop =
		Pimpact::createMultiOperatorBase( op );


	Teuchos::RCP<Belos::OutputManager<S> > myOM = Teuchos::rcp(
			new Belos::OutputManager<S>(Belos::Warnings+Belos::TimingDetails,rcp(&out,false)) );

	bool res = Belos::TestOperatorTraits< S, MTVF, VOpBase >( myOM, mv, bop );

  TEST_EQUALITY( res, true );
}

} // end of namespace
