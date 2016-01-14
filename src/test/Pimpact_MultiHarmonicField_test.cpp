#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Pimpact_Fields.hpp"

#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_FieldFactory.hpp"




namespace {


using S = double;
using O = int;
const int d = 4;
const int dNC = 4;

typedef Pimpact::Space<S,O,d,dNC>                SpaceT;
typedef typename Pimpact::ScalarField<SpaceT>    SF;
typedef typename Pimpact::VectorField<SpaceT>    VF;

bool testMpi = true;
double eps = 1e-6;

int dim = 3;
int domain = 0;

int npx = 1;
int npy = 1;
int npz = 1;
int npf = 8;

int nx = 9;
int ny = 9;
int nz = 9;
int nf = 7;

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

	clp.setOption( "dim", &dim, "dim" );
	clp.setOption( "domain", &domain, "domain" );

	clp.setOption( "npx", &npx, "" );
	clp.setOption( "npy", &npy, "" );
	clp.setOption( "npz", &npz, "" );
	clp.setOption( "npf", &npf, "" );

	clp.setOption( "nx", &nx, "" );
	clp.setOption( "ny", &ny, "" );
	clp.setOption( "nz", &nz, "" );
	clp.setOption( "nf", &nf, "" );

	pl->set( "spectral in time", true );

  pl->set( "nx", nx );
  pl->set( "ny", ny );
  pl->set( "nz", nz );
  pl->set( "nf", nf );

  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, constructor, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	const int m = field->getNumberVecs();

	TEST_EQUALITY( 1, m );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, constructor, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, constructor, VF )





TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, clone, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	auto field2 = field->clone();

	O n1(field->getNumberVecs());
	O n2(field2->getNumberVecs());

	TEST_EQUALITY( 1, n1 );
	TEST_EQUALITY( 1, n2 );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, clone, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, clone, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, InfNorm_and_init, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	S norm;

	// test different float values, assures that initial and norm work smoothly
	for( S i=0.; i<200.1; i+=10.5 ) {
		field->init(i/2.);
		norm = field->norm(Belos::InfNorm);
		TEST_FLOATING_EQUALITY( i/2., norm, eps );
	}

	// one test with infty-norm
	int rank = space->rankST();
	S init;
	int size;
	MPI_Comm_size( space->commST(),&size );
	for( S i = 0.; i<200.1; i+=10.5) {
		init = (size-1)*i-1.;
		init = std::abs(init);
		field->init(rank*i-1.);
		norm = field->norm(Belos::InfNorm);
		TEST_FLOATING_EQUALITY( init, norm, eps );
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, InfNorm_and_init, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, InfNorm_and_init, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, InitTwoNorm, FType  ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	int N = field->getLength();

	for( S i=0.; i< 200.1; ++i ) {
		field->init(i/2.);
    TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*N), field->norm(Belos::TwoNorm), eps );
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, InitTwoNorm, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, InitTwoNorm, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, dot, FType  ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field1 = Pimpact::createMultiHarmonic<FType>( space );

  auto field2 = field1->clone();

	S dot;

	TEST_EQUALITY( field1->getLength(), field2->getLength() )

	int N = field1->getLength();

	field1->init(0.);
	field2->init(1.);
	dot = field1->dot(*field2);
	TEST_EQUALITY( 0, dot );

	field1->init(1.);
	field2->init(1.);
	dot = field2->dot(*field1);
	TEST_EQUALITY( N, dot );

	field1->init(2.);
	field2->init(1.);
	dot = field1->dot(*field2);
	TEST_EQUALITY( 2*N, dot );

	field1->init(1.);
	field2->init(2.);
	dot = field1->dot(*field2);
	TEST_EQUALITY( 2*N, dot );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, dot, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, dot, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, scale, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	int N = field->getLength();
	S norm;

	field->init(1.);
	field->scale(2.);
	norm = field->norm(Belos::TwoNorm);
	TEST_EQUALITY( std::sqrt(4*N), norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, scale, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, scale, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, random, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	int N = field->getLength();
	S norm;

	field->init(1.);
	field->random();
	norm = field->norm(Belos::TwoNorm);
	TEST_INEQUALITY( std::sqrt(N), norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, random, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, random, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, add, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field1 = Pimpact::createMultiHarmonic<FType>( space );
  auto field2 = field1->clone();
  auto field3 = field1->clone();

	TEST_EQUALITY( field1->getLength(), field2->getLength() )
	TEST_EQUALITY( field2->getLength(), field3->getLength() )
	TEST_EQUALITY( field1->getLength(), field3->getLength() )

	S norm;
	int N = field1->getLength();

	field1->init(0.);
	field2->init(1./2.);
	field3->init(1./3.);

	field1->add( 2., *field2, 0., *field3);
	norm = field1->norm(Belos::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

	field1->init(0.);
	field2->init(1./2.);
	field3->init(1./3.);

	field1->add( 0., *field2, 3., *field3);
	norm = field1->norm(Belos::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

	field1->init(0.);
	field2->init(1.);
	field3->init(1.);

	field1->add( 0.5, *field2, 0.5, *field3);
	norm = field1->norm(Belos::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, add, SF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, add, VF ) 


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonic, exchange, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );


  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	field->init( space->getProcGrid()->getRankSlice(3) );
	field->exchange();

  TEST_FLOATING_EQUALITY( field->getConst0FieldPtr()->norm(Belos::InfNorm,false), 0., eps );

//	for( O i=0; i<space->nGlo(3); ++i ) {
//		field->getConstCFieldPtr(i)->print(  );
//		field->getConstSFieldPtr(i)->print(  );
//	}
  TEST_FLOATING_EQUALITY( field->getConstCFieldPtr(space->nGlo(3)-1)->norm(Belos::InfNorm,false), space->getProcGrid()->getNP(3)-1., eps );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonic, exchange, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonic, exchange, VF )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, write, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	field->init( 1. );
//	field->write();

	field->random();
//	field->write( 2 );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, write, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, write, VF )



} // end of namespace
