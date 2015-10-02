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


typedef double S;
typedef int O;
const int d = 4;
const int dNC = 4;

typedef Pimpact::Space<S,O,d,dNC>                SpaceT;
typedef typename Pimpact::ScalarField<SpaceT>    SF;
typedef typename Pimpact::VectorField<SpaceT>    VF;

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
      "dim", &dim,
      "dim" );
  clp.setOption(
      "domain", &domain,
      "domain" );

  pl->set( "nx", 25 );
  pl->set( "ny", 17 );
  pl->set( "nz",  9 );
  pl->set( "nf", 10 );

  pl->set( "npx", 2 );
  pl->set( "npy", 2 );
  pl->set( "npz", 2 );
  pl->set( "npf", 1 );

}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, constructor, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	const int m = field->getNumberVecs();

	TEST_EQUALITY( 1, m );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, constructor, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, constructor, VF )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, TwoNorm_and_init, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	const int n = field->getLength(true);

	// test different float values, assures that initial and norm work smoothly
	for( S i=0.; i< 200.1; ++i ) {
	  field->init(i/2.);
    TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*n), field->norm(Belos::TwoNorm), eps );
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, TwoNorm_and_init, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, TwoNorm_and_init, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicField, clone, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

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

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	S norm;

	// test different float values, assures that initial and norm work smoothly
	for( S i=0.; i< 200.1; ++i ) {
		field->init(i/2.);
		norm = field->norm(Belos::InfNorm);
		TEST_EQUALITY( i/2., norm );
	}

	// one test with infty-norm
	int rank = space->rankST();
	S init;
	int size;
	MPI_Comm_size( space->comm(),&size );
	for( S i = 0.; i<200.1; ++i) {
		init = (size-1)*i-1.;
		init = std::abs(init);
		field->init(rank*i-1.);
		norm = field->norm(Belos::InfNorm);
		TEST_EQUALITY( init, norm );
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, InfNorm_and_init, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicField, InfNorm_and_init, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicScalarField, TwoNorm_and_init, FType  ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	int N = field->getLength();

	for( S i=0.; i< 200.1; ++i ) {
		field->init(i/2.);
    TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*N), field->norm(Belos::TwoNorm), eps );
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, TwoNorm_and_init, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, TwoNorm_and_init, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicScalarField, dot, FType  ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

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

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, dot, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, dot, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicScalarField, scale, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	int N = field->getLength();
	S norm;

	field->init(1.);
	field->scale(2.);
	norm = field->norm(Belos::TwoNorm);
	TEST_EQUALITY( std::sqrt(4*N), norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, scale, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, scale, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicScalarField, random, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	int N = field->getLength();
	S norm;

	field->init(1.);
	field->random();
	norm = field->norm(Belos::TwoNorm);
	TEST_INEQUALITY( std::sqrt(N), norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, random, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, random, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicScalarField, add, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

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

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, add, SF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, add, VF ) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicScalarField, write, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto field = Pimpact::createMultiHarmonic<FType>( space );

	field->init( 1. );
	field->write();

	field->random();
	field->write( 2 );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, write, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicScalarField, write, VF )



} // end of namespace
