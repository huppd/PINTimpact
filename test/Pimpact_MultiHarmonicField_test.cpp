// Pimpact_SalarVectorSpace_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include "Teuchos_Range1D.hpp"
#include <Teuchos_CommHelpers.hpp>

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_FieldFactory.hpp"


#include <iostream>

namespace {

bool testMpi = true;
double errorTolSlack = 1e+1;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used." );
  clp.setOption(
      "error-tol-slack", &errorTolSlack,
      "Slack off of machine epsilon used to check test results" );
}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST( MultiHarmonicScalarField, constructor ) {
	// init impact
  init_impact(0,0);
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

	const int m = field->getNumberVecs();

	TEST_EQUALITY( 1, m );
}



TEUCHOS_UNIT_TEST( MultiHarmonicScalarFieldScalar, TwoNorm_and_init ) {
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

	const int n = field->getLength(true);

	// test different float values, assures that initial and norm work smoothly
	for( S i=0.; i< 200.1; ++i ) {
	  field->init(i/2.);
    TEST_EQUALITY( std::pow(i/2.,2)*n, field->norm(Belos::TwoNorm) );
	}
}



TEUCHOS_UNIT_TEST( MultiHarmonicScalarField, clone ) {
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

	auto field2 = field->clone();

	O n1(field->getNumberVecs());
	O n2(field2->getNumberVecs());

	TEST_EQUALITY( 1, n1 );
	TEST_EQUALITY( 1, n2 );
}



TEUCHOS_UNIT_TEST( MultiHarmonicScalarField, InfNorm_and_init ) {
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

	S norm;

	// test different float values, assures that initial and norm work smoothly
	for( S i=0.; i< 200.1; ++i ) {
		field->init(i/2.);
		norm = field->norm(Belos::InfNorm);
		TEST_EQUALITY( i/2., norm );
	}

	// one test with infty-norm
	int rank;
	S init;
	MPI_Comm_rank(fS->comm_,&rank);
	for( S i = 0.; i<200.1; ++i) {
		init = 3*i-1.;
		init = (init<0)?-init:init;
		field->init(rank*i-1.);
		norm = field->norm(Belos::InfNorm);
		TEST_EQUALITY( init, norm );
	}
}



TEUCHOS_UNIT_TEST( MultiHarmonicScalarField, TwoNorm_and_init ) {
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

	S norm;
	int N = field->getLength();

	// test different float values, assures that initial and norm work smoothly
	for( S i=0.; i< 200.1; ++i ) {
		field->init(i/2.);
		norm = field->norm(Belos::TwoNorm);
		TEST_EQUALITY( std::pow(i/2.,2)*N, norm );
	}
}



TEUCHOS_UNIT_TEST( MultiHarmonicScalarField, dot ) {
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field1 = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);
	auto field2 = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

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



TEUCHOS_UNIT_TEST( MultiHarmonicScalarField, scale ) {
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

	int N = field->getLength();
	S norm;

	field->init(1.);
	field->scale(2.);
	norm = field->norm(Belos::TwoNorm);
	TEST_EQUALITY( 4*N, norm)
}



TEUCHOS_UNIT_TEST( MultiHarmonicScalarField, random ) {
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

	int N = field->getLength();
	S norm;

	field->init(1.);
	field->random();
	norm = field->norm(Belos::TwoNorm);
	TEST_INEQUALITY( N, norm)
}



TEUCHOS_UNIT_TEST( MultiHarmonicScalarField, add ) {
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field1 = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);
	auto field2 = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);
	auto field3 = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

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
	TEST_EQUALITY( N, norm)

	field1->init(0.);
	field2->init(1./2.);
	field3->init(1./3.);

	field1->add( 0., *field2, 3., *field3);
	norm = field1->norm(Belos::TwoNorm);
	TEST_EQUALITY( N, norm)

	field1->init(0.);
	field2->init(1.);
	field3->init(1.);

	field1->add( 0.5, *field2, 0.5, *field3);
	norm = field1->norm(Belos::TwoNorm);
	TEST_EQUALITY( N, norm )
}



TEUCHOS_UNIT_TEST( MultiHarmonicScalarField, write ) {
  typedef double S;
  typedef int O;

	auto fS = Pimpact::createFieldSpace<O>();

	auto field = Pimpact::createMultiHarmonicScalarField<S,O>(fS,10);

	field->init( 1. );
	field->write();

	field->random();
	field->write( 2 );

}


} // end of namespace
