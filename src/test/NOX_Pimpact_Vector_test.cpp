#include <iostream>
#include <vector>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Fields.hpp"

#include "NOX_Pimpact_Vector.hpp"



namespace {

typedef double S;
typedef int O;

typedef typename Pimpact::Space<S,O,3,4> SpaceT;

typedef Pimpact::ScalarField<SpaceT> SF;
typedef Pimpact::VectorField<SpaceT> VF;
typedef Pimpact::ModeField<SF> MSF;
typedef Pimpact::ModeField<VF> MVF;
typedef Pimpact::MultiField< MSF> BMSF;
typedef Pimpact::MultiField< MVF> BMVF;
typedef Pimpact::CompoundField< BMVF, BMSF > CF;
typedef NOX::Pimpact ::Vector<CF> NV;

bool testMpi = true;
double eps = 1e+1;

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

	pl->set( "dim", 3 );
	pl->set( "domain", 0 );

	pl->set( "nx", 25 );
	pl->set( "ny", 17 );
	pl->set( "nz",  9 );

}



TEUCHOS_UNIT_TEST( NOXPimpactVector, create_init_print ) {

	auto space = Pimpact::createSpace( pl );

	auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

	auto xs = Pimpact::createInitMSF( space );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> nx = Teuchos::rcp(new NV(x) );

	nx->init( space->rankST() );
}


TEUCHOS_UNIT_TEST( NOXPimpactVector, InfNorm_and_init ) {

	auto space = Pimpact::createSpace( pl );

	auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

	auto xs = Pimpact::createInitMSF( space );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> um = Teuchos::rcp(new NV(x) );

	double norm;
	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		um->init(i/2.);
		norm = um->norm(NOX::Abstract::Vector::MaxNorm);
		TEST_EQUALITY( i/2., norm );
	}

	// one test with infty-norm
	int rank;
	double init;
	MPI_Comm_rank( space->comm(), &rank );
	for( double i = 0.; i<200.1; ++i) {
		init = 3*i-1.;
		init = (init<0)?-init:init;
		um->init(rank*i-1.);
		norm = um->norm(NOX::Abstract::Vector::MaxNorm);
		TEST_EQUALITY( init, norm );
	}
}



TEUCHOS_UNIT_TEST( NOXPimpactVector, TwoNorm_and_init ) {

	auto space = Pimpact::createSpace( pl );

	auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

	auto xs = Pimpact::createInitMSF( space );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> q = Teuchos::rcp(new NV(x) );

	double norm;
	int N = q->length();

	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		q->init(i/2.);
		norm = q->norm(NOX::Abstract::Vector::TwoNorm);
    TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*N), norm, eps );
	}

}



TEUCHOS_UNIT_TEST( NOXPimpactVector, add ) {

	auto space = Pimpact::createSpace( pl );

  auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

  auto xs = Pimpact::createInitMSF( space );

  auto x  = Pimpact::createCompoundField( xv, xs );

  Teuchos::RCP<NV> vel1 = Teuchos::rcp(new NV(x) );
  Teuchos::RCP<NV> vel2 = Teuchos::rcp(new NV( x->clone(Pimpact::ShallowCopy) ) );

  double dot;

  TEST_EQUALITY( vel1->length(), vel2->length() )

  int N = vel1->length();

  vel1->init(0.);
  vel2->abs( *vel1 );
  dot = vel1->innerProduct(*vel2);
  TEST_EQUALITY( 0, dot );

  vel1->init(1.);
  vel2->abs( *vel1 );
  dot = vel1->innerProduct(*vel2);
  TEST_EQUALITY( N, dot );

  vel1->init(-1.);
  vel2->abs( *vel1 );
  dot = vel2->innerProduct(*vel1);
  TEST_EQUALITY( -N, dot );

}


TEUCHOS_UNIT_TEST( NOXPimpactVector, reciprocal ) {

	auto space = Pimpact::createSpace( pl );

  auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

  auto xs = Pimpact::createInitMSF( space );

  auto x  = Pimpact::createCompoundField( xv, xs );

  Teuchos::RCP<NV> vel1 = Teuchos::rcp(new NV(x) );
  Teuchos::RCP<NV> vel2 = Teuchos::rcp(new NV( x->clone(Pimpact::ShallowCopy) ) );

  double dot;

  TEST_EQUALITY( vel1->length(), vel2->length() )

  int N = vel1->length();

  vel1->init(0.);
  vel2->reciprocal( *vel1 );
  dot = vel1->innerProduct(*vel2);
  TEST_EQUALITY( 0., dot );

  vel1->init(1.);
  vel2->reciprocal( *vel1 );
  dot = vel1->innerProduct(*vel2);
  TEST_EQUALITY( N, dot );

  vel1->init(-1.);
  vel2->reciprocal( *vel1 );
  dot = vel2->innerProduct(*vel1);
  TEST_EQUALITY( N, dot );

  vel1->random();
  vel2->reciprocal( *vel1 );
  dot = vel2->innerProduct(*vel1);
  TEST_EQUALITY( N, dot );

}


TEUCHOS_UNIT_TEST( NOXPimpactVector, norm_weighted ) {

	auto space = Pimpact::createSpace( pl );

  auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

  auto xs = Pimpact::createInitMSF( space );

  auto x  = Pimpact::createCompoundField( xv, xs );

  Teuchos::RCP<NV> vel1 = Teuchos::rcp(new NV(x) );
  Teuchos::RCP<NV> vel2 = Teuchos::rcp(new NV( x->clone(Pimpact::ShallowCopy) ) );

  double norm;

  TEST_EQUALITY( vel1->length(), vel2->length() )

  int N = vel1->length();

  vel1->init(0.);
  vel2->reciprocal( *vel1 );
  norm = vel1->norm(*vel2);
  TEST_EQUALITY( 0., norm );

  vel1->init(1.);
  vel2->reciprocal( *vel1 );
  norm = vel1->norm(*vel2);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );


  vel1->init(-1.);
  vel2->reciprocal( *vel1 );
  norm = vel2->norm(*vel1);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

  vel1->init(2.);
  vel2->reciprocal( *vel1 );
  norm = vel2->norm(*vel1);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

  vel1->random();
  vel2->reciprocal( *vel1 );
  norm = vel2->norm(*vel1);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

}



TEUCHOS_UNIT_TEST( NOXPimpactVector, scale2 ) {

	auto space = Pimpact::createSpace( pl );

  auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

  auto xs = Pimpact::createInitMSF( space );

  auto x  = Pimpact::createCompoundField( xv, xs );

  Teuchos::RCP<NV> vel1 = Teuchos::rcp(new NV(x) );
  Teuchos::RCP<NV> vel2 = Teuchos::rcp(new NV( x->clone(Pimpact::ShallowCopy) ) );

  double dot;

  TEST_EQUALITY( vel1->length(), vel2->length() )

  int N = vel1->length();

  vel1->init(0.);
  vel2->init(1.);
  vel2->scale( *vel1 );
  dot = vel1->innerProduct(*vel2);
  TEST_EQUALITY( 0., dot );

  vel1->init(1.);
  vel2->init(1.);
  vel2->scale( *vel1 );
  dot = vel1->innerProduct(*vel2);
  TEST_EQUALITY( N, dot );

  vel1->init(-1.);
  vel2->init(1.);
  vel2->scale( *vel1 );
  dot = vel2->innerProduct(*vel1);
  TEST_EQUALITY( N, dot );

  vel1->init(2.);
  vel2->init(1.);
  vel2->scale( *vel1 );
  dot = vel2->innerProduct(*vel1);
  TEST_EQUALITY( 4*N, dot );

}


TEUCHOS_UNIT_TEST( NOXPimpactVector, innerProduct ) {

	auto space = Pimpact::createSpace( pl );

	auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

	auto xs = Pimpact::createInitMSF( space );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> vel1 = Teuchos::rcp(new NV(x) );
	Teuchos::RCP<NV> vel2 = Teuchos::rcp(new NV( x->clone() ) );

	double dot;

	TEST_EQUALITY( vel1->length(), vel2->length() )

	int N = vel1->length();

	vel1->init(0.);
	vel2->init(1.);
	dot = vel1->innerProduct(*vel2);
	TEST_EQUALITY( 0, dot );

	vel1->init(1.);
	vel2->init(1.);
	dot = vel2->innerProduct(*vel1);
	TEST_EQUALITY( N, dot );

	vel1->init(2.);
	vel2->init(1.);
	dot = vel1->innerProduct(*vel2);
	TEST_EQUALITY( 2*N, dot );

	vel1->init(1.);
	vel2->init(2.);
	dot = vel1->innerProduct(*vel2);
	TEST_EQUALITY( 2*N, dot );
}


TEUCHOS_UNIT_TEST( NOXPimpactVector, scale ) {

	auto space = Pimpact::createSpace( pl );

	auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

	auto xs = Pimpact::createInitMSF( space );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> q = Teuchos::rcp(new NV(x) );

	int N = q->length();
	double norm;

	q->init(1.);
	q->scale(2.);
	norm = q->norm(NOX::Abstract::Vector::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(4*N), norm, eps );

}


TEUCHOS_UNIT_TEST( NOXPimpactVector, random ) {

	auto space = Pimpact::createSpace( pl );

	auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

	auto xs = Pimpact::createInitMSF( space );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> q = Teuchos::rcp(new NV(x) );

	int N = q->length();
	double norm;

	q->init(1.);
	q->random();
	norm = q->norm(NOX::Abstract::Vector::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

}



TEUCHOS_UNIT_TEST( NOXPimpactVector, update ) {

	auto space = Pimpact::createSpace( pl );

	auto xv = Pimpact::createInitMVF(Pimpact::Zero2DFlow, space );

	auto xs = Pimpact::createInitMSF( space );

	auto x  = Pimpact::createCompoundField( xv, xs );


	Teuchos::RCP<NV> vel1 = Teuchos::rcp(new NV(x) );
	Teuchos::RCP<NV> vel2 = Teuchos::rcp(new NV(x->clone()) );
	Teuchos::RCP<NV> vel3 = Teuchos::rcp_dynamic_cast<NV>(vel2->clone());

	TEST_EQUALITY( vel1->length(), vel2->length() )
	TEST_EQUALITY( vel2->length(), vel3->length() )
	TEST_EQUALITY( vel1->length(), vel3->length() )

	double norm;
	int N = vel1->length();

	vel1->init(0.);
	vel2->init(1./2.);
	vel3->init(1./3.);

	vel1->update( 2., *vel2, 0., *vel3);
	norm = vel1->norm(NOX::Abstract::Vector::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

	vel1->init(0.);
	vel2->init(1./2.);
	vel3->init(1./3.);

	vel1->update( 0., *vel2, 3., *vel3);
	norm = vel1->norm(NOX::Abstract::Vector::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

	vel1->init(0.);
	vel2->init(1.);
	vel3->init(1.);

	vel1->update( 0.5, *vel2, 0.5, *vel3);
	norm = vel1->norm(NOX::Abstract::Vector::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(N), norm, eps );

}


} // end of namespace
