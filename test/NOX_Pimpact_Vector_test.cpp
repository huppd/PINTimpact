#include <iostream>
#include <vector>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "BelosTypes.hpp"

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
//#include "Pimpact_ModeField.hpp"
#include "Pimpact_CompoundField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "NOX_Pimpact_Vector.hpp"


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



TEUCHOS_UNIT_TEST( NOXPimpactVector, create_init_print ) {

	int rank = init_impact(0,0);

	typedef double S;
	typedef int O;
	typedef Pimpact::ScalarField<S,O> SF;
	typedef Pimpact::VectorField<S,O> VF;
	typedef Pimpact::ModeField<SF> MSF;
	typedef Pimpact::ModeField<VF> MVF;
	typedef Pimpact::MultiField< MSF> BMSF;
	typedef Pimpact::MultiField< MVF> BMVF;
	typedef Pimpact::CompoundField< BMVF, BMSF > CF;
	typedef NOX::Pimpact::Vector<CF> NV;


	auto fS  = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

	auto xs = Pimpact::createInitMSF<S,O>( fS );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> nx = Teuchos::rcp(new NV(x) );

	nx->init( rank );
}


TEUCHOS_UNIT_TEST( NOXPimpactVector, InfNorm_and_init ) {

	typedef double S;
	typedef int O;
	typedef Pimpact::ScalarField<S,O> SF;
	typedef Pimpact::VectorField<S,O> VF;
	typedef Pimpact::ModeField<SF> MSF;
	typedef Pimpact::ModeField<VF> MVF;
	typedef Pimpact::MultiField< MSF> BMSF;
	typedef Pimpact::MultiField< MVF> BMVF;
	typedef Pimpact::CompoundField< BMVF, BMSF > CF;
	typedef NOX::Pimpact::Vector<CF> NV;


	auto fS  = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

	auto xs = Pimpact::createInitMSF<S,O>( fS );

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
	MPI_Comm_rank(fS->comm_,&rank);
	for( double i = 0.; i<200.1; ++i) {
		init = 3*i-1.;
		init = (init<0)?-init:init;
		um->init(rank*i-1.);
		norm = um->norm(NOX::Abstract::Vector::MaxNorm);
		TEST_EQUALITY( init, norm );
	}
}



TEUCHOS_UNIT_TEST( NOXPimpactVector, TwoNorm_and_init ) {

	typedef double S;
	typedef int O;
	typedef Pimpact::ScalarField<S,O> SF;
	typedef Pimpact::VectorField<S,O> VF;
	typedef Pimpact::ModeField<SF> MSF;
	typedef Pimpact::ModeField<VF> MVF;
	typedef Pimpact::MultiField< MSF> BMSF;
	typedef Pimpact::MultiField< MVF> BMVF;
	typedef Pimpact::CompoundField< BMVF, BMSF > CF;
	typedef NOX::Pimpact::Vector<CF> NV;


	auto fS  = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

	auto xs = Pimpact::createInitMSF<S,O>( fS );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> q = Teuchos::rcp(new NV(x) );

	double norm;
	int N = q->length();

	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		q->init(i/2.);
		norm = q->norm(NOX::Abstract::Vector::TwoNorm);
		TEST_EQUALITY( std::pow(i/2.,2)*N, norm );
	}
}


TEUCHOS_UNIT_TEST( NOXPimpactVector, add ) {
  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField< MSF> BMSF;
  typedef Pimpact::MultiField< MVF> BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
  typedef NOX::Pimpact::Vector<CF> NV;


  auto fS  = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

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
  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField< MSF> BMSF;
  typedef Pimpact::MultiField< MVF> BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
  typedef NOX::Pimpact::Vector<CF> NV;


  auto fS  = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

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
  typedef double                                 S;
  typedef int                                    O;
  typedef Pimpact::ScalarField<S,O>             SF;
  typedef Pimpact::VectorField<S,O>             VF;
  typedef Pimpact::ModeField<SF>               MSF;
  typedef Pimpact::ModeField<VF>               MVF;
  typedef Pimpact::MultiField< MSF>           BMSF;
  typedef Pimpact::MultiField< MVF>           BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF >  CF;
  typedef NOX::Pimpact::Vector<CF>              NV;


  auto fS  = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

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
  TEST_EQUALITY( N, norm );

  vel1->init(-1.);
  vel2->reciprocal( *vel1 );
  norm = vel2->norm(*vel1);
  TEST_EQUALITY( N, norm );

  vel1->init(2.);
  vel2->reciprocal( *vel1 );
  norm = vel2->norm(*vel1);
  TEST_EQUALITY( N, norm );

  vel1->random();
  vel2->reciprocal( *vel1 );
  norm = vel2->norm(*vel1);
  TEST_EQUALITY( N, norm );
}


TEUCHOS_UNIT_TEST( NOXPimpactVector, scale2 ) {
  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField< MSF> BMSF;
  typedef Pimpact::MultiField< MVF> BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
  typedef NOX::Pimpact::Vector<CF> NV;


  auto fS  = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

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
	typedef double S;
	typedef int O;
	typedef Pimpact::ScalarField<S,O> SF;
	typedef Pimpact::VectorField<S,O> VF;
	typedef Pimpact::ModeField<SF> MSF;
	typedef Pimpact::ModeField<VF> MVF;
	typedef Pimpact::MultiField< MSF> BMSF;
	typedef Pimpact::MultiField< MVF> BMVF;
	typedef Pimpact::CompoundField< BMVF, BMSF > CF;
	typedef NOX::Pimpact::Vector<CF> NV;


	auto fS  = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

	auto xs = Pimpact::createInitMSF<S,O>( fS );

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
	typedef double S;
	typedef int O;
	typedef Pimpact::ScalarField<S,O> SF;
	typedef Pimpact::VectorField<S,O> VF;
	typedef Pimpact::ModeField<SF> MSF;
	typedef Pimpact::ModeField<VF> MVF;
	typedef Pimpact::MultiField< MSF> BMSF;
	typedef Pimpact::MultiField< MVF> BMVF;
	typedef Pimpact::CompoundField< BMVF, BMSF > CF;
	typedef NOX::Pimpact::Vector<CF> NV;


	auto fS  = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

	auto xs = Pimpact::createInitMSF<S,O>( fS );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> q = Teuchos::rcp(new NV(x) );

	int N = q->length();
	double norm;

	q->init(1.);
	q->scale(2.);
	norm = q->norm(NOX::Abstract::Vector::TwoNorm);
	TEST_EQUALITY( 4*N, norm)
}


TEUCHOS_UNIT_TEST( NOXPimpactVector, random ) {
	typedef double S;
	typedef int O;
	typedef Pimpact::ScalarField<S,O> SF;
	typedef Pimpact::VectorField<S,O> VF;
	typedef Pimpact::ModeField<SF> MSF;
	typedef Pimpact::ModeField<VF> MVF;
	typedef Pimpact::MultiField< MSF> BMSF;
	typedef Pimpact::MultiField< MVF> BMVF;
	typedef Pimpact::CompoundField< BMVF, BMSF > CF;
	typedef NOX::Pimpact::Vector<CF> NV;


	auto fS  = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

	auto xs = Pimpact::createInitMSF<S,O>( fS );

	auto x  = Pimpact::createCompoundField( xv, xs );

	Teuchos::RCP<NV> q = Teuchos::rcp(new NV(x) );

	int N = q->length();
	double norm;

	q->init(1.);
	q->random();
	norm = q->norm(NOX::Abstract::Vector::TwoNorm);
	TEST_INEQUALITY( N, norm)
}


TEUCHOS_UNIT_TEST( NOXPimpactVector, update ) {
	typedef double S;
	typedef int O;
	typedef Pimpact::ScalarField<S,O> SF;
	typedef Pimpact::VectorField<S,O> VF;
	typedef Pimpact::ModeField<SF> MSF;
	typedef Pimpact::ModeField<VF> MVF;
	typedef Pimpact::MultiField< MSF> BMSF;
	typedef Pimpact::MultiField< MVF> BMVF;
	typedef Pimpact::CompoundField< BMVF, BMSF > CF;
	typedef NOX::Pimpact::Vector<CF> NV;


	auto fS  = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto xv = Pimpact::createInitMVF<S,O>(Pimpact::ZeroFLow, fS, iIS, fIS );

	auto xs = Pimpact::createInitMSF<S,O>( fS );

	auto x  = Pimpact::createCompoundField( xv, xs );


	Teuchos::RCP<NV> vel1 = Teuchos::rcp(new NV(x) );
	Teuchos::RCP<NV> vel2 = Teuchos::rcp(new NV(x->clone()) );
	Teuchos::RCP<NV> vel3 = Teuchos::rcp_static_cast<NV>(vel2->clone());

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
	TEST_EQUALITY( N, norm)

	vel1->init(0.);
	vel2->init(1./2.);
	vel3->init(1./3.);

	vel1->update( 0., *vel2, 3., *vel3);
	norm = vel1->norm(NOX::Abstract::Vector::TwoNorm);
	TEST_EQUALITY( N, norm)

	vel1->init(0.);
	vel2->init(1.);
	vel3->init(1.);

	vel1->update( 0.5, *vel2, 0.5, *vel3);
	norm = vel1->norm(NOX::Abstract::Vector::TwoNorm);
	TEST_EQUALITY( N, norm )
}


} // namespace

