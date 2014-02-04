// Pimpact_SalarVector_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "BelosTypes.hpp"

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_ScalarField.hpp"

#include <iostream>
#include <cmath>


namespace {

	bool testMpi = true;
	double errorTolSlack = 1e+1;

	TEUCHOS_STATIC_SETUP()
	{
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

//  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
//  {
//    if (testMpi) {
//      return Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
//    }
//    return rcp(new Teuchos::SerialComm<int>());
//  }


	TEUCHOS_UNIT_TEST( ScalarField, create_init_print ) {
    int rank = init_impact(0,0);
		// init impact
		auto sVS = Pimpact::createFieldSpace<int>();
		auto p = Pimpact::createScalarField<double,int>(sVS);

//		sVS->print();
		p->init(rank);
//		p->print();

		TEST_EQUALITY( 0, 0 );
	}


TEUCHOS_UNIT_TEST( ScalarField, InfNorm_and_init ) {
	auto sVS = Pimpact::createFieldSpace<int>();
	auto p = Pimpact::createScalarField<double,int>(sVS);

	double norm;

	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		p->init(i/2.);
		norm = p->norm(Belos::InfNorm);
		TEST_EQUALITY( i/2., norm );
	}

	// one test with infty-norm
	int rank;
	double init;
	MPI_Comm_rank(sVS->comm_,&rank);
	for( double i = 0.; i<200.1; ++i) {
		init = 3*i-1.;
		init = (init<0)?-init:init;
		p->init(rank*i-1.);
		norm = p->norm(Belos::InfNorm);
		TEST_EQUALITY( init, norm );
	}
}


	TEUCHOS_UNIT_TEST( ScalarField, TwoNorm_and_init ) {

		auto sVS = Pimpact::createFieldSpace<int>();
		auto p = Pimpact::createScalarField<double,int>(sVS);


	  // test different float values, assures that initial and norm work smoothly
		for( double i=0.; i< 200.1; ++i ) {
			p->init(i/2.);
			TEST_EQUALITY( std::pow(i/2.,2)*p->getVecLength(), p->norm(Belos::TwoNorm) );
		}
	}


TEUCHOS_UNIT_TEST( ScalarField, dot ) {

		auto fS = Pimpact::createFieldSpace<int>();
		auto p = Pimpact::createScalarField<double,int>(fS);
		auto q = Pimpact::createScalarField<double,int>(fS);


		int Np = p->getVecLength();
		int Nq = q->getVecLength();
		double dot;

		TEST_EQUALITY( Np , Nq );

//		p->init(0.);
//		q->init(1.);
//		dot = p->dot(*q);
//		TEST_EQUALITY( 0, dot );

//		p->init(1.);
//		q->init(1.);
//		dot = p->dot(*q);
//		TEST_EQUALITY( Np, dot );

//		p->init(2.);
//		q->init(1.);
//		dot = p->dot(*q);
//		TEST_EQUALITY( 2*Nq, dot );
//
		p->init(1.);
		q->init(2.);
		dot = p->dot(*q);
		TEST_EQUALITY( 2*Np, dot );

	}


TEUCHOS_UNIT_TEST( ScalarField, scale ) {

		auto sVS = Pimpact::createFieldSpace<int>();
		auto p = Pimpact::createScalarField<double,int>(sVS);

		double norm;
		int N = p->getVecLength();

		p->init(1.);
		p->scale(2.);
		norm = p->norm(Belos::TwoNorm);
		TEST_EQUALITY( 4*N, norm)

	}


TEUCHOS_UNIT_TEST( ScalarField, random ) {

		auto sVS = Pimpact::createFieldSpace<int>();
		auto p = Pimpact::createScalarField<double,int>(sVS);

		double norm;
		int N = p->getVecLength();

		p->init(1.);
		p->random();
		norm = p->norm(Belos::TwoNorm);
		TEST_INEQUALITY( N, norm)

	}


TEUCHOS_UNIT_TEST( ScalarField, add ) {

		auto sVS = Pimpact::createFieldSpace<int>();
//		auto p = Pimpact::createScalarField<double,int>(sVS);
		auto q = Pimpact::createScalarField<double,int>(sVS);
//		auto r = Pimpact::createScalarField<double,int>(sVS);
		auto r(q);
		auto p(q);

		double norm;
		int N = p->getVecLength();


//		p->init(0.);
//		q->init(1./2.);
//		r->init(1./3.);
//
//		p->add( 2., *q, 0, *r);
//		norm = p->norm(Belos::TwoNorm);
//		TEST_EQUALITY( N, norm)

		q->init(1.);
		r->init(1./3.);

		p->add( 0., *q, 3., *r);
		norm = p->norm(Belos::TwoNorm);
		TEST_EQUALITY( N, norm)
	}

TEUCHOS_UNIT_TEST( ScalarField, write ) {

		auto sVS = Pimpact::createFieldSpace<int>();
		auto p = Pimpact::createScalarField<double,int>(sVS);

		p->init(1.);
		p->write();

		p->random();
		p->write(1);

		TEST_EQUALITY( 0, 0)
	}
} // namespace

