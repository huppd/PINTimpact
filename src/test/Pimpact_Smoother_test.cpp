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

#include "Pimpact_Fields.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_TimeOpWrap.hpp"
#include "Pimpact_DtTimeOp.hpp"
//#include "Pimpact_TimeNonlinearJacobianOp.hpp"

#include "Pimpact_TimeStokesBSmoother.hpp"
#include "Pimpact_TimeStokesLSmoother.hpp"

namespace {


typedef double S;
typedef int O;
const int d = 4;
const int dNC = 2;

bool testMpi = true;
bool output = false;
double eps = 3e-1;
int domain = 1;
int dim = 3;

int npx = 1;
int npy = 1;
int npz = 1;
int npf = 1;


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
      "output", "noutput", &output,
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
	clp.setOption(
	    "npx", &npx,
	    "" );
	clp.setOption(
	    "npy", &npy,
	    "" );
	clp.setOption(
	    "npz", &npz,
	    "" );
	clp.setOption(
	    "npf", &npf,
	    "" );

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );
  pl->set( "domain", domain );

  pl->set( "lx", 2. );
  pl->set( "ly", 2. );
  pl->set( "lz", 2. );

  pl->set( "dim", dim );

  pl->set("nx",  33 );
  pl->set("ny",  33 );
  pl->set("nz",  33 );

  pl->set("nf", 32 );

}

TEUCHOS_UNIT_TEST( TimeOperator, TimeStokesOperator ) {
        
	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	typedef Pimpact::TimeStokesOp<SpaceT> OpT;
	auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto op = Pimpact::create<OpT>( space );

	auto x = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<SpaceT> >( space ),
			Pimpact::createTimeField< Pimpact::ScalarField<SpaceT> >( space ));
	auto y = x->clone();
	auto y_cc = x->clone();
	auto error = x->clone();
	auto pulsatile = x->clone();

	double p = 1;
	double alpha = std::sqrt(pl->get<double>("alpha2"));

	//RHS
	Pimpact::initVectorTimeField( y->getVFieldPtr(), Pimpact::ConstVel_inX, p);

	// true solution //pl->get<int> ("Block Size");
	Pimpact::initVectorTimeField( pulsatile->getVFieldPtr(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha );

	// consistecy check
	op->apply(*pulsatile,*y_cc);
	error->add( 1., *y_cc, -1., *y );

	if( output ) {
		y->write();
		y_cc->write(100);
	}

	std::cout << "|| error || = " << error()->norm()/std::sqrt( error->getLength() ) << std::endl;

	TEST_EQUALITY( error()->norm()<eps, true );

}



TEUCHOS_UNIT_TEST( TimeOperator, TimeStokesBSmooth ) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	typedef Pimpact::TimeStokesOp<SpaceT> OpT;
	auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto op = Pimpact::create<OpT>( space );

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeStokesBSmoother<OpT>( op ));

	// test smoothing properties
	auto x = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<SpaceT> >( space ),
					       Pimpact::createTimeField< Pimpact::ScalarField<SpaceT> >( space ));
	auto y = x->clone();
	auto error = x->clone();
	auto true_sol = x->clone();	
	
	double p = 1;
	double alpha = std::sqrt(pl->get<double>("alpha2"));

	// RHS
	Pimpact::initVectorTimeField( y->getVFieldPtr(), Pimpact::ConstVel_inX, p);	

	// true solution //pl->get<int> ("Block Size");
	Pimpact::initVectorTimeField( true_sol->getVFieldPtr(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha );
		
	// test smoothing
	x->random();
	x->scale(10);
	
	// initial error
	error->add( 1., *x, -1., *true_sol );
	if( output )
		error->write();
	std::cout << "err: " << error->norm() << "\n";	
	
	// aprrox. solution
	bSmoother->apply(*y,*x);		

	// final error
	error->add( 1., *x, -1., *true_sol );
	if( output )
		error->write(100);
	std::cout << "err: " << error->norm() << "\n";
	// aprrox. solution2
	bSmoother->apply(*y,*x);

	// final error2
	error->add( 1., *x, -1., *true_sol );
	if( output )
		error->write(200);	
	std::cout << "err: " << error->norm() << "\n";
}

TEUCHOS_UNIT_TEST( TimeOperator, TimeStokesLSmooth ) {

        pl->set("npx", npx) ;
        pl->set("npy", npy) ;
        pl->set("npz", npz) ;
        pl->set("npf", npf) ;

        typedef Pimpact::TimeStokesOp<SpaceT> OpT;
        auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

        auto op = Pimpact::create<OpT>( space );

        auto bSmoother = Teuchos::rcp(new Pimpact::TimeStokesBSmoother<OpT>( op ));      
        auto lSmoother = Teuchos::rcp(new Pimpact::TimeStokesLSmoother<OpT>( op ));
        
       
	auto x = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<SpaceT> >( space ),
                                               Pimpact::createTimeField< Pimpact::ScalarField<SpaceT> >( space ));

	auto error = x->clone(); 
	x->random();
	auto x2 = x->clone();
	auto y = x->clone();
	y->random();
	
	bSmoother->apply(*y,*x);
        lSmoother->apply(*y,*x2,2);                                        
	
	error->add(1.,*x,-1.,*x2);
        //error->write();

        std::cout << "|| error || = " << error()->norm()/std::sqrt( error->getLength() ) << std::endl;

        TEST_EQUALITY( error()->norm()<eps, true );
	}


TEUCHOS_UNIT_TEST( TimeOperator, TimeStokesBSmooth_conv ) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	typedef Pimpact::TimeStokesOp<SpaceT> OpT;
	auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto op = Pimpact::create<OpT>( space );

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeStokesBSmoother<OpT>( op ));

	auto x = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<SpaceT> >( space ),
			Pimpact::createTimeField< Pimpact::ScalarField<SpaceT> >( space ));
	auto y = x->clone();
	auto error = x->clone();
	auto true_sol = x->clone();

	double p = 1;
	double alpha = std::sqrt(pl->get<double>("alpha2"));

	Pimpact::initVectorTimeField( y->getVFieldPtr(), Pimpact::ConstVel_inX, p);

	Pimpact::initVectorTimeField( true_sol->getVFieldPtr(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha );

	x->random();
	x->scale(10);

	for (int i = 1; i < 50; i++){
                error->add( 1., *x, -1., *true_sol );
                std::cout  << error->norm()/std::sqrt( error->getLength() ) << "\n";
		bSmoother->apply(*y,*x);
	
		if (i%5==0 &&  output)
			error->write(300+i*100);
	}

}

TEUCHOS_UNIT_TEST( TimeOperator, TimeStokesLSmooth_conv ) {

        pl->set("npx", npx) ;
        pl->set("npy", npy) ;
        pl->set("npz", npz) ;
        pl->set("npf", npf) ;

        typedef Pimpact::TimeStokesOp<SpaceT> OpT;
        auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

        auto op = Pimpact::create<OpT>( space );

        auto lSmoother = Teuchos::rcp(new Pimpact::TimeStokesLSmoother<OpT>( op ));

        auto x = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<SpaceT> >( space ),
                        Pimpact::createTimeField< Pimpact::ScalarField<SpaceT> >( space ));
        auto y = x->clone();
        auto error = x->clone();
        auto true_sol = x->clone();

        double p = 1;
        double alpha = std::sqrt(pl->get<double>("alpha2"));

        Pimpact::initVectorTimeField( y->getVFieldPtr(), Pimpact::ConstVel_inX, p);

        Pimpact::initVectorTimeField( true_sol->getVFieldPtr(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha );

        x->random();
        x->scale(10);

        for (int i = 1; i < 50; i++){
                error->add( 1., *x, -1., *true_sol );
                std::cout  << error->norm()/std::sqrt( error->getLength() ) << "\n";
                lSmoother->apply(*y,*x,2);

                if (i%5==0 &&  output)
                        error->write(300+i*100);
        }
}

} // end of namespace
                                        

