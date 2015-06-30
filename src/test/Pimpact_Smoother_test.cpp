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
#include "Pimpact_TimeNonlinearJacobianOp.hpp"

#include "Pimpact_TimeStokesBSmoother.hpp"

namespace {


typedef double S;
typedef int O;
const int d = 4;
const int dNC = 2;

bool testMpi = true;
double eps = 3e-1;
int domain = 1;
int dim = 3;

int npx = 2;
int npy = 2;
int npz = 1;
int npf = 2;


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

  pl->set( "lx", 1. );
  pl->set( "ly", 1. );
  pl->set( "lz", 1. );

  pl->set( "dim", dim );

  pl->set("nx", 9 );
  pl->set("ny", 9 );
  pl->set("nz", 9 );

  pl->set("nf", 8 );

}

TEUCHOS_UNIT_TEST( TimeOperator, TimeStokesBSmooth ) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	typedef Pimpact::TimeStokesOp<SpaceT> OpT;
	auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

	auto op = Pimpact::create<OpT>( space );

        //auto x = Pimpact::create<typename OpT::DomainFieldT>( space );
        //auto y = Pimpact::create<typename OpT::RangeFieldT>( space );

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeStokesBSmoother<OpT>( op ));

	//bSmoother->apply(*x,*y);

	// zero test
	//TEST_EQUALITY( y->norm()<eps, true );

	// test smoothing properties
	auto x = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<SpaceT> >( space ),
					       Pimpact::createTimeField< Pimpact::ScalarField<SpaceT> >( space ));
	auto y = x->clone();
	auto y_cc = x->clone(); // for the consistecy check
	auto true_sol = x->clone();
	auto error = x->clone();

	// RHS
	Pimpact::initVectorTimeField( y->getVFieldPtr(), Pimpact::ConstVel_inX,1.,0.,0. );	

	// true solution
	Pimpact::initVectorTimeField( true_sol->getVFieldPtr(), Pimpact::Pulsatile_inX,1.,1.,0. );
	
	// consistecy check
	op->apply(*y_cc,*true_sol);
	error->add( 1., *y_cc, -1., *true_sol );
	error->write();
	
	// initial error
	error->add( 1., *x, -1., *true_sol );
	error->write(100);

	// aprrox. solution
        bSmoother->apply(*x,*y);		

	// final error
	error->add( 1., *x, -1., *true_sol );
	error->write(200);
	
	//true_sol->write();
	//x->write(1);
}


} // end of namespace
