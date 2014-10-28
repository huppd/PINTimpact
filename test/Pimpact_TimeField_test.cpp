// Pimpact_SalarVectorSpace_test.cpp

#include "mpi.h"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include "Teuchos_Range1D.hpp"
#include <Teuchos_CommHelpers.hpp>

#include "BelosOutputManager.hpp"
#include "BelosMVOPTester.hpp"

#include "BelosPimpactAdapter.hpp"


#include "pimpact.hpp"
#include "Pimpact_Space.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

#include "Pimpact_TimeField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_TimeOpWrap.hpp"
#include "Pimpact_DtTimeOp.hpp"
#include "Pimpact_TimeNonlinearJacobianOp.hpp"

#include <iostream>



namespace {


bool testMpi = true;
double errorTolSlack = 3e-1;
//bool isImpactInit=false;

typedef double S;
typedef int O;

typedef Pimpact::VectorField<S,O,4> VF;
typedef Pimpact::ScalarField<S,O,4> SF;

typedef Pimpact::TimeField<VF> TVF;
typedef Pimpact::TimeField<SF> TSF;

typedef Pimpact::MultiField<TVF> MTVF;
typedef Pimpact::MultiField<TSF> MTSF;


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
TEUCHOS_UNIT_TEST( TimeFieldVector, all ) {
  auto pl = Teuchos::parameterList();

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );
  pl->set( "domain", 1 );

  pl->set( "lx", 1. );
  pl->set( "ly", 1. );
  pl->set( "lz", 1. );

  pl->set( "dim", 2 );

  pl->set("nx", 33 );
  pl->set("ny", 33 );
  pl->set("nz", 2 );

  pl->set("nf", 32 );
//  pl->set("nfs", 0 );
//  pl->set("nfe", 0 );

  // processor grid size
  pl->set("npx", 2 );
  pl->set("npy", 2 );
  pl->set("npz", 1 );
  pl->set("npf", 2 );

  auto space = Pimpact::createSpace<S,O,4>( pl, true );

  space->print();
  //---------------------------------------------------
  auto field1 = Pimpact::createTimeField<VF>( space );

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

  auto mvec = Pimpact::createMultiField<TVF>(*field1,5);

  //   Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM =
      Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<S,MTVF>( MyOM, mvec );
  TEST_EQUALITY_CONST(res,true);


  auto fields = Pimpact::createTimeField<SF>( space );

  auto msca = Pimpact::createMultiField<TSF>(*fields,5);


  res = Belos::TestMultiVecTraits<S,MTSF>( MyOM, msca );
  TEST_EQUALITY_CONST(res,true);

  auto field2s = fields->clone();
//  field2s->random();
//  MPI_Barrier( MPI_COMM_WORLD );
//  field2s->norm();
  field2s->init( space->procCoordinate()[3] );
  field2s->exchange();
  field2s->write();


  auto op = Pimpact::createTimeOpWrap<Pimpact::HelmholtzOp<S,O,4> >(
      Pimpact::createHelmholtzOp<S,O,4>(space ),
      field2->getFieldPtr(0)->clone() );

  field1->random();
  op->apply( *field1, *field2 );
  Pimpact::initVectorTimeField<S,O>( field2, Pimpact::Zero2DFlow );
  field2->write();
  Pimpact::initVectorTimeField<S,O>( field2, Pimpact::Poiseuille_inX );
  field2->write(10);
  Pimpact::initVectorTimeField<S,O>( field2, Pimpact::Poiseuille_inY );
  field2->write(20);
  Pimpact::initVectorTimeField<S,O>( field2, Pimpact::Streaming2DFlow );
  field2->write(30);
  Pimpact::initVectorTimeField<S,O>( field2, Pimpact::OscilatingDisc2D );
  field2->write(40);

  auto dt = Pimpact::createDtTimeOp<S,O>();

//  field1->random();
  Pimpact::initVectorTimeField<S,O>( field1, Pimpact::OscilatingDisc2D );
  field2->init(0);

  dt->apply( *field1, *field2 );

  field2->write(50);

  auto adv = Pimpact::createTimeNonlinearJacobian( space );
  adv->assignField( *field1 );
  adv->apply( *field1, *field2 );

}

} // end of namespace
