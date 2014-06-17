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
bool isImpactInit=false;

typedef double S;
typedef int O;

typedef Pimpact::VectorField<S,O,4> VF;
typedef Pimpact::ScalarField<S,O,4> SF;

typedef Pimpact::TimeField<VF> TVF;
typedef Pimpact::TimeField<SF> TSF;

typedef Pimpact::MultiField<TVF> MTVF;
typedef Pimpact::MultiField<TSF> MTSF;
//typedef Pimpact::VectorField<S,O,4> VF;


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
TEUCHOS_UNIT_TEST( TimeFieldVector, creator ) {

  int rank = Pimpact::init_impact_pre();

  auto ds = Pimpact::createDomainSize<S>( 1., 1., 1. );
  ds->print();

  auto bc = Pimpact::createBoudaryConditionsGlobal( Pimpact::Dirichelt2DChannel );
  bc->print();

  auto gs = Pimpact::createGridSizeGlobal<O>( 33, 33, 2, 8 );
  gs->print();

  auto pgs = Pimpact::createProcGridSize<O>(2,1,1,2);
  pgs->print();

  auto lgs = Pimpact::createGridSizeLocal<O,4>( gs,pgs );
  lgs->print();

  Pimpact::init_impact_mid();

  auto pg = Pimpact::createProcGrid<O,4>( lgs, bc, pgs );
  pg->print();

  Pimpact::init_impact_postpost();

  auto fS = Pimpact::createFieldSpace<O,4>();

  auto iS = Pimpact::createScalarIndexSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto space = Pimpact::createSpace<O,4>( fS, iS, iIS, fIS, gs, lgs, pgs, pg );

  space->print();
  //---------------------------------------------------
  auto field1 = Pimpact::createTimeField<VF>( space );

  auto field2 = field1->clone();

  std::cout << "field1: length: " << field1->getLength() << "\n";
  std::cout << "field2: length: " << field2->getLength() << "\n";

  auto mvec = Pimpact::createMultiField<TVF>(*field1,5);

  //   Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM =
      Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<S,MTVF>( MyOM, mvec );
  TEST_EQUALITY_CONST(res,true);

//  field2->random();
//  MPI_Barrier( MPI_COMM_WORLD );
//  field2->init( space->procCoordinate()[3] );
//  field2->norm();
//  field2->exchange();
//  field2->write();

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


  auto op = Pimpact::createTimeOpWrap<Pimpact::HelmholtzOp<S,O,4> >();

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

  auto adv = Pimpact::createTimeNonlinearJacobian();
  adv->assignField( *field1 );
  adv->apply( *field1, *field2 );

}

} // end of namespace
