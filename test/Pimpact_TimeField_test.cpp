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
  // starting with ininializing
//  int rank = Pimpact::init_impact_pre();
//
//  auto ds = Pimpact::createDomainSize<S>( 1., 1., 1. );
//
//  auto bc = Pimpact::createBoudaryConditionsGlobal( Pimpact::Dirichelt2DChannel );
//
//  auto gs = Pimpact::createGridSizeGlobal<O>( 33, 33, 33 );
//
//  auto pgs = Pimpact::createProcGridSize<O>(2,2,1);
//
//  auto lgs = Pimpact::createGridSizeLocal<O,3>( gs,pgs );
//
//  Pimpact::init_impact_mid();
//
//  auto pg = Pimpact::createProcGrid<O,3>( lgs, bc, pgs );
//
//  Pimpact::init_impact_postpost();
  int rank = Pimpact::init_impact_pre();

  auto ds = Pimpact::createDomainSize<S>( 1., 1., 1. );
  ds->print();

  auto bc = Pimpact::createBoudaryConditionsGlobal( Pimpact::Dirichelt2DChannel );
  bc->print();

  auto gs = Pimpact::createGridSizeGlobal<O>( 33, 33, 2, 8 );
  gs->print();

  auto pgs = Pimpact::createProcGridSize<O>(1,1,1,4);
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


  auto op = Pimpact::createTimeOpWrap<Pimpact::Helmholtz<S,O,4> >();

  field1->random();
  op->apply( *field1, *field2 );
  field2->write();

}

//TEUCHOS_UNIT_TEST( MultiFieldScalarScalar, TwoNorm_and_init ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>( space );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  const int m = mv->getNumberVecs();
//  const int n = mv->getLength();
//  std::vector<double> normval(m);
//
//  // test different float values, assures that initial and norm work smoothly
//  for( double i=0.; i< 200.1; ++i ) {
//    mv->init(i/2.);
//    mv->norm(normval,Belos::TwoNorm);
//    for( int j=0; j<m; ++j )
//      TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*n), normval[j], errorTolSlack );
//  }
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, clone ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>( space );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,1);
//
//  auto mv2 = mv->clone(10);
//
//  int n1(mv->getNumberVecs());
//  int n2(mv2->getNumberVecs());
//
//  TEST_EQUALITY( 1, n1 );
//  TEST_EQUALITY( 10, n2 );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, CloneCopy ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>( space );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv->random();
//  auto mv2 = mv->CloneCopy();
//
//  int n1(mv->getNumberVecs());
//  int n2(mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( n1, n2 );
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv->norm(norm1);
//  mv2->norm(norm2);
//  for( int i=0; i<n1; ++i)
//    TEST_EQUALITY( norm1[i], norm2[i] );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, CloneCopy2 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>( space );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv->random();
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv->CloneCopy(index);
//
//  unsigned int n1 = (mv->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, CloneCopy3 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>( space );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv1->random();
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneCopy(index);
//
//  unsigned int n1(mv1->getNumberVecs());
//  unsigned int n2(mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//  for( int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, CloneViewNonConst1 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>( space );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv1->init(0.);
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv1->CloneViewNonConst(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv2->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, CloneViewNonConst2 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>( space );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneViewNonConst(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv2->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, CloneView1 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>( space );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv1->init(0.);
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv1->CloneView(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  //			mv2->Random(); //< this should give compile error
//  mv1->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, CloneViewt2 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>( space );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneView(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  //	mv2->Random(); // has to give compilation error
//  mv1->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, TimesMatAdd ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>(sVS);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneView(index1);
//  auto mv3 = mv1->clone(mv1->getNumberVecs() );
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  //  mv2->Random(); // has to give compilation error
//  mv1->init(1.);
//  //  mv2->Random();
//  //	mv3->Assign(2.);
//
//  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//
//  for( unsigned int j=0; j<n1; ++j)
//    for( unsigned int i=0; i<n1; ++i)
//      B(j,i) = 1./n1;
//
//  mv1->TimesMatAdd( 0.5, *mv2, B, 0.5 );
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_FLOATING_EQUALITY( norm1[i], norm2[i], errorTolSlack  );
//    TEST_FLOATING_EQUALITY( std::sqrt(mv1->getLength()), norm2[i], errorTolSlack );
//  }
//
//  std::vector<double> scales(n1);
//  for( unsigned int j=0; j<n1; ++j){
//    scales[j] = (j+1);
//    for( unsigned int i=0; i<n1; ++i)
//      B(j,i) = 1./n1/(j+1);
//  }
//  mv1->init(1.);
//  mv3->init(1.);
//
//  mv3->scale(scales);
//
//  mv1->TimesMatAdd( 1., *mv3, B, 0. );
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//  }
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, add ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>(sVS);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneViewNonConst(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  mv1->add( 0.5, *mv2, 0.5, *mv3);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i)
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  mv2->scale(0.5);
//
//  mv1->add( 1., *mv2, 1., *mv3 );
//
//  mv1->norm(norm1);
//  mv2->assign(*mv1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_FLOATING_EQUALITY( norm1[i], norm2[i], errorTolSlack  );
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm1[i], errorTolSlack );
//  }
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, dot ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>(sVS);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneViewNonConst(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  std::vector<double> dots(n1);
//
//  mv1->dot( *mv2, dots );
//
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( mv1->getLength(), dots[i] );
//  }
//
//  mv2->init(2.);
//
//  mv3->dot( *mv2, dots  );
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( 4*mv1->getLength(), dots[i] );
//  }
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalar, Trans ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto p = Pimpact::createScalarField<double,int>(sVS);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneView(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  mv1->init(1.);
//
//  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//
//  mv1->Trans( 1., *mv2, B );
//
//  for( unsigned int j=0; j<n1; ++j) {
//    for( unsigned int i=0; i<n1; ++i)
//      TEST_EQUALITY( mv1->getLength(), B(j,i) );
//  }
//
//  std::vector<double> scales(n1);
//
//  for( unsigned int i=0; i<scales.size(); ++i)
//    scales[i] = i*2;
//
//  mv1->scale(scales);
//
//  mv2->Trans( 1., *mv3, B );
//
//  for( unsigned int j=0; j<n1; ++j) {
//    for( unsigned int i=0; i<n1; ++i)
//      TEST_EQUALITY( scales[i]*scales[j]*mv1->getLength(), B(j,i) );
//  }
//}
//
//
//// test shows that nLoc is not consistent with start and end indexes
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, constructor ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  const int m = mv->getNumberVecs();
//
//  TEST_EQUALITY( 10, m );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, TwoNorm_and_init ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  const int m = mv->getNumberVecs();
//  const int n = mv->getLength();
//  std::vector<double> normval(m);
//
//  // test different float values, assures that initial and norm work smoothly
//  for( double i=0.; i< 200.1; ++i ) {
//    mv->init(i/2.);
//    mv->norm(normval,Belos::TwoNorm);
//    for( int j=0; j<m; ++j )
//      TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*n), normval[j], errorTolSlack );
//  }
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, clone ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >( *vel, 1 );
//
//  auto mv2 = mv->clone(10);
//
//  int m1(mv->getNumberVecs());
//  int m2(mv2->getNumberVecs());
//
//  int n1(mv->getLength());
//  int n2(mv2->getLength());
//
//  TEST_EQUALITY( 1, m1 );
//  TEST_EQUALITY( 10, m2 );
//
//  TEST_EQUALITY( n1, n2);
//
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, CloneCopy ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv->random();
//  auto mv2 = mv->CloneCopy();
//
//  int n1(mv->getNumberVecs());
//  int n2(mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( n1, n2 );
//
//  int m1(mv->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv->norm(norm1);
//  mv2->norm(norm2);
//  for( int i=0; i<n1; ++i)
//    TEST_EQUALITY( norm1[i], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, CloneCopy2 ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv->random();
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv->CloneCopy(index);
//
//  unsigned int n1 = (mv->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, CloneCopy3 ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv1->random();
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneCopy(index);
//
//  unsigned int n1(mv1->getNumberVecs());
//  unsigned int n2(mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//  for( int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, CloneViewNonConst1 ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv1->CloneViewNonConst(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  //
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv2->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, CloneViewNonConst2 ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneViewNonConst(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv2->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, CloneView1 ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv1->CloneView(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, CloneViewt2 ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneView(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  //      mv2->Random(); // has to give compilation error
//  mv1->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, TimesMatAdd ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneView(index1);
//  auto mv3 = mv1->clone(mv1->getNumberVecs() );
//
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  int m3(mv3->getLength());
//  TEST_EQUALITY( m1, m2);
//  TEST_EQUALITY( m2, m3);
//
//
//  //      mv2->Random(); // has to give compilation error
//  mv1->init(1.);
//  //        mv2->Random();
//  //        mv3->Assign(2.);
//
//
//  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//  for( unsigned int j=0; j<n1; ++j)
//    for( unsigned int i=0; i<n1; ++i)
//      B(j,i) = 1./n1;
//
//  mv1->TimesMatAdd( 0.5, *mv2, B, 0.5 );
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_FLOATING_EQUALITY( norm1[i], norm2[i], errorTolSlack   );
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack  );
//  }
//
//  std::vector<double> scales(n1);
//  for( unsigned int j=0; j<n1; ++j){
//    scales[j] = (j+1);
//    //          scales[j] = 0;
//    for( unsigned int i=0; i<n1; ++i)
//      B(j,i) = 1./n1/(j+1);
//    //            B(j,i) = 1./5.;
//  }
//  //        std::cout << B;
//  //        scales[5] = 5.;
//  mv1->init(1.);
//  mv3->init(1.);
//
//  mv3->scale(scales);
//
//  mv1->TimesMatAdd( 1., *mv3, B, 0. );
//
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    //          TEST_EQUALITY( norm1[i], norm2[i] );
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//  }
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, add ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  //          auto mv2 = mv1->CloneCopy(index1);
//  auto mv2 = mv1->CloneViewNonConst(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  int m3(mv3->getLength());
//  TEST_EQUALITY( m1, m2);
//  TEST_EQUALITY( m2, m3);
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  mv1->add( 0.5, *mv2, 0.5, *mv3);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i)
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  mv2->scale(0.5);
//
//  mv1->add( 1., *mv2, 1., *mv3 );
//
//
//  mv1->norm(norm1);
//  mv2->assign(*mv1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( norm1[i], norm2[i] );
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm1[i], errorTolSlack );
//  }
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, dot ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  //          auto mv2 = mv1->CloneCopy(index1);
//  auto mv2 = mv1->CloneViewNonConst(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  int m3(mv3->getLength());
//
//  TEST_EQUALITY( m1, m2);
//  TEST_EQUALITY( m2, m3);
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  std::vector<double> dots(n1);
//
//  mv1->dot( *mv2, dots );
//
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( mv1->getLength(), dots[i] );
//  }
//
//  mv2->init(2.);
//
//  mv3->dot( *mv2, dots  );
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( 4*mv1->getLength(), dots[i] );
//  }
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldScalarMode, Trans ) {
//
//  auto sVS = Pimpact::createSpace();
//
//  auto pc = Pimpact::createScalarField<double,int>(sVS);
//  auto ps = Pimpact::createScalarField<double,int>(sVS);
//
//  auto vel = Pimpact::createModeField( pc, ps );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneView(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//
//  mv1->init(1.);
//
//  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//
//  mv1->Trans( 1., *mv2, B );
//
//  for( unsigned int j=0; j<n1; ++j){
//    for( unsigned int i=0; i<n1; ++i)
//      TEST_EQUALITY( mv1->getLength(), B(j,i) );
//  }
//
//  std::vector<double> scales(n1);
//  for( unsigned int i=0; i<scales.size(); ++i)
//    scales[i] = i*2;
//  mv1->scale(scales);
//
//  mv2->Trans( 1., *mv3, B );
//
//  for( unsigned int j=0; j<n1; ++j){
//    for( unsigned int i=0; i<n1; ++i)
//      TEST_EQUALITY( scales[i]*scales[j]*mv1->getLength(), B(j,i) );
//  }
//}
//
//
//
//// test shows that nLoc is not consistent with start and end indexes
//TEUCHOS_UNIT_TEST( MultiFieldVector, constructor ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  const int m = mv->getNumberVecs();
//
//  TEST_EQUALITY( 10, m );
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, TwoNorm_and_init ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  const int m = mv->getNumberVecs();
//  const int n = mv->getLength();
//  std::vector<double> normval(m);
//
//  // test different float values, assures that initial and norm work smoothly
//  for( double i=0.; i< 200.1; ++i ) {
//    mv->init(i/2.);
//    mv->norm(normval,Belos::TwoNorm);
//    for( int j=0; j<m; ++j )
//      TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*n), normval[j], errorTolSlack );
//  }
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, clone ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,1);
//
//  auto mv2 = mv->clone(10);
//
//  int m1(mv->getNumberVecs());
//  int m2(mv2->getNumberVecs());
//
//  int n1(mv->getLength());
//  int n2(mv2->getLength());
//
//  TEST_EQUALITY( 1, m1 );
//  TEST_EQUALITY( 10, m2 );
//
//  TEST_EQUALITY( n1, n2);
//
//}
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, CloneCopy ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv->random();
//  auto mv2 = mv->CloneCopy();
//
//  int n1(mv->getNumberVecs());
//  int n2(mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( n1, n2 );
//
//  int m1(mv->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv->norm(norm1);
//  mv2->norm(norm2);
//  for( int i=0; i<n1; ++i)
//    TEST_EQUALITY( norm1[i], norm2[i] );
//}
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, CloneCopy2 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv->random();
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv->CloneCopy(index);
//
//  unsigned int n1 = (mv->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//}
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, CloneCopy3 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv1->random();
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneCopy(index);
//
//  unsigned int n1(mv1->getNumberVecs());
//  unsigned int n2(mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//  for( int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//
//}
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, CloneViewNonConst1 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv1->init(0.);
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv1->CloneViewNonConst(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  //
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv2->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, CloneViewNonConst2 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneViewNonConst(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv2->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//}
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, CloneView1 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv1->init(0.);
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv1->CloneView(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  //
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  //      mv2->Random(); //< this should give compile error
//
//  mv1->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//}
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, CloneViewt2 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneView(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  //      mv2->Random(); // has to give compilation error
//  mv1->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//}
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, TimesMatAdd ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneView(index1);
//  auto mv3 = mv1->clone(mv1->getNumberVecs() );
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  int m3(mv3->getLength());
//
//  TEST_EQUALITY( m1, m2);
//  TEST_EQUALITY( m2, m3);
//
//  mv1->init(1.);
//
//  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//
//  for( unsigned int j=0; j<n1; ++j)
//    for( unsigned int i=0; i<n1; ++i)
//      B(j,i) = 1./n1;
//
//  mv1->TimesMatAdd( 0.5, *mv2, B, 0.5 );
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( norm1[i], norm2[i] );
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//  }
//
//  std::vector<double> scales(n1);
//  for( unsigned int j=0; j<n1; ++j){
//    scales[j] = (j+1);
//    for( unsigned int i=0; i<n1; ++i)
//      B(j,i) = 1./n1/(j+1);
//  }
//  mv1->init(1.);
//  mv3->init(1.);
//
//  mv3->scale(scales);
//
//  mv1->TimesMatAdd( 1., *mv3, B, 0. );
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i)
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, add ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneViewNonConst(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  int m3(mv3->getLength());
//
//  TEST_EQUALITY( m1, m2);
//  TEST_EQUALITY( m2, m3);
//
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  mv1->add( 0.5, *mv2, 0.5, *mv3);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i)
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  mv2->scale(0.5);
//
//  mv1->add( 1., *mv2, 1., *mv3 );
//
//  mv1->norm(norm1);
//  mv2->assign(*mv1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( norm1[i], norm2[i] );
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm1[i], errorTolSlack );
//  }
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, dot ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  //          auto mv2 = mv1->CloneCopy(index1);
//  auto mv2 = mv1->CloneViewNonConst(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  int m3(mv3->getLength());
//  TEST_EQUALITY( m1, m2);
//  TEST_EQUALITY( m2, m3);
//
//
//
//  //      mv2->Random(); // has to give compilation error
//  mv1->init(1.);
//  mv2->init(1.);
//  //        mv2->Random();
//  //        mv3->Assign(2.);
//
//
//  std::vector<double> dots(n1);
//
//  mv1->dot( *mv2, dots );
//
//
//  for( unsigned int i=0; i<n1; ++i) {
//    //            TEST_EQUALITY( norm1[i], norm2[i] );
//    TEST_EQUALITY( mv1->getLength(), dots[i] );
//  }
//
//  mv2->init(2.);
//
//  mv3->dot( *mv2, dots  );
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( 4*mv1->getLength(), dots[i] );
//  }
//}
//
//TEUCHOS_UNIT_TEST( MultiFieldVector, Trans ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto vel = Pimpact::createVectorField<double,int>(space);
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneView(index1);
//  auto mv3 = mv1->CloneView(index1);
//  //        auto mv3 = mv1->CloneV(mv1->GetNumberVecs() );
//
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//
//  //      mv2->Random(); // has to give compilation error
//  mv1->init(1.);
//  //        mv2->Random();
//  //        mv3->Assign(2.);
//
//
//  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//
//  mv1->Trans( 1., *mv2, B );
//
//  for( unsigned int j=0; j<n1; ++j){
//    for( unsigned int i=0; i<n1; ++i)
//      TEST_EQUALITY( mv1->getLength(), B(j,i) );
//  }
//
//  //        std::vector<double> norm1(n1);
//  //        std::vector<double> norm2(n2);
//  //
//  //        mv1->Norm(norm1);
//  //        mv2->Norm(norm2);
//  //
//  //        for( unsigned int i=0; i<n1; ++i) {
//  //          TEST_EQUALITY( norm1[i], norm2[i] );
//  //          TEST_EQUALITY( mv1->getLength(), norm2[i] );
//  //        }
//
//  std::vector<double> scales(n1);
//  for( unsigned int i=0; i<scales.size(); ++i)
//    scales[i] = i*2;
//  mv1->scale(scales);
//
//  mv2->Trans( 1., *mv3, B );
//
//  for( unsigned int j=0; j<n1; ++j){
//    for( unsigned int i=0; i<n1; ++i)
//      TEST_EQUALITY( scales[i]*scales[j]*mv1->getLength(), B(j,i) );
//  }
//}
//
//
//// test shows that nLoc is not consistent with start and end indexes
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, constructor ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space);
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  const int m = mv->getNumberVecs();
//
//  TEST_EQUALITY( 10, m );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, TwoNorm_and_init ) {
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space);
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  const int m = mv->getNumberVecs();
//  const int n = mv->getLength();
//  std::vector<double> normval(m);
//
//  // test different float values, assures that initial and norm work smoothly
//  for( double i=0.; i< 200.1; ++i ) {
//    mv->init(i/2.);
//    mv->norm(normval,Belos::TwoNorm);
//    for( int j=0; j<m; ++j )
//      TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*n), normval[j], errorTolSlack );
//  }
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, clone ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space);
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,1);
//
//
//  auto mv2 = mv->clone(10);
//
//  int m1(mv->getNumberVecs());
//  int m2(mv2->getNumberVecs());
//
//  int n1(mv->getLength());
//  int n2(mv2->getLength());
//
//  TEST_EQUALITY( 1, m1 );
//  TEST_EQUALITY( 10, m2 );
//
//  TEST_EQUALITY( n1, n2);
//
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, CloneCopy ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space);
//  auto vels = Pimpact::createVectorField<double,int>(space );
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv->random();
//  auto mv2 = mv->CloneCopy();
//
//  int n1(mv->getNumberVecs());
//  int n2(mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( n1, n2 );
//
//  int m1(mv->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv->norm(norm1);
//  mv2->norm(norm2);
//  for( int i=0; i<n1; ++i)
//    TEST_EQUALITY( norm1[i], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, CloneCopy2 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space );
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv->random();
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv->CloneCopy(index);
//
//  unsigned int n1 = (mv->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, CloneCopy3 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space);
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv1->random();
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneCopy(index);
//
//  unsigned int n1(mv1->getNumberVecs());
//  unsigned int n2(mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//  for( int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, CloneViewNonConst1 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space);
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv1->CloneViewNonConst(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  //
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv2->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, CloneViewNonConst2 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space );
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneViewNonConst(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv2->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, CloneView1 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space);
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  std::vector<int> index(5);
//  for(int i=0; i<5; ++i)
//    index[i] = 2*i;
//
//  auto mv2 = mv1->CloneView(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( 5, n2 );
//  TEST_EQUALITY( index.size(), n2 );
//
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[index[i]], norm2[i] );
//}
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, CloneViewt2 ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space);
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index(2,7);
//
//  auto mv2 = mv1->CloneView(index);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//
//  TEST_EQUALITY( 10, n1 );
//  TEST_EQUALITY( index.size(), n2 );
//
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  TEST_EQUALITY( m1, m2);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  //      mv2->random(); // has to give compilation error
//  mv1->random();
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<index.size(); ++i)
//    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, TimesMatAdd ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>(space);
//  auto vels = Pimpact::createVectorField<double,int>(space);
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneView(index1);
//  auto mv3 = mv1->clone(mv1->getNumberVecs() );
//
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  int m3(mv3->getLength());
//  TEST_EQUALITY( m1, m2);
//  TEST_EQUALITY( m2, m3);
//
//
//  mv1->init(1.);
//
//  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//  for( unsigned int j=0; j<n1; ++j)
//    for( unsigned int i=0; i<n1; ++i)
//      B(j,i) = 1./n1;
//
//  mv1->TimesMatAdd( 0.5, *mv2, B, 0.5 );
//
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( norm1[i], norm2[i] );
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//  }
//
//  std::vector<double> scales(n1);
//  for( unsigned int j=0; j<n1; ++j){
//    scales[j] = (j+1);
//    //          scales[j] = 0;
//    for( unsigned int i=0; i<n1; ++i)
//      B(j,i) = 1./n1/(j+1);
//    //            B(j,i) = 1./5.;
//  }
//  //        std::cout << B;
//  //        scales[5] = 5.;
//  mv1->init(1.);
//  mv3->init(1.);
//
//  mv3->scale(scales);
//
//  mv1->TimesMatAdd( 1., *mv3, B, 0. );
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//  }
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, add ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>( space );
//  auto vels = Pimpact::createVectorField<double,int>( space );
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  //          auto mv2 = mv1->CloneCopy(index1);
//  auto mv2 = mv1->CloneViewNonConst(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  int m3(mv3->getLength());
//  TEST_EQUALITY( m1, m2);
//  TEST_EQUALITY( m2, m3);
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  mv1->add( 0.5, *mv2, 0.5, *mv3);
//
//  std::vector<double> norm1(n1);
//  std::vector<double> norm2(n2);
//
//  mv1->norm(norm1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], errorTolSlack );
//  }
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  mv2->scale(0.5);
//
//  mv1->add( 1., *mv2, 1., *mv3 );
//
//  mv1->norm(norm1);
//  mv2->assign(*mv1);
//  mv2->norm(norm2);
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( norm1[i], norm2[i] );
//    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm1[i], errorTolSlack );
//  }
//
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, dot ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>( space );
//  auto vels = Pimpact::createVectorField<double,int>( space );
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  //          auto mv2 = mv1->CloneCopy(index1);
//  auto mv2 = mv1->CloneViewNonConst(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//  int m1(mv1->getLength());
//  int m2(mv2->getLength());
//  int m3(mv3->getLength());
//
//  TEST_EQUALITY( m1, m2);
//  TEST_EQUALITY( m2, m3);
//
//  mv1->init(1.);
//  mv2->init(1.);
//
//  std::vector<double> dots(n1);
//
//  mv1->dot( *mv2, dots );
//
//
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( mv1->getLength(), dots[i] );
//  }
//
//  mv2->init(2.);
//
//  mv3->dot( *mv2, dots  );
//  for( unsigned int i=0; i<n1; ++i) {
//    TEST_EQUALITY( 4*mv1->getLength(), dots[i] );
//  }
//}
//
//
//
//TEUCHOS_UNIT_TEST( MultiFieldVectorMode, Trans ) {
//
//  auto space = Pimpact::createSpace();
//
//  auto velc = Pimpact::createVectorField<double,int>( space );
//  auto vels = Pimpact::createVectorField<double,int>( space );
//
//  auto vel = Pimpact::createModeField( velc, vels );
//
//  auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);
//
//  mv1->init(0.);
//
//  Teuchos::Range1D index1(0,9);
//  std::vector<int> index2(10);
//  for(int i=0; i<10; ++i)
//    index2[i] = i;
//
//  auto mv2 = mv1->CloneView(index1);
//  auto mv3 = mv1->CloneView(index1);
//
//  unsigned int n1 = (mv1->getNumberVecs());
//  unsigned int n2 = (mv2->getNumberVecs());
//  unsigned int n3 = (mv3->getNumberVecs());
//
//  TEST_EQUALITY( n1, n2 );
//  TEST_EQUALITY( n3, n2 );
//  TEST_EQUALITY( n3, n1 );
//
//
//  mv1->init(1.);
//
//  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
//
//  mv1->Trans( 1., *mv2, B );
//
//  for( unsigned int j=0; j<n1; ++j){
//    for( unsigned int i=0; i<n1; ++i)
//      TEST_EQUALITY( mv1->getLength(), B(j,i) );
//  }
//
//  std::vector<double> scales(n1);
//  for( unsigned int i=0; i<scales.size(); ++i)
//    scales[i] = i*2;
//  mv1->scale(scales);
//
//  mv2->Trans( 1., *mv3, B );
//
//  for( unsigned int j=0; j<n1; ++j){
//    for( unsigned int i=0; i<n1; ++i)
//      TEST_EQUALITY( scales[i]*scales[j]*mv1->getLength(), B(j,i) );
//  }
//}


} // end of namespace
