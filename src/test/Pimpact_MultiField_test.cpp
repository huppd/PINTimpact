#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Pimpact_Fields.hpp"




namespace {



typedef double S;
typedef int O;

const int d = 3;
const int dNC=4;

typedef typename Pimpact::Space<S,O,d,dNC>    SpaceT;
typedef typename Pimpact::ScalarField<SpaceT> SF;
typedef typename Pimpact::VectorField<SpaceT> VF;
typedef typename Pimpact::ModeField<SF>       MSF;
typedef typename Pimpact::ModeField<VF>       MVF;


bool testMpi = true;
double eps = 1e-6;
int domain = 0;
int dim = 3;

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

  pl->set("nx", 33 );
  pl->set("ny", 17 );
  pl->set("nz", 9 );

}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, constructor, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>(space);

  auto mv = Pimpact::createMultiField( *p, 10 );

  const int m = mv->getNumberVecs();

  TEST_EQUALITY( 10, m );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, constructor, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, constructor, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, constructor, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, constructor, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, TwoNorm_and_init, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv = Pimpact::createMultiField( *p, 10 );

  const O m = mv->getNumberVecs();
  const O n = mv->getLength();
  std::vector<S> normval(m);

  // test different float values, assures that initial and norm work smoothly
  for( S i=0.; i< 200.1; ++i ) {
    mv->init( i/2. );
    mv->norm( normval, Belos::TwoNorm );
    for( int j=0; j<m; ++j )
      TEST_FLOATING_EQUALITY( std::sqrt(std::pow(i/2.,2)*n), normval[j], eps );
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, TwoNorm_and_init, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, TwoNorm_and_init, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, TwoNorm_and_init, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, TwoNorm_and_init, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, clone, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv = Pimpact::createMultiField( *p, 1 );

  auto mv2 = mv->clone(10);

  int n1(mv->getNumberVecs());
  int n2(mv2->getNumberVecs());

  TEST_EQUALITY( 1, n1 );
  TEST_EQUALITY( 10, n2 );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, clone, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, clone, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, clone, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, clone, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, CloneCopy, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv = Pimpact::createMultiField(*p,10);

  mv->random();
  auto mv2 = mv->CloneCopy();

  int n1(mv->getNumberVecs());
  int n2(mv2->getNumberVecs());

  TEST_EQUALITY( 10, n1 );
  TEST_EQUALITY( n1, n2 );

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv->norm(norm1);
  mv2->norm(norm2);
  for( int i=0; i<n1; ++i)
    TEST_EQUALITY( norm1[i], norm2[i] );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy, MVF )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, CloneCopy2, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv = Pimpact::createMultiField(*p,10);

  mv->random();

  std::vector<int> index(5);
  for(int i=0; i<5; ++i)
    index[i] = 2*i;

  auto mv2 = mv->CloneCopy(index);

  unsigned int n1 = (mv->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY( 10, n1 );
  TEST_EQUALITY( 5, n2 );
  TEST_EQUALITY( index.size(), n2 );

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv->norm(norm1);
  mv2->norm(norm2);

  for( unsigned int i=0; i<index.size(); ++i)
    TEST_EQUALITY( norm1[index[i]], norm2[i] );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy2, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy2, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy2, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy2, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, CloneCopy3, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv1 = Pimpact::createMultiField(*p,10);

  mv1->random();

  Teuchos::Range1D index(2,7);

  auto mv2 = mv1->CloneCopy(index);

  unsigned int n1(mv1->getNumberVecs());
  unsigned int n2(mv2->getNumberVecs());

  TEST_EQUALITY( 10, n1 );
  TEST_EQUALITY( index.size(), n2 );

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv1->norm(norm1);
  mv2->norm(norm2);
  for( int i=0; i<index.size(); ++i)
    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy3, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy3, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy3, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneCopy3, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, CloneViewNonConst1, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv1 = Pimpact::createMultiField(*p,10);

  mv1->init(0.);

  std::vector<int> index(5);
  for(int i=0; i<5; ++i)
    index[i] = 2*i;

  auto mv2 = mv1->CloneViewNonConst(index);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY( 10, n1 );
  TEST_EQUALITY( 5, n2 );
  TEST_EQUALITY( index.size(), n2 );

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv2->random();

  mv1->norm(norm1);
  mv2->norm(norm2);

  for( unsigned int i=0; i<index.size(); ++i)
    TEST_EQUALITY( norm1[index[i]], norm2[i] );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneViewNonConst1, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneViewNonConst1, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneViewNonConst1, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneViewNonConst1, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, CloneViewNonConst2, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv1 = Pimpact::createMultiField(*p,10);

  mv1->init(0.);

  Teuchos::Range1D index(2,7);

  auto mv2 = mv1->CloneViewNonConst(index);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY( 10, n1 );
  TEST_EQUALITY( index.size(), n2 );

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv2->random();

  mv1->norm(norm1);
  mv2->norm(norm2);

  for( unsigned int i=0; i<index.size(); ++i)
    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneViewNonConst2, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneViewNonConst2, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneViewNonConst2, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneViewNonConst2, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, CloneView1, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv1 = Pimpact::createMultiField(*p,10);

  mv1->init(0.);

  std::vector<int> index(5);
  for(int i=0; i<5; ++i)
    index[i] = 2*i;

  auto mv2 = mv1->CloneView(index);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY( 10, n1 );
  TEST_EQUALITY( 5, n2 );
  TEST_EQUALITY( index.size(), n2 );

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  //			mv2->Random(); //< this should give compile error
  mv1->random();

  mv1->norm(norm1);
  mv2->norm(norm2);

  for( unsigned int i=0; i<index.size(); ++i)
    TEST_EQUALITY( norm1[index[i]], norm2[i] );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneView1, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneView1, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneView1, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneView1, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, CloneView2, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv1 = Pimpact::createMultiField(*p,10);

  mv1->init(0.);

  Teuchos::Range1D index(2,7);

  auto mv2 = mv1->CloneView(index);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY( 10, n1 );
  TEST_EQUALITY( index.size(), n2 );

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  //	mv2->Random(); // has to give compilation error
  mv1->random();

  mv1->norm(norm1);
  mv2->norm(norm2);

  for( unsigned int i=0; i<index.size(); ++i)
    TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneView2, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneView2, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneView2, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, CloneView2, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, TimesMatAdd, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv1 = Pimpact::createMultiField( *p, 10 );

  mv1->init(0.);

  Teuchos::Range1D index1(0,9);
  std::vector<int> index2(10);
  for(int i=0; i<10; ++i)
    index2[i] = i;

  auto mv2 = mv1->CloneView(index1);
  auto mv3 = mv1->clone(mv1->getNumberVecs() );

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());
  unsigned int n3 = (mv3->getNumberVecs());

  TEST_EQUALITY( n1, n2 );
  TEST_EQUALITY( n3, n2 );
  TEST_EQUALITY( n3, n1 );

  //  mv2->Random(); // has to give compilation error
  mv1->init(1.);
  //  mv2->Random();
  //	mv3->Assign(2.);

  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);

  for( unsigned int j=0; j<n1; ++j)
    for( unsigned int i=0; i<n1; ++i)
      B(j,i) = 1./n1;

  mv1->TimesMatAdd( 0.5, *mv2, B, 0.5 );

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv1->norm(norm1);
  mv2->norm(norm2);

  for( unsigned int i=0; i<n1; ++i) {
    TEST_FLOATING_EQUALITY( norm1[i], norm2[i], eps  );
    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], 3.e-1 );
  }

  std::vector<double> scales(n1);
  for( unsigned int j=0; j<n1; ++j){
    scales[j] = (j+1);
    for( unsigned int i=0; i<n1; ++i)
      B(j,i) = 1./n1/(j+1);
  }
  mv1->init(1.);
  mv3->init(1.);

  mv3->scale(scales);

  mv1->TimesMatAdd( 1., *mv3, B, 0. );

  mv1->norm(norm1);
  mv2->norm(norm2);

  for( unsigned int i=0; i<n1; ++i) {
    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], eps );
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, TimesMatAdd, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, TimesMatAdd, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, TimesMatAdd, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, TimesMatAdd, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, add, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv1 = Pimpact::createMultiField(*p,10);

  mv1->init(0.);

  Teuchos::Range1D index1(0,9);
  std::vector<int> index2(10);
  for(int i=0; i<10; ++i)
    index2[i] = i;

  auto mv2 = mv1->CloneViewNonConst(index1);
  auto mv3 = mv1->CloneView(index1);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());
  unsigned int n3 = (mv3->getNumberVecs());

  TEST_EQUALITY( n1, n2 );
  TEST_EQUALITY( n3, n2 );
  TEST_EQUALITY( n3, n1 );

  mv1->init(1.);
  mv2->init(1.);

  mv1->add( 0.5, *mv2, 0.5, *mv3);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv1->norm(norm1);
  mv2->norm(norm2);

  for( unsigned int i=0; i<n1; ++i)
    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm2[i], eps );

  mv1->init(1.);
  mv2->init(1.);

  mv2->scale(0.5);

  mv1->add( 1., *mv2, 1., *mv3 );

  mv1->norm(norm1);
  mv2->assign(*mv1);
  mv2->norm(norm2);

  for( unsigned int i=0; i<n1; ++i) {
    TEST_FLOATING_EQUALITY( norm1[i], norm2[i], eps  );
    TEST_FLOATING_EQUALITY( std::sqrt((double)mv1->getLength()), norm1[i], eps );
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, add, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, add, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, add, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, add, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, dot, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv1 = Pimpact::createMultiField(*p,10);

  mv1->init(0.);

  Teuchos::Range1D index1(0,9);
  std::vector<int> index2(10);
  for(int i=0; i<10; ++i)
    index2[i] = i;

  auto mv2 = mv1->CloneViewNonConst(index1);
  auto mv3 = mv1->CloneView(index1);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());
  unsigned int n3 = (mv3->getNumberVecs());

  TEST_EQUALITY( n1, n2 );
  TEST_EQUALITY( n3, n2 );
  TEST_EQUALITY( n3, n1 );

  mv1->init(1.);
  mv2->init(1.);

  std::vector<double> dots(n1);

  mv1->dot( *mv2, dots );


  for( unsigned int i=0; i<n1; ++i) {
    TEST_EQUALITY( mv1->getLength(), dots[i] );
  }

  mv2->init(2.);

  mv3->dot( *mv2, dots  );

  for( unsigned int i=0; i<n1; ++i) {
    TEST_EQUALITY( 4*mv1->getLength(), dots[i] );
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, dot, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, dot, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, dot, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, dot, MVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiField, Trans, FType ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

  // processor grid size
  pl->set("npx", (2==dim)?4:2 );
  pl->set("npy",            2 );
  pl->set("npz", (2==dim)?1:2 );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  auto mv1 = Pimpact::createMultiField(*p,10);

  mv1->init(0.);

  Teuchos::Range1D index1(0,9);
  std::vector<int> index2(10);
  for(int i=0; i<10; ++i)
    index2[i] = i;

  auto mv2 = mv1->CloneView(index1);
  auto mv3 = mv1->CloneView(index1);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());
  unsigned int n3 = (mv3->getNumberVecs());

  TEST_EQUALITY( n1, n2 );
  TEST_EQUALITY( n3, n2 );
  TEST_EQUALITY( n3, n1 );

  mv1->init(1.);

  Teuchos::SerialDenseMatrix<int,double> B(n1,n2);

  mv1->Trans( 1., *mv2, B );

  for( unsigned int j=0; j<n1; ++j) {
    for( unsigned int i=0; i<n1; ++i)
      TEST_EQUALITY( mv1->getLength(), B(j,i) );
  }

  std::vector<double> scales(n1);

  for( unsigned int i=0; i<scales.size(); ++i)
    scales[i] = i*2;

  mv1->scale(scales);

  mv2->Trans( 1., *mv3, B );

  for( unsigned int j=0; j<n1; ++j) {
    for( unsigned int i=0; i<n1; ++i)
      TEST_EQUALITY( scales[i]*scales[j]*mv1->getLength(), B(j,i) );
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, Trans, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, Trans, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, Trans, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiField, Trans, MVF )




} // end of namespace
