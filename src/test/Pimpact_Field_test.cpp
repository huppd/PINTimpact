#include <cmath>
#include <iostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_Fields.hpp"




namespace {


typedef double S;
typedef int O;
const int d = 3;
const int dNC = 4;

typedef Pimpact::Space<S,O,d,dNC>                SpaceT;
typedef typename Pimpact::ScalarField<SpaceT>    SF;
typedef typename Pimpact::VectorField<SpaceT>    VF;
typedef typename Pimpact::ModeField<SF>          MSF;
typedef typename Pimpact::ModeField<VF>          MVF;
typedef typename Pimpact::CompoundField<VF,SF>   CF;
typedef typename Pimpact::CompoundField<MVF,MSF> CMF;


bool testMpi = true;
double eps = 1e-6;

int dim = 3;
int domain = 1;

int npx = 2;
int npy = 2;
int npz = 2;

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
  clp.setOption( "npx", &npx, "" );
  clp.setOption( "npy", &npy, "" );
  clp.setOption( "npz", &npz, "" );

  pl->set( "lx", 12. );
  pl->set( "ly", 6. );
  pl->set( "lz", 3. );

  pl->set("nx", 33 );
  pl->set("ny", 17 );
  pl->set("nz", 9 );


}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, create_init_print, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );

  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  space->print();

  auto p = Pimpact::create<FType>(space);

  p->init( space->rankST() );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, MVF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, CF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, CMF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, InfNorm_and_init, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );

  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>(space);

  double norm;

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
    norm = p->norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( i/2., norm, eps );

  }

  // one test with infty-norm
  int rank;
  int size;
  double init;
  MPI_Comm_rank(space->comm(),&rank);
  MPI_Comm_size(space->comm(),&size);
  for( double i = 0.; i<10.1; ++i) {
    init = (size-1)*i-1.;
    init = (init<0)?-init:init;
    p->init(rank*i-1.);
    norm = p->norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( init, norm, eps );

  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, SF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, VF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, MSF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, MVF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, CF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, CMF ) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, OneNorm_and_init, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );

  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
//    TEST_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm) );
    TEST_FLOATING_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm), eps );

  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, SF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, MSF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, MVF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, CF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, CMF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, TwoNorm_and_init, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
    TEST_FLOATING_EQUALITY( std::sqrt( std::pow(i/2.,2)*p->getLength() ), p->norm(Belos::TwoNorm), eps );
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, MVF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, CF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, CMF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, dot, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);

  int Np = vel1->getLength();
  int Nq = vel2->getLength();
  double dot;

  TEST_EQUALITY( Np , Nq );
  int N = Np;

  vel1->init(0.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_EQUALITY( dot<eps, true );

  vel1->init(1.);
  vel2->init(1.);
  dot = vel2->dot(*vel1);
  TEST_FLOATING_EQUALITY( (S)N, dot, eps );

  vel1->init(2.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY( 2.*N, dot, eps );

  vel1->init(1.);
  vel2->init(2.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY( 2.*N, dot, eps );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MVF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CMF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, scale, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>(space);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->scale(2.);
  norm = p->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(4*N), norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MVF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CMF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, random, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>(space);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->random();
  norm = p->norm(Belos::TwoNorm);
  TEST_INEQUALITY( N, norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MVF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CMF )

	

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TemplateField, add, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);
  auto vel3 = Pimpact::create<FType>(space);

  TEST_EQUALITY( vel1->getLength(), vel2->getLength() )
  TEST_EQUALITY( vel2->getLength(), vel3->getLength() )
  TEST_EQUALITY( vel1->getLength(), vel3->getLength() )

  double norm;
  int N = vel1->getLength();

  vel1->init(0.);
  vel2->init(1./2.);
  vel3->init(1./3.);

  vel1->add( 2., *vel2, 0., *vel3);
  norm = vel1->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm )

  vel1->init(0.);
  vel2->init(1./2.);
  vel3->init(1./3.);

  vel1->add( 0., *vel2, 3., *vel3);
  norm = vel1->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm )

  vel1->init(0.);
  vel2->init(1.);
  vel3->init(1.);

  vel1->add( 0.5, *vel2, 0.5, *vel3);
  norm = vel1->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm )

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MVF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CMF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, write, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>( space );

  p->init(1.);
  p->write();

  p->random();
  p->write(1);

  TEST_EQUALITY( 0, 0)
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MVF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CMF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, write_restart, FType ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::create<FType>(space);

  p->init(1.);
  p->write();

  p->random();
  p->write( 99, true );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write_restart, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write_restart, VF )

	

TEUCHOS_UNIT_TEST( ScalarField, initField ) {

  pl->set( "dim", dim );
  pl->set( "domain", domain );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", npy );
  pl->set("npz", (2==dim)?1:npz );

  auto space = Pimpact::createSpace( pl );

  auto x = Pimpact::createScalarField( space );

  for( int i=0; i<=7; ++i ) {
    x->initField( Pimpact::EScalarField(i) );
    x->write( i );
  }

}


} // namespace

