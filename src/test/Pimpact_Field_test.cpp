#include <cmath>
#include <iostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_Fields.hpp"




namespace {


using ST = double;
using OT = int;

const int sd = 3;
const int d = 3;
const int dNC = 4;

using SpaceT = Pimpact::Space<ST,OT,sd,d,dNC>;

using Space2D = Pimpact::Space<ST,OT,2,d,dNC>;
using Space3D = Pimpact::Space<ST,OT,3,d,dNC>;

using SF2D  = typename Pimpact::ScalarField<Space2D>;
using VF2D  = typename Pimpact::VectorField<Space2D>;
using MSF2D = typename Pimpact::ModeField<SF2D>;
using MVF2D = typename Pimpact::ModeField<VF2D>;
using CF2D  = typename Pimpact::CompoundField<VF2D,SF2D>;
using CMF2D = typename Pimpact::CompoundField<MVF2D,MSF2D>;

using SF3D  = typename Pimpact::ScalarField<Space3D>;
using VF3D  = typename Pimpact::VectorField<Space3D>;
using MSF3D = typename Pimpact::ModeField<SF3D>;
using MVF3D = typename Pimpact::ModeField<VF3D>;
using CF3D  = typename Pimpact::CompoundField<VF3D,SF3D>;
using CMF3D = typename Pimpact::CompoundField<MVF3D,MSF3D>;

bool testMpi = true;
ST eps = 1e-6;

int domain = 1;

OT nx = 33;
OT ny = 17;
OT nz = 9;

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
	clp.setOption( "domain", &domain, "domain" );
  clp.setOption( "nx", &nx, "" );
  clp.setOption( "ny", &ny, "" );
  clp.setOption( "nz", &nz, "" );
  clp.setOption( "npx", &npx, "" );
  clp.setOption( "npy", &npy, "" );
  clp.setOption( "npz", &npz, "" );

  pl->set( "lx", 1. );
  pl->set( "ly", 1. );
  pl->set( "lz", 1. );

  pl->set("nx", nx );
  pl->set("ny", nx );
  pl->set("nz", nx );


}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, print, FType ) {

	Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );

  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

	//p->initField( Pimpact::Grad2D_inZ );
	p->init( space->rankST() );
	p->print();
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, InfNormAndInit, FType ) {

	Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  ST norm;

  // test different float values, assures that initial and norm work smoothly
  for( ST i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
    norm = p->norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( i/2., norm, eps );

  }

  // one test with infty-norm
  int rank;
  int size;
  ST init;
  MPI_Comm_rank(space->commST(),&rank);
  MPI_Comm_size(space->commST(),&size);
  for( ST i = 0.; i<10.1; ++i) {
    init = ( size-1 )*i-1.;
    init = std::abs( init );
    p->init( rank*i-1. );
    norm = p->norm( Belos::InfNorm );
    TEST_FLOATING_EQUALITY( std::max(init,1.), norm, eps );

  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, SF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, VF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, MSF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, MVF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, CF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, CMF2D ) 

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, SF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, VF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, MSF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, MVF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, CF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, CMF3D ) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, OneNormAndInit, FType ) {

	Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( ST i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
//    TEST_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm) );
    TEST_FLOATING_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm), eps );

  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, SF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, MSF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, CF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, SF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, MSF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, CF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, TwoNormAndInit, FType ) {

	Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( ST i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
    TEST_FLOATING_EQUALITY( std::sqrt( std::pow(i/2.,2)*p->getLength() ), p->norm(Belos::TwoNorm), eps );
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, dot, FType ) {

	Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);

  int Np = vel1->getLength();
  int Nq = vel2->getLength();
  ST dot;

  TEST_EQUALITY( Np , Nq );
  int N = Np;

  vel1->init(0.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_EQUALITY( dot<eps, true );

  vel1->init(1.);
  vel2->init(1.);
  dot = vel2->dot(*vel1);
  TEST_FLOATING_EQUALITY( static_cast<ST>(N), dot, eps );

  vel1->init(2.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY( 2.*N, dot, eps );

  vel1->init(1.);
  vel2->init(2.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY( 2.*N, dot, eps );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CMF3D )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, scale, FType ) {

	Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  ST norm;
  int N = p->getLength();

  p->init(1.);
  p->scale(2.);
  norm = p->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(4*N), norm)

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, random, FType ) {

	Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  ST norm;
  int N = p->getLength();

  p->init(1.);
  p->random();
  norm = p->norm(Belos::TwoNorm);
  TEST_INEQUALITY( N, norm)

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CMF3D )
	


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TemplateField, add, FType ) {

	Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);
  auto vel3 = Pimpact::create<FType>(space);

  TEST_EQUALITY( vel1->getLength(), vel2->getLength() )
  TEST_EQUALITY( vel2->getLength(), vel3->getLength() )
  TEST_EQUALITY( vel1->getLength(), vel3->getLength() )

  ST norm;
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

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, write, FType ) {

	Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>( space );

  p->init(1.);
  p->write();

  p->random();
  p->write(1);

  TEST_EQUALITY( 0, 0 )

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, writeRestart, FType ) {

  Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==FType::SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==FType::SpaceT::sdim)?1:npz );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  p->init(1.);
  p->write();

  p->random();
  p->write( 99, true );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, writeRestart, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, writeRestart, VF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, writeRestart, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, writeRestart, VF3D )
	

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ScalarField, initField, SpaceT ) {

  Pimpact::setBoundaryConditions( pl, domain );

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==SpaceT::sdim)?1:npz );

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

  auto x = Pimpact::createScalarField( space );

  for( int i=0; i<=7; ++i ) {
    x->initField( Pimpact::EScalarField(i) );
    x->write( i );
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ScalarField, initField, Space2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ScalarField, initField, Space3D )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ScalarField, level, SpaceT ) {

  Pimpact::setBoundaryConditions( pl, domain );
	// check: 0, 2, 4, 5
	// error: 1, 3

  pl->set("nx", nx );
  pl->set("ny", ny );
  pl->set("nz", nz );
  // processor grid size
  pl->set("npx", npx );
  pl->set("npy", (2==SpaceT::sdim)?npy*npz:npy );
  pl->set("npz", (2==SpaceT::sdim)?1:npz );

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

  auto x = Pimpact::createScalarField( space );

	x->init( 1. );
	x->level();

	ST level = x->norm();
	if( 0==space()->rankST() )
		std::cout << "\nlevel: " << level << "\n";
  TEST_EQUALITY( level<eps , true );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ScalarField, level, Space2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ScalarField, level, Space3D )

} // namespace

