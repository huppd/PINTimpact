#include <iostream>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_ConvectionOp.hpp"

#include "Pimpact_LinSolverParameter.hpp"




namespace {


typedef double S;
typedef int O;
const int d = 3;
const int dNC = 4;


typedef Pimpact::Space<S,O,d,dNC> SpaceT;
typedef typename Pimpact::ScalarField<SpaceT> SF;
typedef typename Pimpact::VectorField<SpaceT> VF;
typedef typename Pimpact::ModeField<SF>       MSF;
typedef typename Pimpact::ModeField<VF>       MVF;

template<class T> using ConvDiffOpT = Pimpact::ConvectionVOp<Pimpact::ConvectionDiffusionSOp<T> >;
template<class T> using ConvDiffSORT = Pimpact::ConvectionVSmoother<T,Pimpact::ConvectionDiffusionSORSmoother >;


bool testMpi = true;
S eps = 1e-10;
S omega = 0.8;
S winds = 1;
int nIter = 1000;

int dim = 2;
int domain = 1;

bool isImpactInit = false;

Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();



TEUCHOS_STATIC_SETUP() {

  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used." );
  clp.setOption(
      "eps", &eps,
      "Slack off of machine epsilon used to check test results" );
  clp.setOption(
      "omega", &omega,
      "Slack off of machine epsilon used to check test results" );
  clp.setOption(
      "wind", &winds,
      "Slack off of machine epsilon used to check test results" );
  clp.setOption(
      "nIter", &nIter,
      "Slack off of machine epsilon used to check test results" );
  clp.setOption(
      "domain", &domain,
      "domain" );
  clp.setOption(
      "dim", &dim,
      "dim" );

	pl->set( "dim", dim );
  pl->set( "domain", domain );

	pl->set("nx", 125 );
	//  pl->set("ny", 125 );
	//  pl->set("nx", 65 );
	pl->set("ny", 65 );
	//  pl->set("nz", 2 );
	pl->set("nz", 33 );
	
	//pl->set("npx", 1 );
	//pl->set("npy", 1 );
	//  pl->set("npz", 1 );
	
}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionSOp ) {


  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );
		
	isImpactInit=true;

  auto u =
      Teuchos::tuple(
          Teuchos::tuple(
              Pimpact::createScalarField( space, Pimpact::U ),
              Pimpact::createScalarField( space, Pimpact::U ),
              Pimpact::createScalarField( space, Pimpact::U )
          ),
          Teuchos::tuple(
              Pimpact::createScalarField( space, Pimpact::V ),
              Pimpact::createScalarField( space, Pimpact::V ),
              Pimpact::createScalarField( space, Pimpact::V )
          ),
          Teuchos::tuple(
              Pimpact::createScalarField( space, Pimpact::W ),
              Pimpact::createScalarField( space, Pimpact::W ),
              Pimpact::createScalarField( space, Pimpact::W )
          )
      );

  auto x = Pimpact::create<Pimpact::VectorField>( space );
  auto y = Pimpact::create<Pimpact::VectorField>( space );

  auto op = Pimpact::create<Pimpact::ConvectionSOp>( space ) ;

  op->print();

  for( int i=0; i<space->dim(); ++i )
    for( int j=0; j<space->dim(); ++j )
      u[i][j]->initField( Pimpact::ConstField, 2. );


  x->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
  x->getFieldPtr(1)->initField( Pimpact::Grad2D_inY );

  if( space->rankST()==0 )
    x->print();

  //x->write(0);
  y->random();

  //  op->apply( u[0], x->getConstField(0), y->getField(0), 1. );
  for( int i=0; i<space->dim(); ++i ) {
    op->apply( u[i], x->getConstField(i), y->getField(i) );
  }

  for( int i=0; i<space->dim(); ++i ) {
    TEST_FLOATING_EQUALITY( 2., y->getFieldPtr(i)->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( 2.* (double) y->getFieldPtr(i)->getLength()  , y->getFieldPtr(i)->norm(Belos::OneNorm), eps );
    TEST_FLOATING_EQUALITY( std::sqrt( 4.* y->getFieldPtr(i)->getLength() ), y->getFieldPtr(i)->norm(Belos::TwoNorm), eps );
  }

  //y->write(1);

}




template<class T> using ConvOpT = Pimpact::ConvectionVOp<Pimpact::ConvectionSOp<T> >;


TEUCHOS_UNIT_TEST( BasicOperator, ConvectionVOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

	isImpactInit=true;

  auto x = Pimpact::create<Pimpact::VectorField>( space );
  auto y = Pimpact::create<Pimpact::VectorField>( space );
  auto z = Pimpact::create<Pimpact::VectorField>( space );
  auto z2 = Pimpact::create<Pimpact::VectorField>( space );


  auto op = Pimpact::create<ConvOpT>( space  ) ;

  auto op2 = Pimpact::createConvectionOp( space ) ;


  x->initField( Pimpact::ConstFlow, 2., 2., 2. );
  y->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
  y->getFieldPtr(1)->initField( Pimpact::Grad2D_inY );
  z->random();

  op->assignField( *x );
  op->apply( *y, *z );
  op2->apply( *x, *y, *z2 );

  for( int i=0; i<space->dim(); ++i ) {
    TEST_FLOATING_EQUALITY( 2., z->getFieldPtr(i)->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( 2.* (double) z->getFieldPtr(i)->getLength()  , z->getFieldPtr(i)->norm(Belos::OneNorm), eps );
    TEST_FLOATING_EQUALITY( std::sqrt( 4.* z->getFieldPtr(i)->getLength() ), z->getFieldPtr(i)->norm(Belos::TwoNorm), eps );
  }
  //z->write(2);
  //z2->write(3);
  z2->add( -1, *z2, 1, *z );

  TEST_FLOATING_EQUALITY( z2->norm(), 0., eps );

  //z2->write(1);

}



TEUCHOS_UNIT_TEST( BasicOperator, ConvectionDiffusionOp  ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

	isImpactInit=true;

  auto x = Pimpact::create<Pimpact::VectorField>( space );
  auto y = Pimpact::create<Pimpact::VectorField>( space );
  auto z = Pimpact::create<Pimpact::VectorField>( space );
  auto z2 = Pimpact::create<Pimpact::VectorField>( space );
  auto sol = y->clone(Pimpact::ShallowCopy);


  auto op = Pimpact::createAdd2Op(
      Pimpact::create< ConvOpT<SpaceT> >(space),
      Pimpact::create<Pimpact::HelmholtzOp>(space) );

  auto op2 =
      Pimpact::create< ConvDiffOpT<SpaceT> >( space );

  op2->print();


  x->initField( Pimpact::ConstFlow, 2., 2., 2. );
  sol->init( Teuchos::tuple(2.,2.,2.) );

  // test in flow dir
  y->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
  y->getFieldPtr(1)->initField( Pimpact::Grad2D_inY );

  op->assignField( *x );
  op->apply( *y, *z );

  // connsistent test
  op2->assignField( *x );
  op2->apply(  *y, *z2 );

  //z->write(2);
  //z2->write(3);
  //sol->write(5);
  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

  TEST_EQUALITY( z2->norm() < eps, true );

  TEST_EQUALITY( sol->norm() < eps, true );

  //z2->write(1);
  //sol->write(4);

  // test against flow dir
  y->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
  y->getFieldPtr(1)->initField( Pimpact::Grad2D_inX );
  sol->init( Teuchos::tuple(2.,2.,2.) );
  z->random();

  op->assignField( *x );
  op->apply( *y, *z );

  // connsistent test
  op2->assignField( *x );
  op2->apply(  *y, *z2 );

  //z->write(2);
  //z2->write(3);
  //sol->write(5);
  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

  TEST_EQUALITY( z2->norm()<eps, true );

  TEST_EQUALITY( sol->norm()<eps, true );

  //z2->write(1);
  //sol->write(4);

  // test against flow dir
  y->initField( Pimpact::PoiseuilleFlow2D_inX );
  sol->init( Teuchos::tuple(0.,0.,0.) );
  z->random();

  op->apply( *y, *z );

  // connsistent test
  op2->assignField( *x );
  op2->apply(  *y, *z2 );

  //z->write(2);
  //z2->write(3);
  //sol->write(5);
  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

  TEST_EQUALITY( z2->norm()<eps, true );

  TEST_EQUALITY( sol->getConstFieldPtr(1)->norm()<eps, true );

  //z2->write(1);
  //sol->write(4);

  // test against flow dir
  y->initField( Pimpact::PoiseuilleFlow2D_inY );
  sol->init( Teuchos::tuple(0.,0.,0.) );
  z->random();

  op->assignField( *x );
  op->apply( *y, *z );

  // connsistent test
  op2->assignField( *x );
  op2->apply(  *y, *z2 );

  //z->write(2);
  //z2->write(3);
  //sol->write(5);
  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

  TEST_EQUALITY( z2->norm()<eps, true );

//  TEST_EQUALITY( sol->getConstFieldPtr(0)->norm()<eps, true );

  //z2->write(1);
  //sol->write(4);

}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionDiffusionSORSmoother ) {

  pl->set<S>("Re",1000);

  auto space = Pimpact::createSpace<S,O,d,2>( pl, !isImpactInit );

	isImpactInit=true;

  auto wind = Pimpact::create<Pimpact::VectorField>( space );
  auto y = Pimpact::create<Pimpact::VectorField>( space );
  auto z = Pimpact::create<Pimpact::VectorField>( space );
  auto z2 = Pimpact::create<Pimpact::VectorField>( space );


  auto op = Pimpact::create<ConvDiffOpT>( space );

  auto pls = Teuchos::parameterList();
  pls->set( "omega", 1. );
  pls->set( "numIter", 1 );
  pls->set( "Ordering", 0 );
  pls->set<short int>( "dir X",  1 );
  pls->set<short int>( "dir Y",  1 );
  pls->set<short int>( "dir Z",  1 );
  //  pls->set( "numIters",10)

  auto smoother = Pimpact::create< ConvDiffSORT >( op, pls );

  smoother->print();

  wind->initField( Pimpact::ConstFlow, winds, winds, winds );

  z->initField( Pimpact::ConstFlow, 0., 0., 0. );

  y->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
  y->getFieldPtr(1)->initField( Pimpact::Grad2D_inX );

  auto sol = y->clone( Pimpact::DeepCopy );
  //  sol->write(2);

  op->assignField( *wind );

  op->apply( *y, *z );

  {
    y->init(0);
    auto bc = z->clone( Pimpact::ShallowCopy );
    op->apply( *y, *bc );
    z->add( 1., *z, -1., *bc );
  }
  z->write(7777);


  y->initField( Pimpact::ConstFlow, 0., 0., 0. );
  //  y->init(0);
  //  y->random();

  std::ofstream ofs;
  if( space()->rankST()==0 )
    ofs.open("GS.txt", std::ofstream::out);

  double n;

  for(int i=0; i<20; ++i) {
    y->write(10+i);
    smoother->apply( *z, *y );

    z2->add( -1, *sol, 1, *y );
    //z2->write(100+i);

    n = z2->norm();
    if( space()->rankST()==0 )
      std::cout << "error: " << n << "\n";

    if( space()->rankST()==0 )
      ofs << n << "\n";

  }

  TEST_EQUALITY_CONST( z2->norm()/z2->getLength()<eps, true );

  if( space()->rankST()==0 )
    ofs.close();


}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionDiffusionJSmoother ) {

  pl->set<S>("Re",1);
  auto space = Pimpact::createSpace<S,O,d,2>( pl, !isImpactInit );

	isImpactInit=true;

  auto wind = Pimpact::create<Pimpact::VectorField>( space );
  auto y = Pimpact::create<Pimpact::VectorField>( space );
  auto z = Pimpact::create<Pimpact::VectorField>( space );
  auto z2 = Pimpact::create<Pimpact::VectorField>( space );

  auto op = Pimpact::create<ConvDiffOpT>( space );


  auto pls = Teuchos::parameterList();
  pls->set( "omega", 0.6 );
  pls->set( "numIters", 10 );

  auto smoother =
      Pimpact::create<
        Pimpact::ConvectionVSmoother<
          ConvDiffOpT<Pimpact::Space<S,O,d,2> > ,
          Pimpact::ConvectionDiffusionJSmoother > > (
              op,
              pls );

  wind->initField( Pimpact::ConstFlow, winds, winds, winds );
  z->initField( Pimpact::ConstFlow, 0., 0., 0. );

  y->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
  y->getFieldPtr(1)->initField( Pimpact::Grad2D_inX );

  auto sol = y->clone( Pimpact::DeepCopy );
  //  sol->write(2);

  op->assignField( *wind );

  op->apply( *y, *z );

  {
    y->init(0);
    auto bc = z->clone( Pimpact::ShallowCopy );
    op->apply( *y, *bc );
    z->add( 1., *z, -1., *bc );
  }

  //z->write(3);

  y->initField( Pimpact::ConstFlow, 0., 0., 0. );
  //  y->init(0);
  //  y->random();

  std::ofstream ofs;
  if( space()->rankST()==0 )
    ofs.open("GS.txt", std::ofstream::out);

  //  ofs << "lorem ipsum";

//  smoother->assignField( *wind );

  for(int i=0; i<200; ++i) {
    //y->write(10+i);
    smoother->apply( *z, *y );

    z2->add( -1, *sol, 1, *y );
    //z2->write(100+i);

    double n = z2->norm();
    if( space()->rankST()==0 )
      std::cout << "error: " << n << "\n";

    if( space()->rankST()==0 )
      ofs << n << "\n";

  }

  TEST_EQUALITY_CONST( z2->norm()/z2->getLength()<eps, true );

  if( space()->rankST()==0 )
    ofs.close();

}




TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicConvectionOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

	isImpactInit=true;

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto mv1 = Pimpact::createMultiHarmonicVectorField( space, 10 );
  auto mv2 = Pimpact::createMultiHarmonicVectorField( space, 10 );

  auto op = Pimpact::createMultiHarmonicConvectionOp( space, 10 );

  op->assignField( *mv1 );
  op->apply( *mv1, *mv2 );

}



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicDtConvectionDiffusionOp ) {

  pl->set<S>( "Re", 1.e1 );
  pl->set<S>( "alpha2", 20. );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

	S re = space->getDomain()->getDomainSize()->getRe();
	S alpha2 = space->getDomain()->getDomainSize()->getAlpha2();
	O nf = 4;

  isImpactInit=true;

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto wind = Pimpact::createMultiHarmonicVectorField( space, nf );
  auto x = Pimpact::createMultiHarmonicVectorField( space, nf );
  auto y1 = Pimpact::createMultiHarmonicVectorField( space, nf );
  auto y2 = Pimpact::createMultiHarmonicVectorField( space, nf );
  auto diff = Pimpact::createMultiHarmonicVectorField( space, nf );


	auto op1 =
		Pimpact::createAdd2Op(
				Pimpact::createMultiDtHelmholtz( space, alpha2/re, 1./re ),
				Pimpact::createMultiHarmonicConvectionOp( space, nf ) );

  auto op2 =
		Pimpact::createMultiDtConvectionDiffusionOp( space, nf );

  // init zero wind 
	wind->get0FieldPtr()->initField( Pimpact::ConstFlow, 1, 1, 0 );
	for( int i=0; i<nf; ++i ) {
		wind->getCFieldPtr(i)->initField( Pimpact::ConstFlow, 1, 1, 1 );
		wind->getSFieldPtr(i)->initField( Pimpact::ConstFlow, 1, 1, 1 );
	}
	// init const x
	//x->get0FieldPtr()->initField( Pimpact::ConstFlow, 1, 1, 1 );
	x->get0FieldPtr()->initField( Pimpact::PoiseuilleFlow2D_inX );
	for( int i=0; i<nf; ++i ) {
		x->getCFieldPtr(i)->initField( Pimpact::PoiseuilleFlow2D_inX );
		x->getSFieldPtr(i)->initField( Pimpact::PoiseuilleFlow2D_inX );
	}

	op1->assignField( *wind );
	op2->assignField( *wind );

  op1->apply( *x, *y1 );
  op2->apply( *x, *y2 );

	diff->add( 1., *y1, -1.,  *y2 );
	y1->write( 10 );
	y2->write( 20 );
	diff->write( 0 );

	std::cout << "error: " << diff->norm()/std::sqrt( diff->getLength() ) << "\n";
  TEST_EQUALITY_CONST( diff->norm()/std::sqrt( diff->getLength() ) <eps, true );

}





} // namespace
