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


bool testMpi = true;
S eps = 1e-10;
S omega = 0.8;
int nIter = 1000;

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
      "nIter", &nIter,
      "Slack off of machine epsilon used to check test results" );
  pl->set( "domain", 1);

  pl->set( "lz", 1. );


//  pl->set("nx", 125 );
//  pl->set("ny", 125 );
//  pl->set("nx", 65 );
//  pl->set("ny", 65 );
//  pl->set("nz", 2 );
//
//  pl->set("npx", 1 );
//  pl->set("npy", 1 );
//  pl->set("npz", 1 );

  // integer coefficients:
//  pl->set("nx", 9 );
//  pl->set("ny", 9 );
//  pl->set( "lx", 8. );
//  pl->set( "ly", 8. );

}



TEUCHOS_UNIT_TEST( BasicOperator, DivOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p   = Pimpact::createScalarField( space );
  auto vel = Pimpact::createVectorField( space );

  // zero test
  vel->initField();
  p->init(2.);

  auto op = Pimpact::createDivOp( space );

  op->print();

  op->apply(*vel,*p);

  p->write(0);

  TEST_EQUALITY( p->norm()<eps, true );


  // random test
  vel->random();
  p->init(2.);

  op->apply(*vel,*p);

  p->write(1);

  TEST_EQUALITY( p->norm()>eps, true );

  // circle test
  vel->initField( Pimpact::Circle2D );
  p->init(2.);

  op->apply(*vel,*p);

  TEST_EQUALITY( p->norm()<eps, true );

  p->write(2);

}



TEUCHOS_UNIT_TEST( BasicOperator, InterpolateV2SOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::createScalarField( space );
  auto vel = Pimpact::createVectorField( space );

  auto op = Pimpact::createInterpolateV2S( space );

  // PoiseuilleFlow2D_inX test
  vel->initField( Pimpact::PoiseuilleFlow2D_inX );

  for( int i=0; i<space->dim(); ++i ) {
    p->random();
    op->apply( vel->getConstField( i ), *p );
    p->write(i);
  }

  // PoiseuilleFlow2D_inY test
  vel->initField( Pimpact::PoiseuilleFlow2D_inY );

  for( int i=0; i<space->dim(); ++i ) {
    p->random();
    op->apply( vel->getConstField( i ), *p );
    p->write( i+space->dim() );
  }

  // PoiseuilleFlow2D_inY test
  vel->init( 1. );

  for( int i=0; i<space->dim(); ++i ) {
    p->random();
    op->apply( vel->getConstField( i ), *p );
    p->write( i+space->dim()*2 );
  }
  //  TEST_EQUALITY( p->norm()<eps, true );


}



TEUCHOS_UNIT_TEST( BasicOperator, InterpolateS2VOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::createScalarField( space );
  auto vel = Pimpact::createVectorField( space );

  auto op = Pimpact::createInterpolateS2V( space );
  op->print();

  // zero test
  p->initField( Pimpact::Poiseuille2D_inX );

  for( int i=0; i<space->dim(); ++i ) {
    vel->getFieldPtr( i )->random();
    op->apply( *p, vel->getField( i ) );
  }
  vel->write(0);

  // zero test
  p->initField( Pimpact::Poiseuille2D_inY );

  for( int i=0; i<space->dim(); ++i ) {
    vel->getFieldPtr( i )->random();
    op->apply( *p, vel->getField( i ) );
  }
  vel->write(1);

  //  TEST_EQUALITY( p->norm()<eps, true );


}



TEUCHOS_UNIT_TEST( BasicOperator, TransferOp ) {

  auto fSpace = Pimpact::createSpace<S,O,d,4>( pl );
  auto cSpace = Pimpact::createSpace<S,O,d,2>( pl );

  auto fx = Pimpact::createScalarField( fSpace );
  auto cx = Pimpact::createScalarField( cSpace );

  auto op = Pimpact::createTransferOp( fSpace, cSpace );

  // test
  fx->initField( Pimpact::Poiseuille2D_inX );
  cx->random();

  op->apply( *fx, *cx );

  TEST_FLOATING_EQUALITY( fx->norm(Belos::OneNorm), cx->norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::TwoNorm), cx->norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::InfNorm), cx->norm(Belos::InfNorm), eps );

  cx->write(0);

  fx->random();

  op->apply( *cx, *fx );

  TEST_FLOATING_EQUALITY( fx->norm(Belos::OneNorm), cx->norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::TwoNorm), cx->norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx->norm(Belos::InfNorm), cx->norm(Belos::InfNorm), eps );

  fx->write(1);

}



TEUCHOS_UNIT_TEST( BasicOperator, GradOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p = Pimpact::createScalarField( space );
  auto v = Pimpact::createVectorField( space );

  auto op = Pimpact::createGradOp( space );
  op->print();

  // op in x test
  p->initField( Pimpact::Grad2D_inX );
  v->random();

  op->apply(*p,*v);

  TEST_EQUALITY( (v->getConstFieldPtr(Pimpact::U)->norm()-std::sqrt( std::pow(1.,2)*p->getLength() ))<eps, true );
  TEST_EQUALITY( v->getConstFieldPtr(Pimpact::V)->norm()<eps, true );

  v->write(0);

  // op in x test
  p->initField( Pimpact::Grad2D_inY );
  v->random();

  op->apply(*p,*v);

  TEST_EQUALITY(  v->getConstFieldPtr(Pimpact::U)->norm()<eps, true );
  TEST_EQUALITY( (v->getConstFieldPtr(Pimpact::V)->norm()-std::sqrt( std::pow(1.,2)*p->getLength() ))<eps, true );

  v->write(1);

  TEST_EQUALITY( 0, 0 );
}



TEUCHOS_UNIT_TEST( BasicOperator, HelmholtzOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto x = Pimpact::createVectorField(space);
  auto bs= Pimpact::createVectorField(space);
  auto b = Pimpact::createVectorField(space);

  double mulI = 5.;
  double mulL = 3.;

  auto op= Pimpact::createHelmholtzOp( space, mulI, mulL );

  op->print();

  // test in x direction
  x->initField( Pimpact::PoiseuilleFlow2D_inX );
  bs->init( Teuchos::tuple( 8./std::pow(space->getDomain()->getDomainSize()->getSize(Pimpact::Y),2), 0., 0. ) );
  bs->add( mulI, *x, mulL, *bs );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);
  b->write();
  bs->write(1);

  std::cout << "error: " << bs->norm() << "\n";
  TEST_EQUALITY( bs->norm()<eps, true );

  // test in y direction
  x->initField( Pimpact::PoiseuilleFlow2D_inY );
  bs->init( Teuchos::tuple( 0., 8./std::pow(space->getDomain()->getDomainSize()->getSize(Pimpact::X),2), 0. ) );
  bs->add( mulI, *x, mulL, *bs );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);
  b->write(2);
  bs->write(3);

  std::cout << "error: " << bs->norm() << "\n";
  TEST_EQUALITY( bs->norm()<eps, true );

  // the circle test
  x->initField( Pimpact::Circle2D );
  bs->initField( Pimpact::Circle2D );
  bs->scale( mulI );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b );
  x->write(4);
  b->write(5);
  bs->write(6);

  std::cout << "error: " << bs->norm() << "\n";
  TEST_EQUALITY( bs->norm()<eps, true );

}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionSOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing


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

  auto x = Pimpact::createVectorField( space );
  auto y = Pimpact::createVectorField( space );

  auto op = Pimpact::createConvectionSOp( space ) ;

  op->print();

  for( int i=0; i<space->dim(); ++i )
    for( int j=0; j<space->dim(); ++j )
      u[i][j]->init( 2. );


  x->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
  x->getFieldPtr(1)->initField( Pimpact::Grad2D_inY );
  //  x->scale();

  x->write(0);
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

  y->write(1);

}



TEUCHOS_UNIT_TEST( BasicOperator, ConvectionVOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing


  auto x = Pimpact::createVectorField( space );
  auto y = Pimpact::createVectorField( space );
  auto z = Pimpact::createVectorField( space );
  auto z2 = Pimpact::createVectorField( space );

  auto op = Pimpact::createConvectionVOp( space ) ;
  auto op2 = Pimpact::createConvectionOp( space ) ;


  for(int i=0; i<10; ++i )
    x->write(i);

  //  x->write(0);
  x->initField( Pimpact::ConstFlow, 2., 2., 2. );
  y->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
  y->getFieldPtr(1)->initField( Pimpact::Grad2D_inY );
  z->random();

  op->apply( *x, *y, *z );
  op2->apply( *x, *y, *z2 );

  for( int i=0; i<space->dim(); ++i ) {
    TEST_FLOATING_EQUALITY( 2., z->getFieldPtr(i)->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( 2.* (double) z->getFieldPtr(i)->getLength()  , z->getFieldPtr(i)->norm(Belos::OneNorm), eps );
    TEST_FLOATING_EQUALITY( std::sqrt( 4.* z->getFieldPtr(i)->getLength() ), z->getFieldPtr(i)->norm(Belos::TwoNorm), eps );
  }
  z2->add( -1, *z2, 1, *z );

  TEST_FLOATING_EQUALITY( z2->norm(), 0., eps );

  z2->write(1);

}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionJacobianVOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing


  auto x = Pimpact::createVectorField( space );
  auto y = Pimpact::createVectorField( space );
  auto z = Pimpact::createVectorField( space );
  auto z2 = Pimpact::createVectorField( space );

  auto op = Pimpact::createConvectionJacobianOp( space, false ) ;
  auto op2 = Pimpact::createConvectionOp( space ) ;

  x->initField( Pimpact::ConstFlow, 2., 2., 2. );

  y->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
  y->getFieldPtr(1)->initField( Pimpact::Grad2D_inY );

  z->random();

  op->assignField( *x );
  op->apply( *y, *z );

  op2->apply( *x, *y, *z2 );

  // test is not exact because the boundaries of x are zero, and
  for( int i=0; i<space->dim(); ++i ) {
    TEST_FLOATING_EQUALITY( 2., z->getFieldPtr(i)->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( 2.* (double) z->getFieldPtr(i)->getLength()  , z->getFieldPtr(i)->norm(Belos::OneNorm), eps );
    TEST_FLOATING_EQUALITY( std::sqrt( 4.* z->getFieldPtr(i)->getLength() ), z->getFieldPtr(i)->norm(Belos::TwoNorm), eps );
  }
  z2->add( -1, *z2, 1, *z );

  TEST_FLOATING_EQUALITY( z2->norm(), 0., eps );

  z2->write(2);
  z->write(1);

}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef Pimpact::VectorField<SpaceT> VF;
  typedef Pimpact::MultiField<VF> MVF;


  auto vel = Pimpact::createVectorField( space );

  auto x = Pimpact::createMultiField(*vel->clone(),10);
  auto y = Pimpact::createMultiField(*vel->clone(),10);

  auto op = Pimpact::createMultiOperatorBase<MVF>( Pimpact::createConvectionOp( space ) );

  for( O i=0; i<10; ++i ) {
    x->getFieldPtr(i)->initField(Pimpact::Circle2D );
  }
  //  x->random();
  x->getFieldPtr(0)->write();

  op->apply( *x, *y);

  y->getFieldPtr(0)->write(99);
}


TEUCHOS_UNIT_TEST( BasicOperator, DivGradO2Op ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto x = Pimpact::createScalarField( space );
  auto b = Pimpact::createScalarField( space );

  // zero test
  b->initField();
  x->init(2.);

  auto op = Pimpact::createDivGradO2Op( space );

  op->print();

  op->apply(*b,*x);

  x->write(0);

  TEST_EQUALITY( x->norm()<eps, true );


  // random test
  b->random();
  x->init(2.);

  op->apply(*b,*x);

  x->write(1);

  TEST_EQUALITY( x->norm()>eps, true );

  // circle test
  b->initField( Pimpact::Grad2D_inX );
  x->init(2.);

  op->apply(*b,*x);

  TEST_EQUALITY( x->norm()<eps, true );

  x->write(2);

}



TEUCHOS_UNIT_TEST( BasicOperator, DivGradO2JSmoother ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto x = Pimpact::createScalarField( space );
  auto b = Pimpact::createScalarField( space );

  auto op = Pimpact::createDivGradO2Op( space );
  auto smoother = Pimpact::createDivGradO2JSmoother<SpaceT>( op, omega, nIter );

  smoother->print();

  // Grad in x
  x->initField( Pimpact::Grad2D_inX );
  x->write(0);

  op->apply(*x,*b);
  b->write(1);

  x->init();

  smoother->apply( *b, *x );

  x->write(2);

  b->initField( Pimpact::Grad2D_inX );
  x->add( -1, *x, 1., *b );
  std::cout << "error: " << x->norm()/b->norm() << "\n";
  TEST_EQUALITY( x->norm()/b->norm()<0.5, true );

  // Grad in y
  x->initField( Pimpact::Grad2D_inY );
  x->write(3);

  op->apply(*x,*b);
  b->write(4);

  x->init();

  smoother->apply( *b, *x );

  x->write(5);

  b->initField( Pimpact::Grad2D_inY );
  x->add( -1, *x, 1., *b );
  std::cout << "error: " << x->norm()/b->norm() << "\n";
  TEST_EQUALITY( x->norm()/b->norm()<0.5, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, ForcingOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto vel = Pimpact::createVectorField(space);
  auto force = Pimpact::createVectorField(space);
  auto res = Pimpact::createVectorField(space);

  vel->init(1.);
  force->initField( Pimpact::EFlowField(11) );
  res->random();

  auto op = Pimpact::createForcingOp( force, 1. );

  op->apply(*vel,*res);
  vel->add( 1., *res, -1., *force );

  TEST_EQUALITY( vel->norm(), 0 );
}




TEUCHOS_UNIT_TEST( BasicOperator, Add2Op ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;


  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing


  auto vel = Pimpact::createVectorField(space);


  auto x = Pimpact::createMultiField(*vel->clone(),10);
  auto y = Pimpact::createMultiField(*vel->clone(),10);

  auto op = Pimpact::createOperatorBase(
      Pimpact::createMultiOpWrap(
          Pimpact::createAdd2Op(
              Pimpact::createHelmholtzOp( space ),
              Pimpact::createConvectionOp( space ),
              vel->clone() ) ) );

  for( O i=0; i<10; ++i ) {
    x->getFieldPtr(i)->initField(Pimpact::RankineVortex2D );
  }

  x->getFieldPtr(0)->write();

  op->apply( *x, *y);

  y->getFieldPtr(0)->write(99);

}




TEUCHOS_UNIT_TEST( ModeOperator, HelmholtzOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;


  auto velc = Pimpact::createVectorField(space);
  auto vels = Pimpact::createVectorField(space);

  auto vel = Pimpact::createModeField( velc, vels );

  auto mv = Pimpact::createMultiField(*vel,10);

  auto mv2 = mv->clone(11);

  auto op = Pimpact::createMultiOpWrap(
      Pimpact::createModeOpWrap(
          Pimpact::createHelmholtzOp(space) ) );

  op->apply(*mv,*mv2);

}




TEUCHOS_UNIT_TEST( ModeOperator, DtModeOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  typedef Pimpact::MultiField< Pimpact::ModeField< Pimpact::VectorField<SpaceT> > > MF;
  typedef Pimpact::OperatorBase<MF> BOp;


  auto velc = Pimpact::createVectorField(space);
  auto vels = Pimpact::createVectorField(space);

  auto vel = Pimpact::createModeField( velc, vels );

  auto mv = Pimpact::createMultiField(*vel,1);

  mv->getFieldPtr(0)->getCFieldPtr()->init(1.);
  mv->getFieldPtr(0)->getSFieldPtr()->init(0.);

  auto mv2 = mv->clone(1);
  mv2->init(0.);

  auto op = Pimpact::createMultiOperatorBase<MF>(
      Pimpact::createDtModeOp<SpaceT>() );

  op->apply(*mv,*mv2);

  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm() );
  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm() );
  Belos::OperatorTraits<S,MF,BOp>::Apply(*op,*mv,*mv2);

}




TEUCHOS_UNIT_TEST( ModeOperator, DtLapOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto velc = Pimpact::createVectorField(space);
  auto vels = Pimpact::createVectorField(space);

  auto vel = Pimpact::createModeField( velc, vels );

  auto mv = Pimpact::createMultiField(*vel,1);

  mv->random();

  auto mv2 = mv->clone(1);
  auto mv3 = mv->clone(1);
  mv2->init(0.);
  mv3->init(0.);

  auto A = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp( space, 1., 0. ) );

  A->apply(*mv,*mv2);

  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm() );
  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm() );
  //	TEST_EQUALITY( mv->getConstField(0).getConstFieldC()->norm(), mv2->getConstField(0).getConstFieldS()->norm() );
  //	TEST_EQUALITY( mv->getConstField(0).getConstFieldS()->norm(), mv2->getConstField(0).getConstFieldC()->norm() );

  mv2->init(0.);

  auto A2 = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp( space, 0., 1. ) );
  auto A3 = Pimpact::createMultiOpWrap( Pimpact::createModeOpWrap( Pimpact::createHelmholtzOp( space, 0., 1. ) ) );

  A2->apply(*mv,*mv2);
  A3->apply(*mv,*mv3);

  std::vector<S> norm2( 1, 0. );
  std::vector<S> norm3( 1, 0. );

  mv2->norm(norm2);
  mv3->norm(norm3);

  TEST_EQUALITY( norm2[0] , norm3[0] );

}




TEUCHOS_UNIT_TEST( ModeOperator, TripleCompostion ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;


  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > > MVF;


  auto X = Pimpact::createMultiModeScalarField( space );
  auto B = Pimpact::createMultiModeScalarField( space );

  auto temp = Pimpact::createMultiModeVectorField( space );

  X->init(0.);
  B->random();


  auto H =
      Pimpact::createMultiOperatorBase<MVF>(
          Pimpact::createDtLapOp( space, 0., 10. ) );

  H->apply( *temp, *temp );

  // Make an empty new parameter list.
  auto solverParams = Teuchos::parameterList();

  // Create the Pimpact::LinearSolver solver.
  auto Hprob = Pimpact::createLinearProblem<MVF>( H, temp, temp, solverParams,"GMRES" );
  auto Hinv  = Pimpact::createInverseOperator( Hprob );

  auto schur = Pimpact::createTripleCompositionOp(
      temp->clone(),
      temp->clone(),
      Pimpact::createMultiModeOpWrap( Pimpact::createDivOp( space ) ),
      Hinv,
      Pimpact::createMultiModeOpWrap( Pimpact::createGradOp( space ))
  );

  schur->apply( *B, *X );

}




TEUCHOS_UNIT_TEST( ModeOperator, InverseOperator ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > > BVF;


  auto X = Pimpact::createMultiModeVectorField( space );
  auto B = Pimpact::createMultiModeVectorField( space );

  X->init(0.);
  B->random();

  auto op = Pimpact::createMultiModeOperatorBase<BVF>(
      Pimpact::createHelmholtzOp(space) );

  // Make an empty new parameter list.
  auto solverParams = Pimpact::createLinSolverParameter( "CG", 1.e-1 );

  // Create the Pimpact::LinearSolver solver.
  auto prob =
      Pimpact::createLinearProblem<BVF>(
          op,
          X,X,
          //          Pimpact::createMultiField(X->getFieldPtr(0)->getCFieldPtr()),
          //          Pimpact::createMultiField(B->getFieldPtr(0)->getCFieldPtr()),
          solverParams );

  auto opinv = Pimpact::createInverseOperator<BVF>( prob );
  auto opp = Pimpact::createOperatorBase( opinv );

}





TEUCHOS_UNIT_TEST( ModeOperator, EddyPrec ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > > MVF;
  typedef Pimpact::HelmholtzOp<SpaceT>  Op;
  typedef Pimpact::EddyPrec<SpaceT>  Op2;

  auto temp = Pimpact::createMultiModeVectorField( space );

  auto X = Pimpact::createMultiModeVectorField( space );
  auto B = Pimpact::createMultiModeVectorField( space );

  X->init(0.);
  B->random();

  // Make an empty new parameter list.
  auto solverParams = Pimpact::createLinSolverParameter("CG",1.e-1);

  // Create the Pimpact::LinearSolver solver.
  auto A = Pimpact::createMultiModeOperatorBase<MVF,Op>( Pimpact::createHelmholtzOp( space, 14.,1.) );

  A->apply( *temp, *temp );

  auto prob = Pimpact::createLinearProblem<MVF>( A, temp, temp, solverParams,"CG" );

  prob->solve(temp,temp);

  auto op = Pimpact::createEddyPrec<SpaceT>( temp, Pimpact::createInverseOperatorBase(prob) ) ;

  op->apply( X->getField(0), B->getField(0) );

  auto schur = Pimpact::createMultiOperatorBase<MVF,Op2>( op );

  schur->apply( *B, *X );

}




TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicConvectionOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto vel = Pimpact::createVectorField( space );

  auto mv1 = Pimpact::createMultiHarmonicVectorField( space, 10 );
  auto mv2 = Pimpact::createMultiHarmonicVectorField( space, 10 );

  auto op = Pimpact::createMultiHarmonicConvectionOp( space );

  op->apply( *mv1, *mv2 );
}



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicOpWrap ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto vel = Pimpact::createVectorField( space );

  auto mv1 = Pimpact::createMultiHarmonicVectorField( space, 10 );
  auto mv2 = Pimpact::createMultiHarmonicVectorField( space, 10 );

  auto op = Pimpact::createMultiHarmonicOpWrap< Pimpact::HelmholtzOp<SpaceT> >(
      Pimpact::createHelmholtzOp(space) );

  op->apply( *mv1, *mv2 );
}



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiDtHelmholtz) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto mv1 = Pimpact::createMultiHarmonicVectorField( space, 10 );
  auto mv2 = Pimpact::createMultiHarmonicVectorField( space, 10 );

  auto op = Pimpact::createMultiDtHelmholtz( space );

  op->apply( *mv1, *mv2 );
}



TEUCHOS_UNIT_TEST( CompoundOperator, CompoundOpWrap ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  typedef Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> > VF;
  typedef Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> > SF;
  typedef Pimpact::CompoundField< VF, SF> CF;
  typedef Pimpact::MultiField<CF> MF;


  int nfs = 4;

  auto x    = Pimpact::createMultiField( Pimpact::createCompoundField(
      Pimpact::createMultiHarmonicVectorField( space, nfs ),
      Pimpact::createMultiHarmonicScalarField( space, nfs )) );
  auto fu   = x->clone();
  x->init(1.);
  x->random();

  auto opV2V =
      Pimpact::createAdd2Op(
          Pimpact::createMultiDtHelmholtz( space, 1., 1. ),
          Pimpact::createMultiHarmonicConvectionOp( space ),
          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone()
      );

  auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::createGradOp( space ) );
  auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::createDivOp( space ) );

  auto op =
      Pimpact::createMultiOperatorBase<MF>(
          Pimpact::createCompoundOpWrap(
              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(), opV2V, opS2V, opV2S ) );

  // vector to vector operator
  fu->init(0.);
  TEST_EQUALITY( 0., fu->norm() );

  opV2V->apply( x->getConstFieldPtr(0)->getConstVField(), fu->getFieldPtr(0)->getVField() );

  TEST_INEQUALITY( 0., fu->norm() );


  // scalar to vector operator
  fu->init(0.);
  TEST_EQUALITY( 0., fu->norm() );

  opS2V->apply( x->getConstFieldPtr(0)->getConstSField(), fu->getFieldPtr(0)->getVField() );

  TEST_INEQUALITY( 0., fu->norm() );

  // vector to scalar operator
  fu->init(0.);
  TEST_EQUALITY( 0., fu->norm() );

  opV2S->apply( x->getConstFieldPtr(0)->getConstVField(), fu->getFieldPtr(0)->getSField() );

  TEST_INEQUALITY( 0., fu->norm() );

  // compound operator
  fu->init(0.);
  TEST_EQUALITY( 0., fu->norm() );
  op->apply(*x,*fu);
  TEST_INEQUALITY( 0., fu->norm() );

}


} // namespace
