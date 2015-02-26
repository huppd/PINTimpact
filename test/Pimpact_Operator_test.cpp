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

#include "Pimpact_VectorFieldOpWrap.hpp"

#include "Pimpact_LinSolverParameter.hpp"




namespace {


typedef double S;
typedef int O;
const int d = 3;
const int dNC = 4;
//const int dNC = 2;

typedef Pimpact::Space<S,O,d,dNC>             SpaceT;

typedef typename Pimpact::ScalarField<SpaceT> SF;
typedef typename Pimpact::VectorField<SpaceT> VF;
typedef typename Pimpact::ModeField<SF>       MSF;
typedef typename Pimpact::ModeField<VF>       MVF;



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

  pl->set( "domain", domain);

  pl->set( "dim", dim );

	pl->set("nx", 33 );
	pl->set("ny", 17 );
	pl->set("nz", 9 );
	pl->set("nf", 4 );
  //
  pl->set("npx", 2 );
  pl->set("npy", 2 );
  //  pl->set("npz", 1 );

}



TEUCHOS_UNIT_TEST( BasicOperator, DivOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto p   = Pimpact::create<Pimpact::ScalarField>( space );
  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  // zero test
  vel->initField();
  p->init(2.);

  auto op = Pimpact::create<Pimpact::DivOp>( space );

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

  auto p = Pimpact::create<Pimpact::ScalarField>( space );
  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto op = Pimpact::createInterpolateV2S( space );
  //  auto op = Pimpact::create<Pimpact::InterpolateV2S>( space );

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

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto p = Pimpact::create<Pimpact::ScalarField>( space );
  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto op = Pimpact::create<Pimpact::InterpolateS2V>( space );
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

  typedef Pimpact::Space<S,O,d,4> FSpaceT;
  typedef Pimpact::Space<S,O,d,2> CSpaceT;

  auto fSpace = Pimpact::createSpace<S,O,d,4>( pl );
  auto cSpace = Pimpact::createSpace<S,O,d,2>( pl );

  auto fx = Pimpact::create<Pimpact::ScalarField>( fSpace );
  auto cx = Pimpact::create<Pimpact::ScalarField>( cSpace );

  auto op = Pimpact::create< Pimpact::TransferOp<FSpaceT,CSpaceT> >( fSpace, cSpace );

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



TEUCHOS_UNIT_TEST( BasicOperator, VectorFieldOpWrap ) {
  typedef Pimpact::Space<S,O,d,4> FSpaceT;
  typedef Pimpact::Space<S,O,d,2> CSpaceT;

  auto fSpace = Pimpact::createSpace<S,O,d,4>( pl );
  auto cSpace = Pimpact::createSpace<S,O,d,2>( pl );

  auto fx = Pimpact::create<Pimpact::VectorField>( fSpace );
  auto cx = Pimpact::create<Pimpact::VectorField>( cSpace );

  auto op = Pimpact::create< Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<FSpaceT,CSpaceT> > >( fSpace, cSpace );

  // test
  fx->initField( Pimpact::PoiseuilleFlow2D_inX );
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

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto p = Pimpact::create<Pimpact::ScalarField>( space );
  auto v = Pimpact::create<Pimpact::VectorField>( space );

  auto op = Pimpact::create<Pimpact::GradOp>( space );
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


  S mulI = 5.;
  S mulL = 3.;

  pl->set( "alpha2", mulI );
  pl->set( "Re", 1./mulL );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto x = Pimpact::create<Pimpact::VectorField>(space);
  auto bs= Pimpact::create<Pimpact::VectorField>(space);
  auto b = Pimpact::create<Pimpact::VectorField>(space);

	S n = std::sqrt( (S)x->getLength() );

  auto op= Pimpact::create<Pimpact::HelmholtzOp>( space );

  op->print();

  // test in x direction
  x->initField( Pimpact::PoiseuilleFlow2D_inX );
  bs->init( Teuchos::tuple( 8./std::pow(space->getDomain()->getDomainSize()->getSize(Pimpact::Y),2), 0., 0. ) );
  bs->add( mulI, *x, mulL, *bs );

	auto para = Teuchos::parameterList();
	para->set<S>( "mulI", mulI );
	op->setParameter( para );
  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);
  b->write();
  bs->write(1);

  std::cout << "error: " << bs->norm()/n << "\n";
  TEST_EQUALITY( bs->norm()/n<eps, true );

  // test in y direction
  x->initField( Pimpact::PoiseuilleFlow2D_inY );
  bs->init( Teuchos::tuple( 0., 8./std::pow(space->getDomain()->getDomainSize()->getSize(Pimpact::X),2), 0. ) );
  bs->add( mulI, *x, mulL, *bs );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);
  b->write(2);
  bs->write(3);

  std::cout << "error: " << bs->norm()/n << "\n";
  TEST_EQUALITY( bs->norm()/n<eps, true );

  // the circle test
  x->initField( Pimpact::Circle2D );
  bs->initField( Pimpact::Circle2D );
  bs->scale( mulI );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b );
  x->write(4);
  b->write(5);
  bs->write(6);

  std::cout << "error: " << bs->norm()/n << "\n";
  TEST_EQUALITY( bs->norm()/n<eps, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, DivGradO2Op ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto x = Pimpact::create<Pimpact::ScalarField>( space );
  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  // zero test
  b->initField( Pimpact::ConstField, 0. );
  x->init(2.);

  auto op = Pimpact::createDivGradO2Op( space );

  op->print();

  op->apply(*b,*x);

  x->write(0);

	S n = std::sqrt( x->getLength() );
	std::cout << "error: " << x->norm()/n << "\n";
  TEST_EQUALITY( x->norm()/n<eps, true );


  // random test
  b->random();
  x->init(2.);

  op->apply(*b,*x);

  x->write(1);

	std::cout << "error: " << x->norm()/n << "\n";
  TEST_EQUALITY( x->norm()/n>eps, true );

  // circle test
  //b->initField( Pimpact::Grad2D_inX );
  b->initField( Pimpact::Grad2D_inX );
  x->init(2.);

  op->apply(*b,*x);

	//std::cout << "error: " << x->norm()/n << "\n";
  //TEST_EQUALITY( x->norm()/n<eps, true );

  x->write(2);

}



TEUCHOS_UNIT_TEST( BasicOperator, DivGradO2JSmoother ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto x = Pimpact::create<Pimpact::ScalarField>( space );
  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  auto op = Pimpact::createDivGradO2Op( space );

  auto ppl = Teuchos::parameterList();

  ppl->set("omega", 0.80 );
  ppl->set("numIters", 400);

  auto smoother = Pimpact::create< Pimpact::DivGradO2JSmoother >( op, ppl );

  smoother->print();

  // Grad in x
  x->initField( Pimpact::Grad2D_inX );
  //x->write(0);

  op->apply(*x,*b);
	
  //b->write(1);

  x->init();

  smoother->apply( *b, *x );

  //x->write(2);

  b->initField( Pimpact::Grad2D_inX );
  x->add( -1, *x, 1., *b );
  std::cout << "error: " << x->norm()/b->norm() << "\n";
  TEST_EQUALITY( x->norm()/b->norm()<1., true );

  // Grad in y
  x->initField( Pimpact::Grad2D_inY );
  //x->write(3);

  op->apply(*x,*b);
  //b->write(4);

  x->init();

  smoother->apply( *b, *x );

  //x->write(5);

  b->initField( Pimpact::Grad2D_inY );
  x->add( -1, *x, 1., *b );
  std::cout << "error: " << x->norm()/b->norm() << "\n";
  TEST_EQUALITY( x->norm()/b->norm()<1., true );

}



TEUCHOS_UNIT_TEST( BasicOperator, ForcingOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto vel = Pimpact::create<Pimpact::VectorField>(space);
  auto force = Pimpact::create<Pimpact::VectorField>(space);
  auto res = Pimpact::create<Pimpact::VectorField>(space);

  vel->init(1.);
  force->initField( Pimpact::EFlowField(11) );
  res->random();

  auto op = Pimpact::createForcingOp( force, 1. );

  op->apply(*vel,*res);
  vel->add( 1., *res, -1., *force );

  TEST_EQUALITY( vel->norm(), 0 );
}





template<class T> using MSFT = Pimpact::MultiField< Pimpact::ScalarField<T> >;
template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;

TEUCHOS_UNIT_TEST( MultiOperator, InverseOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto x = Pimpact::create<Pimpact::ScalarField>( space );
  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );
  //  auto opm = Pimpact::createMultiOpWrap( op );
  //  auto op = Pimpact::create<MOP>( space );

  //  auto opinv = Pimpact::create<Pimpact::InverseOp>( opm );
  auto opinv = Pimpact::create<MOP>( op );

  // Grad in x
  //  x->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
  x->initField( Pimpact::Grad2D_inX );
  x->write(0);

  op->apply(*x,*b);
  b->write(1);

  x->init();

  opinv->apply( *b, *x );

  x->write(2);

  //  b->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
  b->initField( Pimpact::Grad2D_inX );
  x->add( -1, *x, 1., *b );
  std::cout << "error: " << x->norm()/b->norm() << "\n";
  TEST_EQUALITY( x->norm()/b->norm()<0.5, true );

  // Grad in y
  //  x->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
  x->initField( Pimpact::Grad2D_inY );
  x->write(3);

  op->apply(*x,*b);
  b->write(4);

  x->init();

  opinv->apply( *b, *x );

  x->write(5);

  //  b->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
  b->initField( Pimpact::Grad2D_inY );
  x->add( -1, *x, 1., *b );
  std::cout << "error: " << x->norm()/b->norm() << "\n";
  TEST_EQUALITY( x->norm()/b->norm()<0.5, true );

}



TEUCHOS_UNIT_TEST( MultiOperator, Add2Op ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;


  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing


  auto velx = Pimpact::create<Pimpact::VectorField>(space);


  auto x = Pimpact::createMultiField(*velx->clone(),10);
  auto y = Pimpact::createMultiField(*velx->clone(),10);

  auto op = Pimpact::createOperatorBase(
      Pimpact::createMultiOpWrap(
          Pimpact::createAdd2Op(
              Pimpact::create<Pimpact::HelmholtzOp>( space ),
              Pimpact::createConvectionOp( space ) ) ) );


  for( O i=0; i<10; ++i ) {
    x->getFieldPtr(i)->initField(Pimpact::RankineVortex2D );
  }

  x->getFieldPtr(0)->write();

  op->apply( *x, *y);

  y->getFieldPtr(0)->write(99);

}


template<class T> using WUP = Pimpact::MultiOpUnWrap<Pimpact::MultiOpWrap<T> >;


TEUCHOS_UNIT_TEST( MultiOperator, MulitOpUnWrap ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto x = Pimpact::create<Pimpact::VectorField>(space);
  auto y = Pimpact::create<Pimpact::VectorField>(space);

  auto op = Pimpact::createOperatorBase(
      Pimpact::createMultiOpWrap(
          Pimpact::createAdd2Op(
              Pimpact::create<Pimpact::HelmholtzOp>( space ),
              Pimpact::createConvectionOp( space ) ) ) );

  auto opUW = Pimpact::create<Pimpact::MultiOpUnWrap>( op );

  auto opUW2 =
      Pimpact::create<WUP>(
          //      Pimpact::create<Pimpact::MultiOpWrap >(
          Pimpact::create<Pimpact::HelmholtzOp>(space)
          //          ;
      );


  x->initField(Pimpact::RankineVortex2D );

  x->write();

  opUW2->apply( *x, *y);

  y->write(99);

}



TEUCHOS_UNIT_TEST( MultiModeOperator, HelmholtzOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;


  auto velc = Pimpact::create<Pimpact::VectorField>(space);
  auto vels = Pimpact::create<Pimpact::VectorField>(space);

  auto vel = Pimpact::create<MVF>( space );

  auto mv = Pimpact::createMultiField(*vel,10);

  auto mv2 = mv->clone(10);

  auto op = Pimpact::createMultiOpWrap(
      Pimpact::createModeOpWrap(
          Pimpact::create<Pimpact::HelmholtzOp>(space) ) );

  op->apply(*mv,*mv2);

}




TEUCHOS_UNIT_TEST( MultiModeOperator, DtModeOp ) {

	pl->set<S>("Re", 1., "Reynolds number");
	pl->set<S>("alpha2", 1.,
			"Womersley square \alpha^2");
  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  typedef Pimpact::MultiField< Pimpact::ModeField< Pimpact::VectorField<SpaceT> > > MF;
  typedef Pimpact::OperatorBase<MF> BOp;


  auto velc = Pimpact::create<Pimpact::VectorField>(space);
  auto vels = Pimpact::create<Pimpact::VectorField>(space);

  auto vel = Pimpact::create<MVF>( space );

  auto mv = Pimpact::createMultiField(*vel,1);

  mv->getFieldPtr(0)->getCFieldPtr()->init(1.);
  mv->getFieldPtr(0)->getSFieldPtr()->init(0.);

  auto mv2 = mv->clone(1);
  mv2->init(0.);

  auto op = Pimpact::createMultiOperatorBase(
      Pimpact::create<Pimpact::DtModeOp>(space) );

  op->apply(*mv,*mv2);

  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm() );
  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm() );
  Belos::OperatorTraits<S,MF,BOp>::Apply(*op,*mv,*mv2);

}




TEUCHOS_UNIT_TEST( MultiModeOperator, DtLapOp ) {

	pl->set( "Re", 1. );
  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  //if( !isImpactInit ) isImpactInit=true;

  auto velc = Pimpact::create<Pimpact::VectorField>(space);
  auto vels = Pimpact::create<Pimpact::VectorField>(space);

  auto vel = Pimpact::create<MVF>( space );

  auto mv = Pimpact::createMultiField(*vel,1);

  mv->random();

  auto mv2 = mv->clone(1);
  auto mv3 = mv->clone(1);
  mv2->init(0.);
  mv3->init(0.);

  auto A = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp( space, 1., 0. ) );

	mv->random();
  A->apply(*mv,*mv2);

  TEST_FLOATING_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), eps );
  TEST_FLOATING_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), eps );

  mv2->init(0.);

  auto A2 = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp( space, 0., 1. ) );
  auto A3 = Pimpact::createMultiOpWrap( Pimpact::createModeOpWrap( Pimpact::create<Pimpact::HelmholtzOp>( space ) ) );

  A2->apply(*mv,*mv2);
  A3->apply(*mv,*mv3);

  std::vector<S> norm2( 1, 0. );
  std::vector<S> norm3( 1, 0. );

  mv2->norm(norm2);
  mv3->norm(norm3);

  TEST_EQUALITY( norm2[0] , norm3[0] );

}




TEUCHOS_UNIT_TEST( MultiModeOperator, TripleCompostion ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;


  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > > MVF;
  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<SpaceT> > > MSF;


  auto X = Pimpact::create<MSF>( space );
  auto B = Pimpact::create<MSF>( space );

  auto temp = Pimpact::create<MVF>( space );

  X->init(0.);
  B->random();


  auto H =
      Pimpact::createMultiOperatorBase(
          Pimpact::createDtLapOp( space, 0., 10. ) );

  H->apply( *temp, *temp );

  // Make an empty new parameter list.
  //auto solverParams = Teuchos::parameterList();
  auto solverParams = Pimpact::createLinSolverParameter( "GMRES", 1.e-1 );

  // Create the Pimpact::LinearSolver solver.
  auto Hprob = Pimpact::createLinearProblem<MVF>( H, temp, temp, solverParams,"GMRES" );
  auto Hinv  = Pimpact::createInverseOperator( Hprob );

  auto schur = Pimpact::createTripleCompositionOp(
      Pimpact::createMultiModeOpWrap( Pimpact::create<Pimpact::DivOp>( space ) ),
      Hinv,
      Pimpact::createMultiModeOpWrap( Pimpact::create<Pimpact::GradOp>( space ))
  );

  schur->apply( *B, *X );

}




TEUCHOS_UNIT_TEST( MultiModeOperator, InverseOperator ) {

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
      Pimpact::create<Pimpact::HelmholtzOp>(space) );

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





TEUCHOS_UNIT_TEST( MultiModeOperator, EddyPrec ) {

  //pl->set("alpha2",1.);
	//pl->set("Re",1./14.);
  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;


  auto temp = Pimpact::createMultiModeVectorField( space );

  auto X = Pimpact::createMultiModeVectorField( space );
  auto B = Pimpact::createMultiModeVectorField( space );

  X->init(0.);
  B->random();

  // Make an empty new parameter list.
  auto solverParams = Pimpact::createLinSolverParameter("CG",1.e-1);

  // Create the Pimpact::LinearSolver solver.
	auto A =
		Pimpact::createMultiOperatorBase(
				Pimpact::create<Pimpact::HelmholtzOp>( space )
				);

	auto prob =
		Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::VectorField<SpaceT> > >(
				A,
				Teuchos::null,
				Teuchos::null,
				solverParams,
				"CG" );

	auto op =
		Pimpact::createEddyPrec(
				Pimpact::createMultiOpUnWrap(
					Pimpact::createInverseOperatorBase(prob) ) );

  op->apply( X->getField(0), B->getField(0) );

	auto schur = Pimpact::createMultiOperatorBase( op );

	schur->apply( *B, *X );

}




TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicConvectionOp ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto mv1 = Pimpact::createMultiHarmonicVectorField( space, 10 );
  auto mv2 = Pimpact::createMultiHarmonicVectorField( space, 10 );

  auto op = Pimpact::createMultiHarmonicConvectionOp( space, 10 );

  op->assignField( *mv1 );
  op->apply( *mv1, *mv2 );

}



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicOpWrap ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );

  if( !isImpactInit ) isImpactInit=true;

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto mv1 = Pimpact::createMultiHarmonicVectorField( space, 10 );
  auto mv2 = Pimpact::createMultiHarmonicVectorField( space, 10 );

  auto op = Pimpact::createMultiHarmonicOpWrap< Pimpact::HelmholtzOp<SpaceT> >(
      Pimpact::create<Pimpact::HelmholtzOp>(space) );

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

  auto x =
		Pimpact::createMultiField( Pimpact::createCompoundField(
					Pimpact::createMultiHarmonic<VF>( space ),
      		Pimpact::createMultiHarmonic<SF>( space )) );
  auto fu   = x->clone();
  x->init(1.);
  x->random();

  auto opV2V =
      Pimpact::createAdd2Op(
          Pimpact::createMultiDtHelmholtz( space, 1., 1. ),
          Pimpact::createMultiHarmonicConvectionOp(space)	);

  auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
  auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );

  auto op =
      Pimpact::createMultiOperatorBase(
          Pimpact::createCompoundOpWrap(
              opV2V,
							opS2V,
							opV2S )
					);

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
