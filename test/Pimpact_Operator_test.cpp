// Pimpact_SalarVector_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include "BelosTypes.hpp"

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_Space.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinSolverParameter.hpp"

#include <iostream>
#include <cmath>


namespace {

typedef double S;
typedef int O;

bool testMpi = true;
S errorTolSlack = 1e-6;

bool isImpactInit = false;



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
  //  init_impact(0,0);
}


TEUCHOS_UNIT_TEST( BasicOperator, Div ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto p = Pimpact::createScalarField<S,O>(space);
  auto vel = Pimpact::createVectorField<S,O>(space);

  vel->initField();
  p->init(2.);

  Pimpact::Div<S,O> div;

  div.apply(*vel,*p);

  p->write(1);

  TEST_FLOATING_EQUALITY( 0., p->norm(), errorTolSlack*10.e6 )

}


TEUCHOS_UNIT_TEST( BasicOperator, Div2 ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto p = Pimpact::createScalarField<S,O>(space);
  auto vel = Pimpact::createVectorField<S,O>(space);

  vel->random();
  p->init(2.);

  Pimpact::Div<S,O> div;

  div.apply(*vel,*p);

  p->write(2);

  TEST_INEQUALITY( 0., p->norm() )
}


TEUCHOS_UNIT_TEST( BasicOperator, Grad ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto p = Pimpact::createScalarField<S,O>(space);
  auto vel = Pimpact::createVectorField<S,O>(space);

  Pimpact::Grad<S,O> grad;

  grad.apply(*p,*vel);


  TEST_EQUALITY( 0, 0 );
}



TEUCHOS_UNIT_TEST( BasicOperator, HelmholtzOp ) {

  //  if( !isImpactInit ) {
  //    init_impact(0,0);
  //    isImpactInit=true;
  //  }
  auto pl = Teuchos::parameterList();

  pl->set( "domain", 1);

  pl->set( "lx", 2. );
  pl->set( "ly", 2. );
  pl->set( "lz", 1. );


  pl->set("nx", 49 );
  pl->set("ny", 49 );
  pl->set("nz", 2 );

  // processor grid size
  pl->set("npx", 3 );
  pl->set("npy", 3 );
  pl->set("npz", 1 );

  auto space = Pimpact::createSpace( pl );

  auto x = Pimpact::createVectorField(space);
  auto bs= Pimpact::createVectorField(space);
  auto b = Pimpact::createVectorField(space);

  double mulI = 0.;
  double mulL = 1.;

  auto op= Pimpact::createHelmholtzOp( space, mulI, mulL );

//  op->print();

  // test in x direction
  x->initField( Pimpact::Poiseuille2D_inX );
  bs->init( Teuchos::tuple( 8./std::pow(space->getDomain()->getDomainSize()->getSize(Pimpact::Y),2), 0., 0. ) );
  bs->add( mulI, *x, mulL, *bs );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);
  b->write();
  bs->write(1);

  std::cout << "error: " << bs->norm() << "\n";
  TEST_EQUALITY( bs->norm()<errorTolSlack, true );

  // test in y direction
  x->initField( Pimpact::Poiseuille2D_inY );
  bs->init( Teuchos::tuple( 0., 8./std::pow(space->getDomain()->getDomainSize()->getSize(Pimpact::X),2), 0. ) );
  bs->add( mulI, *x, mulL, *bs );

  op->apply( *x, *b );

  bs->add( 1., *bs, -1, *b);
  b->write(2);
  bs->write(3);

  std::cout << "error: " << bs->norm() << "\n";
  TEST_EQUALITY( bs->norm()<errorTolSlack, true );

  // the circle test
  x->initField( Pimpact::Circle2D );
  bs->init( 0. );

  op->apply( *x, *b );

//  bs->add( 1., *bs, -1, *b );
  x->write(4);
  b->write(5);

  std::cout << "error: " << bs->norm() << "\n";
  TEST_EQUALITY( b->norm()<errorTolSlack, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, InverseHelmholtz ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto rhs = Pimpact::createVectorField(space);
  auto sol = Pimpact::createVectorField(space);
  auto x   = Pimpact::createVectorField(space);

  sol->random();
  //sol->initField( Pimpact::EFlowProfile(1) );

  auto hel = Pimpact::createHelmholtzOp( space, 1.,1.);

  hel->apply( *sol, *rhs );

  auto op = Pimpact::createInverseHelmholtzOp( 1., 1., 1.e-4, 1e4, true, false, false );

  op->apply( *rhs, *x );

  x->add(1.,*x, -1, *sol);

  x->write(77);

  //  TEST_EQUALITY( x->norm()<=1.e-3, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, MGVHelmholtzOp ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto rhs = Pimpact::createVectorField( space );
  auto sol = Pimpact::createVectorField( space );
  auto x   = Pimpact::createVectorField( space );

  //  sol->random();
  sol->initField( Pimpact::EFlowProfile(8) );
  x->initField( Pimpact::EFlowProfile(8) );

  auto hel = Pimpact::createHelmholtzOp( space, 1.,1.);

  hel->apply( *sol, *rhs );

  auto op = Pimpact::createMGVHelmholtzOp<S,O>(1.,1.,false);

  //  rhs->init( 1. );
  x->init( 0. );

  for( int i=0; i<100; ++i )
    op->apply( *rhs, *x );


  x->add(1.,*x, -1, *sol);

  x->write(777);

  //  TEST_EQUALITY( x->norm()<=1.e-2, true );

}


TEUCHOS_UNIT_TEST( BasicOperator, MGVDivGradOp ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto u = Pimpact::createVectorField<S,O>(space);
  auto temp = Pimpact::createVectorField<S,O>(space);

  auto rhs = Pimpact::createScalarField<S,O>(space);
  auto sol = Pimpact::createScalarField<S,O>(space);

  auto op = Pimpact::createMGVDivGradOp<S,O>(false);

  auto lap = Pimpact::createHelmholtzOp<S,O>( space, 0., 1. );

  auto div = Pimpact::createDivOp<S,O>();

  u->initField( Pimpact::Poiseuille2D_inX );
  u->write();


  lap->apply( *u, *temp );
  temp->write(1);
  div->apply( *temp, *rhs );
  rhs->write(2);

  sol->init( 0. );

  for( int i=0; i<100; ++i )
    op->apply( *rhs,*sol);

  sol->write(777);

}


TEUCHOS_UNIT_TEST( BasicOperator, ForcingOp ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto vel = Pimpact::createVectorField<S,O>(space);
  auto force = Pimpact::createVectorField<S,O>(space);
  auto res = Pimpact::createVectorField<S,O>(space);

  vel->init(1.);
  force->initField( Pimpact::EFlowProfile(11) );
  res->random();

  auto op = Pimpact::createForcingOp( force, 1. );

  op->apply(*vel,*res);
  vel->add( 1., *res, -1., *force );

  TEST_EQUALITY( vel->norm(), 0 );
}



TEUCHOS_UNIT_TEST( ModeOperator, HelmholtzOp ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto velc = Pimpact::createVectorField<S,O>(space);
  auto vels = Pimpact::createVectorField<S,O>(space);

  auto vel = Pimpact::createModeField( velc, vels );

  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<S,O> > >(*vel,10);

  auto mv2 = mv->clone(11);

  Pimpact::HelmholtzOp<S,O> helmholtz(space);
  auto helm = Pimpact::createMultiOpWrap<Pimpact::ModeOpWrap<Pimpact::HelmholtzOp<S,O> > >(
      Pimpact::createModeOpWrap( Pimpact::createHelmholtzOp<S,O>(space) ) );

  helm->apply(*mv,*mv2);
}



TEUCHOS_UNIT_TEST( BasicOperator, Nonlinear ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::MultiField<VF> MVF;
  //  typedef Pimpact::Nonlinear<S,O>  Op;

  auto space = Pimpact::createSpace();

  auto vel = Pimpact::createVectorField<S,O>(space);

  auto x = Pimpact::createMultiField<VF>(*vel->clone(),10);
  auto y = Pimpact::createMultiField<VF>(*vel->clone(),10);

  auto op = Pimpact::createMultiOperatorBase<MVF>( Pimpact::createNonlinear() );

  for( O i=0; i<10; ++i ) {
    x->getFieldPtr(i)->initField(Pimpact::Circle2D );
  }
  //  x->random();
  x->getFieldPtr(0)->write();

  op->apply( *x, *y);

  y->getFieldPtr(0)->write(99);
}



TEUCHOS_UNIT_TEST( BasicOperator, Add2Op ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef Pimpact::VectorField<S,O>     VF;
  typedef Pimpact::MultiField<VF>      MVF;
  typedef Pimpact::HelmholtzOp<S,O>      Op1;
  typedef Pimpact::Nonlinear<S,O>      Op2;
  typedef Pimpact::Add2Op<Op1,Op2> COp;
  typedef Pimpact::MultiOpWrap<COp>    Op;

  auto space = Pimpact::createSpace();

  auto vel = Pimpact::createVectorField<S,O>(space);


  auto x = Pimpact::createMultiField<VF>(*vel->clone(),10);
  auto y = Pimpact::createMultiField<VF>(*vel->clone(),10);

  //  auto op = Pimpact::createOperatorBase<MVF,COP>();
  auto op = Pimpact::createOperatorBase<MVF,Op>(
      Pimpact::createMultiOpWrap(
          Pimpact::createAdd2Op<Op1,Op2>(
              //              Teuchos::rcp( new Pimpact::HelmholtzOp<S,O>() ),
              Pimpact::createHelmholtzOp<S,O>( space ),
              Pimpact::createNonlinear<S,O>(),
              vel->clone() ) ) );
  //  auto op = Pimpact::createOperatorMV<O(void);

  for( O i=0; i<10; ++i ) {
    x->getFieldPtr(i)->initField(Pimpact::RankineVortex2D );
  }

  x->getFieldPtr(0)->write();

  op->apply( *x, *y);

  y->getFieldPtr(0)->write(99);

}




TEUCHOS_UNIT_TEST( ModeOperator, DtModeOp ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  typedef Pimpact::MultiField< Pimpact::ModeField< Pimpact::VectorField<S,O> > > MF;
  typedef Pimpact::OperatorBase<MF> BOp;

  auto space = Pimpact::createSpace();

  auto velc = Pimpact::createVectorField<S,O>(space);
  auto vels = Pimpact::createVectorField<S,O>(space);

  auto vel = Pimpact::createModeField( velc, vels );

  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<S,O> > >(*vel,1);

  mv->getFieldPtr(0)->getCFieldPtr()->init(1.);
  mv->getFieldPtr(0)->getSFieldPtr()->init(0.);

  auto mv2 = mv->clone(1);
  mv2->init(0.);

  auto A = Pimpact::createMultiOperatorBase<MF,Pimpact::DtModeOp<S,O> >(
      Pimpact::createDtModeOp<S,O>() );

  A->apply(*mv,*mv2);

  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm() );
  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm() );
  Belos::OperatorTraits<S,MF,BOp>::Apply(*A,*mv,*mv2);
}



TEUCHOS_UNIT_TEST( ModeOperator, DtLapOp ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto velc = Pimpact::createVectorField<S,O>(space);
  auto vels = Pimpact::createVectorField<S,O>(space);

  auto vel = Pimpact::createModeField( velc, vels );

  auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<S,O> > >(*vel,1);

  mv->random();

  auto mv2 = mv->clone(1);
  auto mv3 = mv->clone(1);
  mv2->init(0.);
  mv3->init(0.);

  auto A = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp<S,O>( space, 1., 0. ) );

  A->apply(*mv,*mv2);

  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm() );
  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm() );
  //	TEST_EQUALITY( mv->getConstField(0).getConstFieldC()->norm(), mv2->getConstField(0).getConstFieldS()->norm() );
  //	TEST_EQUALITY( mv->getConstField(0).getConstFieldS()->norm(), mv2->getConstField(0).getConstFieldC()->norm() );

  mv2->init(0.);

  auto A2 = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp<S,O>( space, 0., 1. ) );
  auto A3 = Pimpact::createMultiOpWrap( Pimpact::createModeOpWrap( Pimpact::createHelmholtzOp<S,O>( space, 0., 1. ) ) );

  A2->apply(*mv,*mv2);
  A3->apply(*mv,*mv3);

  std::vector<S> norm2( 1, 0. );
  std::vector<S> norm3( 1, 0. );

  mv2->norm(norm2);
  mv3->norm(norm3);

  TEST_EQUALITY( norm2[0] , norm3[0] );

  //	mv2->init(0.);
  //
  //	A2 = Pimpact::createMultiOpWrap( Pimpact::createDtLapOp<S,O>( 0., 0. ) );
  //	A3 = Pimpact::createMultiOpWrap( Pimpact::createModeOpWrap( Pimpact::createHelmholtzOp<S,O>( 1., 0. ) ) );
  //
  //	A2->apply(*mv,*mv2);
  //	A3->apply(*mv,*mv3);
  //
  //	mv2->norm(norm2);
  //	mv3->norm(norm3);
  //
  //	TEST_EQUALITY( norm2[0] , norm3[0] );
}



TEUCHOS_UNIT_TEST( ModeOperator, DivDtLinvGrad ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<S,O> > > MVF;
  typedef Pimpact::DtLapOp<S,O> Op;

  auto space = Pimpact::createSpace();

  auto temp = Pimpact::createMultiModeVectorField<S,O>();

  auto X = Pimpact::createMultiModeScalarField<S,O>();
  auto B = Pimpact::createMultiModeScalarField<S,O>();

  X->init(0.);
  B->random();

  //	Pimpact::createOperatorB
  auto A = Pimpact::createMultiOperatorBase<MVF,Op>(
      ( Pimpact::createDtLapOp<S,O>( space, 0., 10. ) ) );

  A->apply( *temp, *temp );

  // Make an empty new parameter list.
  RCP<ParameterList> solverParams = parameterList();

  // Set some Belos parameters.
  //
  // "Num Blocks" = Maximum number of Krylov vectors to store.  This
  // is also the restart length.  "Block" here refers to the ability
  // of this particular solver (and many other Belos solvers) to solve
  // multiple linear systems at a time, even though we are only solving
  // one linear system in this example.
  solverParams->set ("Num Blocks", 40);
  solverParams->set ("Maximum Iterations", 400);
  solverParams->set ("Convergence Tolerance", 1.0e-1);
  solverParams->set ("Output Frequency", 50);
  solverParams->set ("Output Style", 1);
  solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings +
      Belos::TimingDetails + Belos::StatusTestDetails);

  // Create the Pimpact::LinearSolver solver.
  auto H_prob = Pimpact::createLinearProblem<MVF>( A, temp, temp, solverParams,"GMRES" );

  auto schur = Pimpact::createMultiOpWrap( Pimpact::createDivDtLinvGrad<S,O>( temp, H_prob ) );
  //
  schur->apply( *B, *X );
}



TEUCHOS_UNIT_TEST( ModeOperator, TripleCompostion) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<S,O> > > MVF;
  typedef Pimpact::DtLapOp<S,O> Op;

  auto temp = Pimpact::createMultiModeVectorField<S,O>();

  auto X = Pimpact::createMultiModeScalarField<S,O>();
  auto B = Pimpact::createMultiModeScalarField<S,O>();

  X->init(0.);
  B->random();

  //  Pimpact::createOperatorB
  auto space = Pimpact::createSpace();

  auto A = Pimpact::createMultiOperatorBase<MVF,Op>(
      ( Pimpact::createDtLapOp<S,O>( space, 0., 10. ) ) );

  A->apply( *temp, *temp );

  // Make an empty new parameter list.
  auto solverParams = parameterList();

  // Create the Pimpact::LinearSolver solver.
  auto Hprob = Pimpact::createLinearProblem<MVF>( A, temp, temp, solverParams,"GMRES" );
  auto Hinv  = Pimpact::createInverseOperator( Hprob );

  auto schur = Pimpact::createTripleCompositionOp(
      Pimpact::createMultiField( temp->getConstFieldPtr(0)->clone() ),
      Pimpact::createMultiField( temp->getConstFieldPtr(0)->clone() ),
      Pimpact::createMultiModeOpWrap( Pimpact::createGradOp<S,O>()),
      Hinv,
      Pimpact::createMultiModeOpWrap( Pimpact::createDivOp<S,O>() )
  );

  schur->apply( *B, *X );
}



TEUCHOS_UNIT_TEST( ModeOperator, InverseOperator ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<S,O> > > BVF;
  typedef Pimpact::HelmholtzOp<S,O> Op;
  typedef Pimpact::InverseOperator<BVF> Op2;

  auto space = Pimpact::createSpace();

  auto X = Pimpact::createMultiModeVectorField<S,O>();
  auto B = Pimpact::createMultiModeVectorField<S,O>();

  X->init(0.);
  B->random();

  auto op = Pimpact::createMultiModeOperatorBase<BVF,Op>(
      Pimpact::createHelmholtzOp<S,O>(space) );

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
  auto opp = Pimpact::createOperatorBase<BVF,Op2 >( opinv );

}



TEUCHOS_UNIT_TEST( ModeOperator, DivOpGrad ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<S,O> > > MVF;
  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<S,O> > > MSF;
  typedef Pimpact::DtLapOp<S,O>  Op;
  typedef Pimpact::DivOpGrad<S,O>  Op2;

  auto temp = Pimpact::createMultiModeVectorField<S,O>();

  auto X = Pimpact::createMultiModeScalarField<S,O>();
  auto B = Pimpact::createMultiModeScalarField<S,O>();

  X->init(0.);
  B->random();

  // Make an empty new parameter list.
  RCP<ParameterList> solverParams = Pimpact::createLinSolverParameter("GMRES",1.e-1);

  // Create the Pimpact::LinoearSolver solver.
  auto space = Pimpact::createSpace();

  auto A = Pimpact::createMultiOperatorBase<MVF,Op>( Pimpact::createDtLapOp<S,O>( space, 14., 1. ) );

  A->apply( *temp, *temp );

  auto prob = Pimpact::createLinearProblem<MVF>( A, temp, temp, solverParams,"GMRES" );

  prob->solve(temp,temp);

  auto op = Pimpact::createDivOpGrad<S,O>( temp, prob ) ;

  op->apply( X->getField(0), B->getField(0) );

  auto schur = Pimpact::createMultiOperatorBase<MSF,Op2>( op );

  schur->apply( *B, *X );
}


TEUCHOS_UNIT_TEST( ModeOperator, EddyPrec ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<S,O> > > MVF;
  typedef Pimpact::HelmholtzOp<S,O>  Op;
  typedef Pimpact::EddyPrec<S,O>  Op2;

  auto temp = Pimpact::createMultiModeVectorField<S,O>();

  auto X = Pimpact::createMultiModeVectorField<S,O>();
  auto B = Pimpact::createMultiModeVectorField<S,O>();

  X->init(0.);
  B->random();

  // Make an empty new parameter list.
  RCP<ParameterList> solverParams = Pimpact::createLinSolverParameter("CG",1.e-1);

  // Create the Pimpact::LinearSolver solver.
  auto space = Pimpact::createSpace();

  auto A = Pimpact::createMultiModeOperatorBase<MVF,Op>( Pimpact::createHelmholtzOp<S,O>( space, 14.,1.) );

  A->apply( *temp, *temp );

  auto prob = Pimpact::createLinearProblem<MVF>( A, temp, temp, solverParams,"CG" );

  prob->solve(temp,temp);

  auto op = Pimpact::createEddyPrec<S,O>( temp, Pimpact::createInverseOperatorBase(prob) ) ;

  op->apply( X->getField(0), B->getField(0) );

  auto schur = Pimpact::createMultiOperatorBase<MVF,Op2>( op );

  schur->apply( *B, *X );
}




TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicNonlinear ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto vel = Pimpact::createVectorField<S,O>( space );

  auto mv1 = Pimpact::createMultiHarmonicVectorField<S,O>( space, 10 );
  auto mv2 = Pimpact::createMultiHarmonicVectorField<S,O>( space, 10 );

  auto op = Pimpact::createMultiHarmonicNonlinear<S,O>();

  op->apply( *mv1, *mv2 );
}



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicOpWrap ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto vel = Pimpact::createVectorField<S,O>( space );

  auto mv1 = Pimpact::createMultiHarmonicVectorField<S,O>( space, 10 );
  auto mv2 = Pimpact::createMultiHarmonicVectorField<S,O>( space, 10 );

  auto op = Pimpact::createMultiHarmonicOpWrap< Pimpact::HelmholtzOp<S,O> >(
      Pimpact::createHelmholtzOp<S,O>(space) );

  op->apply( *mv1, *mv2 );
}



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiDtHelmholtz) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto space = Pimpact::createSpace();

  auto mv1 = Pimpact::createMultiHarmonicVectorField<S,O>( space, 10 );
  auto mv2 = Pimpact::createMultiHarmonicVectorField<S,O>( space, 10 );

  auto op = Pimpact::createMultiDtHelmholtz<S,O>( space );

  op->apply( *mv1, *mv2 );
}



TEUCHOS_UNIT_TEST( CompoundOperator, CompoundOpWrap ) {

  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  typedef Pimpact::MultiHarmonicField< Pimpact::VectorField<S,O> > VF;
  typedef Pimpact::MultiHarmonicField< Pimpact::ScalarField<S,O> > SF;
  typedef Pimpact::CompoundField< VF, SF> CF;
  typedef Pimpact::MultiField<CF> MF;


  auto space = Pimpact::createSpace();

  int nfs = 4;

  auto x    = Pimpact::createMultiField( Pimpact::createCompoundField(
      Pimpact::createMultiHarmonicVectorField<S,O>( space, nfs ),
      Pimpact::createMultiHarmonicScalarField<S,O>( space, nfs )) );
  auto fu   = x->clone();
  x->init(1.);
  x->random();

  auto opV2V =
      Pimpact::createAdd2Op(
          Pimpact::createMultiDtHelmholtz<S,O>( space, 1., 1. ),
          Pimpact::createMultiHarmonicNonlinear<S,O>(),
          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone()
      );

  auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::createGradOp<S,O>() );
  auto opV2S = Pimpact::createMultiHarmonicOpWrap< Pimpact::Div<S,O> >( Pimpact::createDivOp<S,O>() );

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
