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

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_FieldFactory.hpp"

//#include "Pimpact_EddyPrec.hpp"
//#include "Pimpact_Nonlinear.hpp"
//#include "Pimpact_DivOpGrad.hpp"
//#include "Pimpact_CompoundOp.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorMV.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinSolverParameter.hpp"
#include "BelosPimpactAdapter.hpp"

#include <iostream>
#include <cmath>


namespace {

bool testMpi = true;
double errorTolSlack = 1e+1;

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


TEUCHOS_UNIT_TEST( Operator, Div ) {
  // init impact
  init_impact(0,0);

  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto p = Pimpact::createScalarField<double,int>(fS);
  auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  vel->init(1.);
	p->init(2.);

	Pimpact::Div<double,int> div;

	div.apply(*vel,*p);

	p->write(1);

//		TEST_EQUALITY( 0., p->norm() );
}


TEUCHOS_UNIT_TEST( Operator, Div2 ) {
  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto p = Pimpact::createScalarField<double,int>(fS);
  auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  vel->initField();
  p->init(2.);

  Pimpact::Div<double,int> div;

  div.apply(*vel,*p);

  p->write(2);

  TEST_INEQUALITY( 0., p->norm() );
}


TEUCHOS_UNIT_TEST( Operator, Grad ) {

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto p = Pimpact::createScalarField<double,int>(fS);
	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	Pimpact::Grad<double,int> grad;

	grad.apply(*p,*vel);


	TEST_EQUALITY( 0, 0 );
}


TEUCHOS_UNIT_TEST( Operator, Helmholtz ) {
  auto fS = Pimpact::createFieldSpace<int>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

  auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
  auto res = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

  Pimpact::Helmholtz<double,int> helmholtz;

  helmholtz.apply(*vel,*res);
}



TEUCHOS_UNIT_TEST( Operator, HelmholtzMVMode ) {
	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,10);

	auto mv2 = mv->clone(11);

	Pimpact::Helmholtz<double,int> helmholtz;
	auto helm = Pimpact::createMultiOpWrap<Pimpact::ModeOpWrap<Pimpact::Helmholtz<double,int> > >();

	helm->apply(*mv,*mv2);
}


TEUCHOS_UNIT_TEST( Operator, Dt ) {

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,1);

	mv->getFieldPtr(0)->getCFieldPtr()->init(1.);
	mv->getFieldPtr(0)->getSFieldPtr()->init(0.);

	auto mv2 = mv->clone(1);
	mv2->init(0.);

	auto A = Pimpact::createOperatorMV<Pimpact::Dt<double,int> >();

	A->apply(*mv,*mv2);

	TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm() );
	TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm() );
	Belos::OperatorTraits<double,Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >, Pimpact::OperatorMV<Pimpact::Dt<double,int> > >::Apply(*A,*mv,*mv2);
}


TEUCHOS_UNIT_TEST( Operator, DtL ) {

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >(*vel,1);

	mv->random();

	auto mv2 = mv->clone(1);
	auto mv3 = mv->clone(1);
	mv2->init(0.);
	mv3->init(0.);

	auto A = Pimpact::createDtL<double,int>( 1., 0., 0. );

	A->apply(*mv,*mv2);

  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstCFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstSFieldPtr()->norm() );
  TEST_EQUALITY( mv->getConstFieldPtr(0)->getConstSFieldPtr()->norm(), mv2->getConstFieldPtr(0)->getConstCFieldPtr()->norm() );
//	TEST_EQUALITY( mv->getConstField(0).getConstFieldC()->norm(), mv2->getConstField(0).getConstFieldS()->norm() );
//	TEST_EQUALITY( mv->getConstField(0).getConstFieldS()->norm(), mv2->getConstField(0).getConstFieldC()->norm() );

	Belos::OperatorTraits<double,Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >, Pimpact::OperatorMV<Pimpact::DtL<double,int> > >::Apply(*A,*mv,*mv2);

	mv2->init(0.);

	auto A2 = Pimpact::createDtL<double,int>( 0., 0., 1. );
	auto A3 = Pimpact::createHelmholtzdep<double,int>( 0., 1. );

	A2->apply(*mv,*mv2);
	A3->apply(*mv,*mv3);

	std::vector<double> norm2( 1, 0. );
	std::vector<double> norm3( 1, 0. );

	mv2->norm(norm2);
	mv3->norm(norm3);

	TEST_EQUALITY( norm2[0] , norm3[0] );

	mv2->init(0.);

	A2 = Pimpact::createDtL<double,int>( 0., 1., 0. );
	A3 = Pimpact::createHelmholtzdep<double,int>( 1., 0. );

	A2->apply(*mv,*mv2);
	A3->apply(*mv,*mv3);

	mv2->norm(norm2);
	mv3->norm(norm3);

	TEST_EQUALITY( norm2[0] , norm3[0] );

}


TEUCHOS_UNIT_TEST( Operator, Div_DtLinv_Grad ) {
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing

	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MVF;
	typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > > MSF;
	typedef Pimpact::OperatorMV< Pimpact::DtL<double,int> >  OP;

	auto temp = Pimpact::createMultiModeVectorField<double,int>();

	auto X = Pimpact::createMultiModeScalarField<double,int>();
	auto B = Pimpact::createMultiModeScalarField<double,int>();

	X->init(0.);
	B->random();

	auto A = Pimpact::createDtL<double,int>(0.,0.,10.);

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
	auto H_prob = Pimpact::createLinearProblem<double,MVF,OP>( A, temp, temp, solverParams,"GMRES" );

	auto schur = Pimpact::createDivDtLinvGrad<double,int>( temp, H_prob );

	schur->apply( *B, *X );
}


TEUCHOS_UNIT_TEST( Operator, Linv ) {
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef double S;
  typedef int    O;
  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<S,O> > > BVF;
  typedef Pimpact::MultiField<Pimpact::VectorField<S,O> > MVF;
  typedef Pimpact::OperatorMV< Pimpact::Helmholtz<S,O> >  OP;
  typedef Pimpact::OperatorBase<MVF> BOP;


  auto X = Pimpact::createMultiModeVectorField<S,O>();
  auto B = Pimpact::createMultiModeVectorField<S,O>();

  X->init(0.);
  B->random();

  auto op = Pimpact::createHelmholtzdep<S,O>( 0., 1. );

  // Make an empty new parameter list.
  auto solverParams = Pimpact::createLinSolverParameter( "CG", 1.e-1 )->get();


// Create the Pimpact::LinearSolver solver.
  auto prob =
      Pimpact::createLinearProblem<S,MVF,OP>(
          op,
          Pimpact::createMultiField(X->getFieldPtr(0)->getCFieldPtr()),
          Pimpact::createMultiField(B->getFieldPtr(0)->getCFieldPtr()), solverParams );


  auto opp = Pimpact::createOperatorBaseMV<BVF,Pimpact::Linv<S,O> >( Pimpact::createLinv<S,O>( prob ) );

  opp->apply( *B, *X );
}


TEUCHOS_UNIT_TEST( Operator, DivOpGrad ) {
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef double S;
  typedef int O;
  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MVF;
  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > > MSF;
  typedef Pimpact::DtL<double,int>  OP;
  typedef Pimpact::OperatorBase<MVF>  BOP;
  typedef Pimpact::OperatorBase<MSF>  BSOP;
  typedef Pimpact::DivOpGrad<S,O>  OP2;

  auto temp = Pimpact::createMultiModeVectorField<double,int>();

  auto X = Pimpact::createMultiModeScalarField<S,O>();
  auto B = Pimpact::createMultiModeScalarField<S,O>();

  X->init(0.);
  B->random();


  // Make an empty new parameter list.
  RCP<ParameterList> solverParams = Pimpact::createLinSolverParameter("GMRES",1.e-1)->get();

// Create the Pimpact::LinearSolver solver.
  auto A = Pimpact::createOperatorBaseMV<MVF,OP>( Pimpact::createDtL<S,O>(14.,0.,1.) );

  A->apply( *temp, *temp );

  auto prob = Pimpact::createLinearProblem<S,MVF,BOP>( A, temp, temp, solverParams,"GMRES" );

  prob->solve(temp,temp);

  auto op = Pimpact::createDivOpGrad<S,O>( temp, prob ) ;

  op->apply( X->getField(0), B->getField(0) );

  auto schur = Pimpact::createOperatorBasedep<MSF,OP2>( op );

  schur->apply( *B, *X );
}


TEUCHOS_UNIT_TEST( Operator, EddyPrec ) {
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef double S;
  typedef int O;
  typedef Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > > MVF;
  typedef Pimpact::Helmholtz<S,O>  OP;
  typedef Pimpact::EddyPrec<S,O>  OP2;
  typedef Pimpact::OperatorBase<MVF>  BOP;

  auto temp = Pimpact::createMultiModeVectorField<S,O>();

  auto X = Pimpact::createMultiModeVectorField<S,O>();
  auto B = Pimpact::createMultiModeVectorField<S,O>();

  X->init(0.);
  B->random();

  // Make an empty new parameter list.
  RCP<ParameterList> solverParams = Pimpact::createLinSolverParameter("CG",1.e-1)->get();

  // Create the Pimpact::LinearSolver solver.
  auto A = Pimpact::createOperatorBaseMV<MVF,OP>( Pimpact::createHelmholtzdep<S,O>(14.,1.) );

  A->apply( *temp, *temp );

  auto prob = Pimpact::createLinearProblem<S,MVF,BOP>( A, temp, temp, solverParams,"CG" );

  prob->solve(temp,temp);

  auto op = Pimpact::createEddyPrec<S,O>( temp, prob ) ;

  op->apply( X->getField(0), B->getField(0) );

  auto schur = Pimpact::createOperatorBasedep<MVF,OP2>( op );

  schur->apply( *B, *X );
}


TEUCHOS_UNIT_TEST( Operator, nonlinear ) {
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef double S;
  typedef int O;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::MultiField<VF> MVF;
  typedef Pimpact::Nonlinear<S,O>  OP;
  typedef Pimpact::OperatorBase<MVF>  BOP;

  auto fS = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto vel = Pimpact::createVectorField<S,O>(fS,iIS,fIS);


  auto x = Pimpact::createMultiField<VF>(*vel->clone(),10);
  auto y = Pimpact::createMultiField<VF>(*vel->clone(),10);

  auto op = Pimpact::createOperatorBasedep<MVF,OP>();
//  auto op = Pimpact::createOperatorMV<O(void);

  for( int i=0; i<10; ++i ) {
    x->getFieldPtr(i)->initField(Pimpact::Circle2D );
  }
//  x->random();
x->getFieldPtr(0)->write();

  op->apply( *x, *y);

  y->getFieldPtr(0)->write(99);

}


TEUCHOS_UNIT_TEST( Operator, compoundOp ) {
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef double S;
  typedef int O;
  typedef Pimpact::VectorField<S,O>     VF;
  typedef Pimpact::MultiField<VF>      MVF;
  typedef Pimpact::Helmholtz<S,O>      Op1;
  typedef Pimpact::Nonlinear<S,O>      Op2;
  typedef Pimpact::CompoundOp<Op1,Op2> COp;
  typedef Pimpact::MultiOpWrap<COp>    Op;
  typedef Pimpact::OperatorBase<MVF>   BOp;

  auto fS = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto vel = Pimpact::createVectorField<S,O>(fS,iIS,fIS);


  auto x = Pimpact::createMultiField<VF>(*vel->clone(),10);
  auto y = Pimpact::createMultiField<VF>(*vel->clone(),10);

//  auto op = Pimpact::createOperatorBase<MVF,COP>();
  auto op = Pimpact::createOperatorBase<MVF,Op>(
      Pimpact::createMultiOpWrap(
          Pimpact::createCompoundOp<Op1,Op2>(
//              Teuchos::rcp( new Pimpact::Helmholtz<S,O>() ),
                  Pimpact::createHelmholtz<S,O>(),
                  Pimpact::createNonlinear<S,O>(),
                  vel->clone() ) ) );
//  auto op = Pimpact::createOperatorMV<O(void);

  for( int i=0; i<10; ++i ) {
    x->getFieldPtr(i)->initField(Pimpact::RankineVortex2D );
  }

  x->getFieldPtr(0)->write();

  op->apply( *x, *y);

  y->getFieldPtr(0)->write(99);

}


} // namespace
