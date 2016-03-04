#include <cmath>
#include <iostream>
#include <vector>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "BelosTypes.hpp"

#include "NOX.H"

#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Group.hpp"
#include "NOX_Pimpact_StatusTest.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"




namespace {


using ST = double;
using OT = int;

const int d = 3;
const int dNC = 2;

using SpaceT = Pimpact::Space<ST,OT,d,dNC>;


bool testMpi = true;
ST eps = 1.e0;
auto pl = Teuchos::parameterList();

template<class T> using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;


TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used." );
  clp.setOption(
      "eps", &eps,
      "epsilon used to check test results" );

  pl->set( "dim", 3 );
  pl->set( "domain", 0 );

  pl->set("nx", 65 );
  pl->set("ny", 65 );
  pl->set("nz", 9 );

  pl->set("lx", 2. );
  pl->set("ly", 2. );

  pl->set("Re", 1. );
}



//TEUCHOS_UNIT_TEST( NOXPimpact_LinearStokes, createLinearStokes ) {
//
//
//  typedef Pimpact::ScalarField<ST,OT> SF;
//  typedef Pimpact::VectorField<ST,OT> VF;
//  typedef Pimpact::ModeField<SF> MSF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField<MSF> BMSF;
//  typedef Pimpact::MultiField<MVF> BMVF;
//  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
//  typedef NOX::Pimpact::Vector<CF> NV;
//  typedef NOX::Pimpact::LinearStokes Interface;
//
//  auto fS  = Pimpact::createStencilWidths<OT>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<OT>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<OT>();
//
//  auto xv = Pimpact::createInitMVF<ST,OT>(Pimpact::Zero2DFlow, fST, iIS, fIS );
//
//  auto xs = Pimpact::createInitMSF<ST,OT>( fS );
//
//  auto x  = Pimpact::createCompoundField( xv, xs );
//
//  auto dtL = Pimpact::createMultiOperatorBase<BMVF,Pimpact::DtL<ST,OT> >( Pimpact::createDtL<ST,OT>(1.,0.,1.) );
//
//  // Make an empty new parameter list.
//  auto solverName = "GMRES";
//  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
//  solverParams->set ("Verbosity",  Belos::Errors );
//
//// Create the Pimpact::LinearSolver solver.
//  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
//      dtL, xv->clone(), xv->clone(), solverParams, solverName );
//
//  auto schur = Pimpact::createMultiOperatorBase< BMSF, Pimpact::DivDtLinvGrad<ST,OT> >(
//      Pimpact::createDivDtLinvGrad<ST,OT>( xv->clone(), lp_DTL ) );
//
//  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
//      schur, xs->clone(), xs->clone(), solverParams, solverName );
//
//  auto stockie = NOX::Pimpact::createLinearStokes(xv->clone(),xs->clone(),lp_DTL,lp_Schur);
//}



//TEUCHOS_UNIT_TEST( NOXPimpact_LinearStokes, computeF ) {
//
//  typedef double S;
//  typedef int O;
//  typedef Pimpact::ScalarField<ST,OT> SF;
//  typedef Pimpact::VectorField<ST,OT> VF;
//  typedef Pimpact::ModeField<SF> MSF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField< MSF> BMSF;
//  typedef Pimpact::MultiField< MVF> BMVF;
//  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
//  typedef NOX::Pimpact::Vector<CF> NV;
//  typedef NOX::Pimpact::LinearStokes Interface;
//
//  auto fS  = Pimpact::createStencilWidths<OT>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<OT>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<OT>();
//
//  auto xv = Pimpact::createInitMVF<ST,OT>(Pimpact::Zero2DFlow, fS, iIS, fIS );
//
//  auto xs = Pimpact::createInitMSF<ST,OT>( fS );
//
//  auto x  = Pimpact::createCompoundField( xv, xs );
//
//  Teuchos::RCP<NV> X = Teuchos::rcp(new NV(x) );
//  Teuchos::RCP<NV> F = Teuchos::rcp_dynamic_cast<NV>( X->clone() );
//
//  auto dtL = Pimpact::createMultiOperatorBase<BMVF,Pimpact::DtL<ST,OT> >( Pimpact::createDtL<ST,OT>(1.,0.,1.) );
//
//  // Make an empty new parameter list.
//  auto solverName = "GMRES";
//  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
//  solverParams->set ("Verbosity",  Belos::Errors );
//
//// Create the Pimpact::LinearSolver solver.
//  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
//      dtL, xv->clone(), xv->clone(), solverParams, solverName );
//
//  auto schur = Pimpact::createMultiOperatorBase< BMSF, Pimpact::DivDtLinvGrad<ST,OT> >(
//      Pimpact::createDivDtLinvGrad<ST,OT>( xv->clone(), lp_DTL ) );
//
//  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
//      schur, xs->clone(), xs->clone(), solverParams, solverName );
//
//  auto stockie = NOX::Pimpact::createLinearStokes(xv->clone(),xs->clone(),lp_DTL,lp_Schur);
//
//
//  NOX::Abstract::Group::ReturnType succes = stockie->computeF( X->getConstField(), F->getField() );
//
//  TEST_EQUALITY( NOX::Abstract::Group::Ok, succes);
//}



//TEUCHOS_UNIT_TEST( NOXPimpact_LinearStokes, computeJacobian ) {
//
//  typedef double S;
//  typedef int O;
//  typedef Pimpact::ScalarField<ST,OT> SF;
//  typedef Pimpact::VectorField<ST,OT> VF;
//  typedef Pimpact::ModeField<SF> MSF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField< MSF> BMSF;
//  typedef Pimpact::MultiField< MVF> BMVF;
//  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
//  typedef NOX::Pimpact::Vector<CF> NV;
//  typedef NOX::Pimpact::LinearStokes Interface;
//
//  auto fS  = Pimpact::createStencilWidths<OT>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<OT>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<OT>();
//
//  auto xv = Pimpact::createInitMVF<ST,OT>(Pimpact::Zero2DFlow, fS, iIS, fIS );
//
//  auto xs = Pimpact::createInitMSF<ST,OT>( fS );
//
//  auto x  = Pimpact::createCompoundField( xv, xs );
//
//  Teuchos::RCP<NV> X = Teuchos::rcp(new NV(x) );
//  Teuchos::RCP<NV> F = Teuchos::rcp_dynamic_cast<NV>( X->clone() );
//
//  auto dtL = Pimpact::createMultiOperatorBase<BMVF,Pimpact::DtL<ST,OT> >( Pimpact::createDtL<ST,OT>(1.,0.,1.) );
//
//  // Make an empty new parameter list.
//  auto solverName = "GMRES";
//  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
//  solverParams->set ("Verbosity",  Belos::Errors );
//
//// Create the Pimpact::LinearSolver solver.
//  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
//      dtL, xv->clone(), xv->clone(), solverParams, solverName );
//
//  auto schur = Pimpact::createMultiOperatorBase< BMSF, Pimpact::DivDtLinvGrad<ST,OT> >(
//      Pimpact::createDivDtLinvGrad<ST,OT>( xv->clone(), lp_DTL ) );
//
//  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
//      schur, xs->clone(), xs->clone(), solverParams, solverName );
//
//  auto stockie = NOX::Pimpact::createLinearStokes(xv->clone(), xs->clone(),lp_DTL,lp_Schur);
//
//
//  NOX::Abstract::Group::ReturnType succes = stockie->computeJacobian( X->getConstField() );
//
//  TEST_EQUALITY( NOX::Abstract::Group::Ok, succes);
//}



//TEUCHOS_UNIT_TEST( NOXPimpact_LinearStokes, applyJacobian ) {
//
//  typedef double S;
//  typedef int O;
//  typedef Pimpact::ScalarField<ST,OT> SF;
//  typedef Pimpact::VectorField<ST,OT> VF;
//  typedef Pimpact::ModeField<SF> MSF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField< MSF> BMSF;
//  typedef Pimpact::MultiField< MVF> BMVF;
//  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
//  typedef NOX::Pimpact::Vector<CF> NV;
//  typedef NOX::Pimpact::LinearStokes Interface;
//
//  auto fS  = Pimpact::createStencilWidths<OT>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<OT>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<OT>();
//
//  auto xv = Pimpact::createInitMVF<ST,OT>(Pimpact::Zero2DFlow, fS, iIS, fIS );
//
//  auto xs = Pimpact::createInitMSF<ST,OT>( fS );
//
//  auto x  = Pimpact::createCompoundField( xv, xs );
//
//  Teuchos::RCP<NV> X = Teuchos::rcp(new NV(x) );
//  Teuchos::RCP<NV> F = Teuchos::rcp_dynamic_cast<NV>( X->clone() );
//
//  auto dtL = Pimpact::createMultiOperatorBase<BMVF,Pimpact::DtL<ST,OT> >( Pimpact::createDtL<ST,OT>(1.,0.,1.) );
//
//  // Make an empty new parameter list.
//  auto solverName = "GMRES";
//  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
//  solverParams->set ("Verbosity",  Belos::Errors );
//
//// Create the Pimpact::LinearSolver solver.
//  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
//      dtL, xv->clone(), xv->clone(), solverParams, solverName );
//
//  auto schur = Pimpact::createMultiOperatorBase< BMSF, Pimpact::DivDtLinvGrad<ST,OT> >(
//      Pimpact::createDivDtLinvGrad<ST,OT>( xv->clone(), lp_DTL ) );
//
//  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
//      schur, xs->clone(), xs->clone(), solverParams, solverName );
//
//  auto stockie = NOX::Pimpact::createLinearStokes(xv->clone(),xs->clone(),lp_DTL,lp_Schur);
//
//
//  NOX::Abstract::Group::ReturnType succes = stockie->applyJacobian( X->getConstField(), F->getField() );
//
//  TEST_EQUALITY( NOX::Abstract::Group::Ok, succes);
//}


//TEUCHOS_UNIT_TEST( NOXPimpact_LinearStokes, applyJacobianInverse ) {
//
//  typedef double S;
//  typedef int O;
//  typedef Pimpact::ScalarField<ST,OT> SF;
//  typedef Pimpact::VectorField<ST,OT> VF;
//  typedef Pimpact::ModeField<SF> MSF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField< MSF> BMSF;
//  typedef Pimpact::MultiField< MVF> BMVF;
//  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
//  typedef NOX::Pimpact::Vector<CF> NV;
//  typedef NOX::Pimpact::LinearStokes Interface;
//
//  auto fS  = Pimpact::createStencilWidths<OT>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<OT>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<OT>();
//
//  auto xv = Pimpact::createInitMVF<ST,OT>(Pimpact::Zero2DFlow, fS, iIS, fIS );
//
//  auto xs = Pimpact::createInitMSF<ST,OT>( fS );
//
//  auto x  = Pimpact::createCompoundField( xv, xs );
//
//  Teuchos::RCP<NV> X = Teuchos::rcp(new NV(x) );
//  Teuchos::RCP<NV> F = Teuchos::rcp_dynamic_cast<NV>( X->clone() );
//
//  auto dtL = Pimpact::createMultiOperatorBase<BMVF,Pimpact::DtL<ST,OT> >( Pimpact::createDtL<ST,OT>(1.,0.,1.) );
//
//  // Make an empty new parameter list.
//  auto solverName = "GMRES";
//  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
//  solverParams->set ("Verbosity",  Belos::Errors );
//
//// Create the Pimpact::LinearSolver solver.
//  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
//      dtL, xv->clone(), xv->clone(), solverParams, solverName );
//
//  auto schur = Pimpact::createMultiOperatorBase< BMSF, Pimpact::DivDtLinvGrad<ST,OT> >(
//      Pimpact::createDivDtLinvGrad<ST,OT>( xv->clone(), lp_DTL ) );
//
//  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
//      schur, xs->clone(), xs->clone(), solverParams, solverName );
//
//  auto stockie = NOX::Pimpact::createLinearStokes(xv->clone(),xs->clone(),lp_DTL,lp_Schur);
//
//
////  NOX::Abstract::Group::ReturnType succes = stockie->applyJacobianInverse(Teuchos::parameterList(), X->getConstField(), F->getField() );
//
////  TEST_EQUALITY( NOX::Abstract::Group::Ok, succes);
//}


//TEUCHOS_UNIT_TEST( NOXPimpact_SimpleNonlinear, computeF ) {
//
//  using Teuchos::ParameterList;
//  using Teuchos::parameterList;
//  using Teuchos::RCP;
//  using Teuchos::rcp; // Save some typing
//
//  typedef double S;
//  typedef int O;
//  typedef Pimpact::VectorField<ST,OT> VF;
//  typedef Pimpact::MultiField<VF> MVF;
//  typedef Pimpact::ConvectionOp<ST,OT>  OP;
//  typedef Pimpact::OperatorBase<MVF>  BOP;
//
//  auto fS = Pimpact::createStencilWidths<OT>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<OT>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<OT>();
//
//  auto vel = Pimpact::createVectorField<ST,OT>(fS,iIS,fIS);
//
//
//  auto x = Pimpact::createMultiField<VF>(*vel->clone(),10);
//  auto y = Pimpact::createMultiField<VF>(*vel->clone(),10);
//
//  auto op = Pimpact::createMultiOperatorBase<MVF,OP>();
//
//  for( int i=0; i<10; ++i ) {
//    x->getFieldPtr(i)->initField(Pimpact::Circle2D );
//  }
////  x->random();
//x->getFieldPtr(0)->write();
//
//  op->apply( *x, *y);
//
//  y->getFieldPtr(0)->write(99);
//
//  auto lp = Pimpact::createLinearProblem<MVF>(
//        op, x->clone(), y->clone(), Pimpact::createLinSolverParameter("GMRES",1.e-12), "GMRES");
//  auto inter = NOX::Pimpact::createSimpleNonlinear( y, op, lp );
//
////  bool succes = inter->computeF( *x, *f );
//
////  TEST_EQUALITY( true, succes);
//}



//TEUCHOS_UNIT_TEST( NOXPimpact_Group, createGroup ) {
//
//  typedef double S;
//  typedef int O;
//  typedef Pimpact::ScalarField<ST,OT> SF;
//  typedef Pimpact::VectorField<ST,OT> VF;
//  typedef Pimpact::ModeField<SF> MSF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField< MSF> BMSF;
//  typedef Pimpact::MultiField< MVF> BMVF;
//  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
//  typedef NOX::Pimpact::Vector<CF> NV;
//  typedef NOX::Pimpact::Interface<> Interface;
//
//  auto fS  = Pimpact::createStencilWidths<OT>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<OT>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<OT>();
//
//  auto xv = Pimpact::createInitMVF<ST,OT>(Pimpact::Zero2DFlow, fS, iIS, fIS );
//
//  auto xs = Pimpact::createInitMSF<ST,OT>( fS );
//
//  auto x  = Pimpact::createCompoundField( xv, xs );
//
//  auto dtL = Pimpact::createMultiOperatorBase<BMVF,Pimpact::DtL<ST,OT> >( Pimpact::createDtL<ST,OT>(1.,0.,1.) );
//
//  // Make an empty new parameter list.
//  auto solverName = "GMRES";
//  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
//  solverParams->set ("Verbosity",  Belos::Errors );
//
//// Create the Pimpact::LinearSolver solver.
//  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
//      dtL, xv->clone(), xv->clone(), solverParams, solverName );
//
//  auto schur = Pimpact::createMultiOperatorBase< BMSF, Pimpact::DivDtLinvGrad<ST,OT> >(
//      Pimpact::createDivDtLinvGrad<ST,OT>( xv->clone(), lp_DTL ) );
//
//  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
//      schur, xs->clone(), xs->clone(), solverParams, solverName );
//
//  auto stockie = NOX::Pimpact::createLinearStokes(xv->clone(),xs->clone(),lp_DTL,lp_Schur);
//
//  Teuchos::RCP<NV> nx = Teuchos::rcp(new NV(x) );
//
//  auto bla = Teuchos::parameterList();
//  NOX::Pimpact::Group<NOX::Pimpact::LinearStokes> asdfasdf( *bla, stockie, *nx );
//
//  auto blabla = NOX::Pimpact::createGroup<NOX::Pimpact::LinearStokes>( bla, stockie, nx );
//}


//TEUCHOS_UNIT_TEST( NOXPimpact_Group, computeF ) {
//
//  typedef double S;
//  typedef int O;
//  typedef Pimpact::ScalarField<ST,OT> SF;
//  typedef Pimpact::VectorField<ST,OT> VF;
//  typedef Pimpact::ModeField<SF> MSF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField< MSF> BMSF;
//  typedef Pimpact::MultiField< MVF> BMVF;
//  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
//  typedef NOX::Pimpact::Vector<CF> NV;
//  typedef NOX::Pimpact::LinearStokes Interface;
//
//  auto fS  = Pimpact::createStencilWidths<OT>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<OT>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<OT>();
//
//  auto xv = Pimpact::createInitMVF<ST,OT>(Pimpact::Zero2DFlow, fS, iIS, fIS );
//
//  auto xs = Pimpact::createInitMSF<ST,OT>( fS );
//
//  auto x  = Pimpact::createCompoundField( xv, xs );
//
//  auto dtL = Pimpact::createMultiOperatorBase<BMVF,Pimpact::DtL<ST,OT> >( Pimpact::createDtL<ST,OT>(1.,0.,1.) );
//
//  // Make an empty new parameter list.
//  auto solverName = "GMRES";
//  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
//  solverParams->set ("Verbosity",  Belos::Errors );
//
//// Create the Pimpact::LinearSolver solver.
//  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
//      dtL, xv->clone(), xv->clone(), solverParams, solverName );
//
//  auto schur = Pimpact::createMultiOperatorBase< BMSF, Pimpact::DivDtLinvGrad<ST,OT> >(
//      Pimpact::createDivDtLinvGrad<ST,OT>( xv->clone(), lp_DTL ) );
//
//  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
//      schur, xs->clone(), xs->clone(), solverParams, solverName );
//
//  auto stockie = NOX::Pimpact::createLinearStokes(xv->clone(),xs->clone(),lp_DTL,lp_Schur);
//
//  Teuchos::RCP<NV> nx = Teuchos::rcp(new NV(x) );
//
//  auto bla = Teuchos::parameterList();
//  NOX::Pimpact::Group<NOX::Pimpact::LinearStokes> group( *bla, stockie, *nx );
//
//  group.computeF();
//
//}


TEUCHOS_UNIT_TEST( NOXPimpact_Group, SimpleNonlinear ) {

  typedef Pimpact::VectorField<SpaceT> VF;

	typedef Pimpact::MultiField<VF> MVF;


  typedef NOX::Pimpact::Interface<MVF> Inter;
  typedef NOX::Pimpact::Vector<typename Inter::Field> NV;


  auto space = Pimpact::createSpace<ST,OT,d,dNC>( pl );

  int rank = space->rankST();

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  vel->initField( Pimpact::ZeroFlow );

  auto x = Pimpact::createMultiField( *vel->clone(), 1 );
  auto f = Pimpact::createMultiField( *vel->clone(), 1 );


 auto sop = Pimpact::create<ConvDiffOpT>( space ); ;

 auto op =
	 Pimpact::createOperatorBase(
			 Pimpact::createMultiOpWrap(
				 sop
				 )
			 );

 auto smoother =
	 Pimpact::createOperatorBase(
			 Pimpact::createMultiOpWrap(
				 Pimpact::create<
				 Pimpact::NonlinearSmoother<
				 ConvDiffOpT<Pimpact::Space<ST,OT,d,dNC> > ,
				 Pimpact::ConvectionDiffusionSORSmoother > > (
					 sop
					 )
				 )
			 );

 auto jop = op;

 // init Fields, init and rhs
 x->getFieldPtr(0)->initField( Pimpact::RankineVortex2D );
 f->getFieldPtr(0)->initField( Pimpact::ZeroFlow );

// x->write( 97 );
 op->assignField(*x);
 op->apply( *x, *f );
// f->write( 98 );

 x->init( 0. );

 auto para = Pimpact::createLinSolverParameter("GMRES",1.e-6);

 auto lp_ = Pimpact::createLinearProblem<MVF>(
		 jop, x->clone(), f->clone(), para , "GMRES" );

 lp_->setRightPrec( smoother );
//  lp_->setLeftPrec( smoother );

 auto lp = Pimpact::createInverseOperatorBase<MVF>( lp_ );

 auto inter = NOX::Pimpact::createInterface<MVF>( f, op, lp );

 auto nx = NOX::Pimpact::createVector( x );

 auto bla = Teuchos::parameterList();

 auto group = NOX::Pimpact::createGroup( bla, inter, nx );

 // Set up the status tests
 auto statusTest = NOX::Pimpact::createStatusTest();

 // Create the list of solver parameters
 Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr =
	 Teuchos::rcp(new Teuchos::ParameterList);


 // Create the line search parameters sublist
 Teuchos::ParameterList& lineSearchParameters = solverParametersPtr->sublist("Line Search");
 lineSearchParameters.set("Method","Backtrack");
 lineSearchParameters.sublist("Backtrack").set("Recovery Step",1.e-12/2.);


 // Create the solver
 Teuchos::RCP<NOX::Solver::Generic> solver =
	 NOX::Solver::buildSolver( group, statusTest, solverParametersPtr);

 // Solve the nonlinear system
 NOX::StatusTest::StatusType status = solver->solve();

 // Print the parameter list
 if(rank==0) std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
 if(rank==0) std::cout << "\n" << status << "\n";
 if(rank==0) solver->getList().print(std::cout);

 // Get the answer
 *group = solver->getSolutionGroup();

 // Print the answer
 if(rank==0) std::cout << "\n" << "-- Final Solution From Solver --" << "\n";
// Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->write(99);

 auto sol = Teuchos::rcp_const_cast<VF>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->getConstFieldPtr(0) );
// sol->write(888);

 vel->initField( Pimpact::RankineVortex2D );

 auto er = vel->clone();
 er->initField( Pimpact::ZeroFlow );
 er->add( 1., *sol, -1., *vel );
// er->write(100);

 std::cout << " error: " << er->norm() << "\n";
 TEST_EQUALITY( er->norm() < 1.e-4, true );

// Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(999);

}



} // namespace
