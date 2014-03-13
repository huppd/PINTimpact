#include <iostream>
#include <vector>
#include <cmath>

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
//#include "Pimpact_ModeField.hpp"
#include "Pimpact_CompoundField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_EddyPrec.hpp"
#include "Pimpact_Nonlinear.hpp"
#include "Pimpact_DivOpGrad.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorMV.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"

#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_SimpleNonlinear.hpp"

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



TEUCHOS_UNIT_TEST( NOXPimpact_Interface, defaultctor ) {

  init_impact(0,0);

  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField< MSF> BMSF;
  typedef Pimpact::MultiField< MVF> BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
  typedef NOX::Pimpact ::Vector<CF> NV;

  NOX::Pimpact::Interface asdf;

}


TEUCHOS_UNIT_TEST( NOXPimpact_Interface, createInterface ) {


  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField< MSF> BMSF;
  typedef Pimpact::MultiField< MVF> BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
  typedef NOX::Pimpact::Vector<CF> NV;
  typedef NOX::Pimpact::Interface Interface;

  auto fS  = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::Zero2DFlow, fS, iIS, fIS );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

  auto x  = Pimpact::createCompoundField( xv, xs );

  auto dtL = Pimpact::createDtL<S,O>(1.,0.,1.);

  // Make an empty new parameter list.
  auto solverName = "GMRES";
  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
  solverParams->get()->set ("Verbosity",  Belos::Errors );

// Create the Pimpact::LinearSolver solver.
  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
      dtL, xv->clone(), xv->clone(), solverParams->get(), solverName );

  auto schur = Pimpact::createDivDtLinvGrad<S,O>( xv->clone(), lp_DTL );

  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
      schur, xs->clone(), xs->clone(), solverParams->get(), solverName );

  auto stockie = NOX::Pimpact::createInterface(xv->clone(),xs->clone(),lp_DTL,lp_Schur);
}



TEUCHOS_UNIT_TEST( NOXPimpact_Interface, computeF ) {

  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField< MSF> BMSF;
  typedef Pimpact::MultiField< MVF> BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
  typedef NOX::Pimpact::Vector<CF> NV;
  typedef NOX::Pimpact::Interface Interface;

  auto fS  = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::Zero2DFlow, fS, iIS, fIS );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

  auto x  = Pimpact::createCompoundField( xv, xs );

  Teuchos::RCP<NV> X = Teuchos::rcp(new NV(x) );
  Teuchos::RCP<NV> F = Teuchos::rcp_dynamic_cast<NV>( X->clone() );

  auto dtL = Pimpact::createDtL<S,O>(1.,0.,1.);

  // Make an empty new parameter list.
  auto solverName = "GMRES";
  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
  solverParams->get()->set ("Verbosity",  Belos::Errors );

// Create the Pimpact::LinearSolver solver.
  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
      dtL, xv->clone(), xv->clone(), solverParams->get(), solverName );

  auto schur = Pimpact::createDivDtLinvGrad<S,O>( xv->clone(), lp_DTL );

  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
      schur, xs->clone(), xs->clone(), solverParams->get(), solverName );

  auto stockie = NOX::Pimpact::createInterface(xv->clone(),xs->clone(),lp_DTL,lp_Schur);


  bool succes = stockie->computeF( X->getConstField(), F->getField() );

  TEST_EQUALITY( true, succes);
}



TEUCHOS_UNIT_TEST( NOXPimpact_Interface, computeJacobian ) {

  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField< MSF> BMSF;
  typedef Pimpact::MultiField< MVF> BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
  typedef NOX::Pimpact::Vector<CF> NV;
  typedef NOX::Pimpact::Interface Interface;

  auto fS  = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::Zero2DFlow, fS, iIS, fIS );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

  auto x  = Pimpact::createCompoundField( xv, xs );

  Teuchos::RCP<NV> X = Teuchos::rcp(new NV(x) );
  Teuchos::RCP<NV> F = Teuchos::rcp_dynamic_cast<NV>( X->clone() );

  auto dtL = Pimpact::createDtL<S,O>(1.,0.,1.);

  // Make an empty new parameter list.
  auto solverName = "GMRES";
  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
  solverParams->get()->set ("Verbosity",  Belos::Errors );

// Create the Pimpact::LinearSolver solver.
  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
      dtL, xv->clone(), xv->clone(), solverParams->get(), solverName );

  auto schur = Pimpact::createDivDtLinvGrad<S,O>( xv->clone(), lp_DTL );

  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
      schur, xs->clone(), xs->clone(), solverParams->get(), solverName );

  auto stockie = NOX::Pimpact::createInterface(xv->clone(), xs->clone(),lp_DTL,lp_Schur);


  bool succes = stockie->computeJacobian( X->getConstField() );

  TEST_EQUALITY( true, succes);
}



TEUCHOS_UNIT_TEST( NOXPimpact_Interface, applyJacobian ) {

  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField< MSF> BMSF;
  typedef Pimpact::MultiField< MVF> BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
  typedef NOX::Pimpact::Vector<CF> NV;
  typedef NOX::Pimpact::Interface Interface;

  auto fS  = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::Zero2DFlow, fS, iIS, fIS );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

  auto x  = Pimpact::createCompoundField( xv, xs );

  Teuchos::RCP<NV> X = Teuchos::rcp(new NV(x) );
  Teuchos::RCP<NV> F = Teuchos::rcp_dynamic_cast<NV>( X->clone() );

  auto dtL = Pimpact::createDtL<S,O>(1.,0.,1.);

  // Make an empty new parameter list.
  auto solverName = "GMRES";
  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
  solverParams->get()->set ("Verbosity",  Belos::Errors );

// Create the Pimpact::LinearSolver solver.
  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
      dtL, xv->clone(), xv->clone(), solverParams->get(), solverName );

  auto schur = Pimpact::createDivDtLinvGrad<S,O>( xv->clone(), lp_DTL );

  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
      schur, xs->clone(), xs->clone(), solverParams->get(), solverName );

  auto stockie = NOX::Pimpact::createInterface(xv->clone(),xs->clone(),lp_DTL,lp_Schur);


  bool succes = stockie->applyJacobian( X->getConstField(), F->getField() );

  TEST_EQUALITY( true, succes);
}


TEUCHOS_UNIT_TEST( NOXPimpact_Interface, applyJacobianInverse ) {

  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField< MSF> BMSF;
  typedef Pimpact::MultiField< MVF> BMVF;
  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
  typedef NOX::Pimpact::Vector<CF> NV;
  typedef NOX::Pimpact::Interface Interface;

  auto fS  = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::Zero2DFlow, fS, iIS, fIS );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

  auto x  = Pimpact::createCompoundField( xv, xs );

  Teuchos::RCP<NV> X = Teuchos::rcp(new NV(x) );
  Teuchos::RCP<NV> F = Teuchos::rcp_dynamic_cast<NV>( X->clone() );

  auto dtL = Pimpact::createDtL<S,O>(1.,0.,1.);

  // Make an empty new parameter list.
  auto solverName = "GMRES";
  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
  solverParams->get()->set ("Verbosity",  Belos::Errors );

// Create the Pimpact::LinearSolver solver.
  auto lp_DTL = Pimpact::createLinearProblem<Interface::BVF>(
      dtL, xv->clone(), xv->clone(), solverParams->get(), solverName );

  auto schur = Pimpact::createDivDtLinvGrad<S,O>( xv->clone(), lp_DTL );

  auto lp_Schur = Pimpact::createLinearProblem<Interface::BSF>(
      schur, xs->clone(), xs->clone(), solverParams->get(), solverName );

  auto stockie = NOX::Pimpact::createInterface(xv->clone(),xs->clone(),lp_DTL,lp_Schur);


//  NOX::Abstract::Group::ReturnType succes = stockie->applyJacobianInverse(Teuchos::parameterList(), X->getConstField(), F->getField() );

//  TEST_EQUALITY( NOX::Abstract::Group::Ok, succes);
}


TEUCHOS_UNIT_TEST( NOXPimpact_SimpleNonlinear, computeF ) {

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

  auto op = Pimpact::createOperatorBaseMV<MVF,OP>();
//  auto op = Pimpact::createOperatorMV<O(void);

  for( int i=0; i<10; ++i ) {
    x->getFieldPtr(i)->initField(Pimpact::Circle2D );
  }
//  x->random();
x->getFieldPtr(0)->write();

  op->apply( *x, *y);

  y->getFieldPtr(0)->write(99);

  auto lp = Pimpact::createLinearProblem<MVF>(
        op, x->clone(), y->clone(), Pimpact::createLinSolverParameter("GMRES",1.e-12)->get(), "GMRES");
  auto inter = NOX::Pimpact::createSimpleNonlinear( y, op, lp );

//  bool succes = inter->computeF( *x, *f );

//  TEST_EQUALITY( true, succes);
}

} // namespace

