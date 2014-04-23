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

#include "Pimpact_LinSolverParameter.hpp"

#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_LinearStokes.hpp"
#include "NOX_Pimpact_Group.hpp"

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




TEUCHOS_UNIT_TEST( NOXPimpact_Group, createGroup ) {
  init_impact(0,0);

//  typedef double S;
//  typedef int O;
//  typedef Pimpact::ScalarField<S,O> SF;
//  typedef Pimpact::VectorField<S,O> VF;
//  typedef Pimpact::ModeField<SF> MSF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField< MSF> BMSF;
//  typedef Pimpact::MultiField< MVF> BMVF;
//  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
//  typedef NOX::Pimpact::Vector<CF> NV;
//  typedef NOX::Pimpact::Interface<> Interface;
//
//  auto fS  = Pimpact::createFieldSpace<O>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();
//
//  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::Zero2DFlow, fS, iIS, fIS );
//
//  auto xs = Pimpact::createInitMSF<S,O>( fS );
//
//  auto x  = Pimpact::createCompoundField( xv, xs );
//
//  auto dtL = Pimpact::createMultiOperatorBase<BMVF,Pimpact::DtL<S,O> >( Pimpact::createDtL<S,O>(1.,0.,1.) );
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
//  auto schur = Pimpact::createMultiOperatorBase< BMSF, Pimpact::DivDtLinvGrad<S,O> >(
//      Pimpact::createDivDtLinvGrad<S,O>( xv->clone(), lp_DTL ) );
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
}


//TEUCHOS_UNIT_TEST( NOXPimpact_Group, computeF ) {
//
//  typedef double S;
//  typedef int O;
//  typedef Pimpact::ScalarField<S,O> SF;
//  typedef Pimpact::VectorField<S,O> VF;
//  typedef Pimpact::ModeField<SF> MSF;
//  typedef Pimpact::ModeField<VF> MVF;
//  typedef Pimpact::MultiField< MSF> BMSF;
//  typedef Pimpact::MultiField< MVF> BMVF;
//  typedef Pimpact::CompoundField< BMVF, BMSF > CF;
//  typedef NOX::Pimpact::Vector<CF> NV;
//  typedef NOX::Pimpact::LinearStokes Interface;
//
//  auto fS  = Pimpact::createFieldSpace<O>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();
//
//  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::Zero2DFlow, fS, iIS, fIS );
//
//  auto xs = Pimpact::createInitMSF<S,O>( fS );
//
//  auto x  = Pimpact::createCompoundField( xv, xs );
//
//  auto dtL = Pimpact::createMultiOperatorBase<BMVF,Pimpact::DtL<S,O> >( Pimpact::createDtL<S,O>(1.,0.,1.) );
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
//  auto schur = Pimpact::createMultiOperatorBase< BMSF, Pimpact::DivDtLinvGrad<S,O> >(
//      Pimpact::createDivDtLinvGrad<S,O>( xv->clone(), lp_DTL ) );
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


} // namespace
