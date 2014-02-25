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

#include "Pimpact_LinSolverParameter.hpp"

#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_SimpleNonlinear.hpp"
#include "NOX_Pimpact_Group.hpp"

#include "NOX.H"

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
  int rank = init_impact(0,0);

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

  auto xv = Pimpact::createInitMVF<S,O>(Pimpact::Streaming2DFlow, fS, iIS, fIS );
  xv->init( 0. );

  auto xs = Pimpact::createInitMSF<S,O>( fS );

  auto x  = Pimpact::createCompoundField( xv, xs );

  auto dtL = Pimpact::createDtL<S,O>(1.,0.,1.);

  // Make an empty new parameter list.
  auto solverName = "GMRES";
  auto solverParams = Pimpact::createLinSolverParameter( solverName, 1.e-1 );
  solverParams->get()->set ("Verbosity",  Belos::Errors );

// Create the Pimpact::LinearSolver solver.
  auto lp_DTL = Pimpact::createLinearProblem<S,Interface::BVF,Interface::DTL>(
      dtL, xv->clone(), xv->clone(), solverParams->get(), solverName );

  auto schur = Pimpact::createDivDtLinvGrad<S,O>( xv->clone(), lp_DTL );

  auto lp_Schur = Pimpact::createLinearProblem<S,Interface::BSF,Interface::Schur>(
      schur, xs->clone(), xs->clone(), solverParams->get(), solverName );

//  xv->random();
//  xs->random();
  auto stockie = NOX::Pimpact::createInterface( xv->clone(),xs->clone(),lp_DTL,lp_Schur);

  Teuchos::RCP<NV> nx = Teuchos::rcp(new NV(x) );

  auto bla = Teuchos::parameterList();

  auto group = NOX::Pimpact::createGroup<NOX::Pimpact::Interface>( bla, stockie, nx );

  // Set up the status tests
  Teuchos::RCP<NOX::StatusTest::NormF> statusTestNormF =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-4));
  Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestMaxIters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(10));
  Teuchos::RCP<NOX::StatusTest::Combo> statusTestsCombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                                            statusTestNormF,
                                            statusTestMaxIters));

  // Create the list of solver parameters
  Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr =
    Teuchos::rcp(new Teuchos::ParameterList);

  // Select the solver (this is the default)
  solverParametersPtr->set("Nonlinear Solver", "Line Search Based");

  // Create the directions parameters sublist
  solverParametersPtr->sublist("Direction").set("Method","NonlinearCG");
  solverParametersPtr->sublist("Direction").sublist("NonlinearCG").set( "Restart Frequency", 20 );

  // Create the line search parameters sublist
  solverParametersPtr->sublist("Line Search").set("Method","Polynomial");
//  Teuchos::ParameterList& lineSearchParameters = solverParametersPtr->sublist("Line Search");

  // Set the line search method
//  lineSearchParameters.set("Method","More'-Thuente");


  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver( group, statusTestsCombo, solverParametersPtr);

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
//  grp.printSolution();
//  group->getXPtr()->write();

}

TEUCHOS_UNIT_TEST( NOXPimpact_Group, SimpleNonlinear ) {
  typedef double S;
  typedef int O;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::MultiField<VF> MVF;
  typedef Pimpact::Nonlinear<S,O>  OP;
  typedef Pimpact::OperatorBase<MVF>  BOP;

  typedef NOX::Pimpact::SimpleNonlinear Interface;
  typedef NOX::Pimpact::Vector<typename Interface::Field> NV;

  int rank = 0;

  auto fS = Pimpact::createFieldSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto vel = Pimpact::createVectorField<S,O>(fS,iIS,fIS);


  auto x = Pimpact::createMultiField<VF>( *vel->clone(), 1 );
  auto f = Pimpact::createMultiField<VF>( *vel->clone(), 1 );


//  auto op = Pimpact::createOperatorMV<OP>();
  auto op = Pimpact::createOperatorBase<MVF,OP>();

  x->GetFieldPtr(0)->initField( Pimpact::Circle2D );

  f->GetFieldPtr(0)->initField( Pimpact::ZeroProf );
//  f->init(1.);
//  f->GetFieldPtr(0)->initField( Pimpact::Circle );

  x->GetFieldPtr(0)->write();
//  x->random();

  op->apply( *x, *f);
  f->scale(-1.);
  x->random();
  x->scale(0.1);

//  x->GetFieldPtr(0)->initField( Pimpact::ZeroProf );

//  f->GetFieldPtr(0)->write(99);



  auto inter = NOX::Pimpact::createSimpleNonlinear( f, op );

  Teuchos::RCP<NV> nx = Teuchos::rcp(new NV(x) );

  auto bla = Teuchos::parameterList();

  auto group = NOX::Pimpact::createGroup<Interface>( bla, inter, nx );

  // Set up the status tests
  Teuchos::RCP<NOX::StatusTest::NormF> statusTestNormF =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-4));
  Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestMaxIters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters( 10 ));
  Teuchos::RCP<NOX::StatusTest::Combo> statusTestsCombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                                            statusTestNormF,
                                            statusTestMaxIters));

  // Create the list of solver parameters
  Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr =
    Teuchos::rcp(new Teuchos::ParameterList);

  // Select the solver (this is the default)
  solverParametersPtr->set("Nonlinear Solver", "Line Search Based");

  // Create the directions parameters sublist
  solverParametersPtr->sublist("Direction").set("Method","NonlinearCG");
//  solverParametersPtr->sublist("Direction").sublist("NonlinearCG").set( "Restart Frequency", nx->length()/10 );
  solverParametersPtr->sublist("Direction").sublist("NonlinearCG").set( "Restart Frequency", 1 );

  // Create the line search parameters sublist
  solverParametersPtr->sublist("Line Search").set("Method","Polynomial");
//  Teuchos::ParameterList& lineSearchParameters = solverParametersPtr->sublist("Line Search");

  // Set the line search method
//  lineSearchParameters.set("Method","More'-Thuente");


  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver( group, statusTestsCombo, solverParametersPtr);

  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver->solve();

  // Print the parameter list
  if(rank==0) std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
  if(rank==0) std::cout << "\n" << status << "\n";
  if(rank==0) solver->getList().print(std::cout);

  // Get the answer
  *group = solver->getSolutionGroup();

  Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->write(99);

  // Print the answer
  if(rank==0) std::cout << "\n" << "-- Final Solution From Solver --" << "\n";
//  grp.printSolution();
//  group->getXPtr()->write();

}



} // namespace
