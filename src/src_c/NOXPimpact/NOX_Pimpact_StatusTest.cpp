#include "NOX_Pimpact_StatusTest.hpp"




Teuchos::RCP< NOX::StatusTest::Generic >
NOX::Pimpact::createStatusTest( int maxI, double tolF, double tolUpdate ) {

  Teuchos::ParameterList stl;
  stl.set( "Test Type", "Combo" );
  stl.set( "Combo Type", "AND" );
  stl.set( "Number of Tests", 2 );
  Teuchos::ParameterList& conv = stl.sublist( "Test 0" );
  Teuchos::ParameterList& nstep = stl.sublist( "Test 1" );
  //  Teuchos::ParameterList& divergence = stl.sublist( "Test 3" );
  //  Teuchos::ParameterList& fv = stl.sublist( "Test 4" );
  conv.set( "Test Type", "Combo" );
  conv.set( "Combo Type", "OR" );
  conv.set( "Number of Tests", 4 );
  Teuchos::ParameterList& normF = conv.sublist( "Test 0" );
  Teuchos::ParameterList& normUpdate = conv.sublist( "Test 1" );
  Teuchos::ParameterList& maxiters = conv.sublist( "Test 2" );
  Teuchos::ParameterList& stagnation = conv.sublist( "Test 3" );
  //  Teuchos::ParameterList& normWRMS = conv.sublist( "Test 2" );
  //  Teuchos::ParameterList& userDefined = conv.sublist( "Test 3" );
  normF.set( "Test Type", "NormF" );
  normF.set( "Tolerance", tolF );
  normF.set( "Norm Type", "Two Norm" );
  normF.set( "Scale Type", "Scaled" );
  //  normWRMS.set( "Test Type", "NormWRMS" );
  //  normWRMS.set( "Absolute Tolerance", 1.0e-8 );
  //  normWRMS.set( "Relative Tolerance", 1.0e-5 );
  //  normWRMS.set( "Tolerance", 1.0 );
  //  normWRMS.set( "BDF Multiplier", 1.0 );
  //  normWRMS.set( "Alpha", 1.0 );
  //  normWRMS.set( "Beta", 0.5 );
  normUpdate.set( "Test Type", "NormUpdate" );
  normUpdate.set( "Tolerance", tolUpdate );
  normUpdate.set( "Norm Type", "Two Norm" );
  normUpdate.set( "Scale Type", "Scaled" );
  //  userDefined.set("Test Type", "User Defined");
  //  Teuchos::RCP<NOX::StatusTest::Generic> myTest =
  //      Teuchos::rcp(new MyTest(1.0e-3));
  //  userDefined.set("User Status Test", myTest);
  //  fv.set("Test Type", "FiniteValue");
  //  fv.set("Vector Type", "F Vector");
  //  fv.set("Norm Type", "Two Norm");
  //  divergence.set("Test Type", "Divergence");
  //  divergence.set("Tolerance", 1.0e+20);
  //  divergence.set("Consecutive Iterations", 3);
  stagnation.set( "Test Type", "Stagnation" );
  stagnation.set( "Tolerance", 1.0 );
  stagnation.set( "Consecutive Iterations", 10 );
  maxiters.set( "Test Type", "MaxIters" );
  maxiters.set( "Maximum Iterations", maxI );
  nstep.set( "Test Type", "NStep" );
  nstep.set( "N", 2 );
  //  Teuchos::RCP<NOX::StatusTest::Generic>
  Teuchos::RCP< NOX::StatusTest::Generic > status_tests =
      NOX::StatusTest::buildStatusTests( stl, NOX::Utils() );

  return( status_tests );

} // end of createStatusTest



Teuchos::RCP<Teuchos::ParameterList>
NOX::Pimpact::createNOXSolverParameter(
    const std::string& solverName,
    const std::string& lineSearchName ) {

  auto solverParametersPtr = Teuchos::parameterList( solverName );
  solverParametersPtr->set("Nonlinear Solver", "Line Search Based");

  // Create the directions parameters sublist
  Teuchos::ParameterList&  sl = solverParametersPtr->sublist("Direction");
  sl.set( "Method", solverName );

  if( solverName=="NonlinearCG" ) {
    Teuchos::ParameterList&  sll = sl.sublist("Nonlinear CG");
    sll.set( "Precondition", "On" );
    //    sll.set( "Restart Frequency", 10  );
  }


  Teuchos::ParameterList& lineSearchParameters = solverParametersPtr->sublist("Line Search");
  lineSearchParameters.set( "Method", lineSearchName );
  if( lineSearchName=="Backtrack" ) {
    lineSearchParameters.sublist("Backtrack").set( "Recovery Step", 1.e-6 );
  }
  if( lineSearchName=="Polynomial" ) {
    Teuchos::ParameterList& sll = lineSearchParameters.sublist("Polynomial");
    sll.set( "Interpolation Type", "Quadratic3" );
    sll.set( "Recovery Step Type", "Last Computed Step" );
  }
  return( solverParametersPtr );

} // end of createNOXSolverParameter
