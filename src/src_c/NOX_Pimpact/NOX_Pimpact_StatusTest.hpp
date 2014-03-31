#pragma once
#ifndef NOX_PIMPACT_STATUSTEST_HPP
#define NOX_PIMPACT_STATUSTEST_HPP

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

//#include "BelosTypes.hpp"
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Factory.H"



namespace NOX{
namespace Pimpact {


Teuchos::RCP<NOX::StatusTest::Generic> createStatusTest( int maxI=10, double tolF=1.e-6, double tolUpdate=1.e-4 ) {

  Teuchos::ParameterList stl;
  stl.set( "Test Type", "Combo" );
  stl.set( "Combo Type", "OR" );
  stl.set( "Number of Tests", 2 );
  Teuchos::ParameterList& conv = stl.sublist( "Test 0" );
  Teuchos::ParameterList& maxiters = stl.sublist( "Test 1" );
//  Teuchos::ParameterList& fv = stl.sublist( "Test 2" );
//  Teuchos::ParameterList& divergence = stl.sublist( "Test 3" );
//  Teuchos::ParameterList& stagnation = stl.sublist( "Test 4" );
  conv.set( "Test Type", "Combo" );
  conv.set( "Combo Type", "AND" );
  conv.set( "Number of Tests", 2 );
  Teuchos::ParameterList& normF = conv.sublist( "Test 0" );
  Teuchos::ParameterList& normUpdate = conv.sublist( "Test 1" );
//  Teuchos::ParameterList& normWRMS = conv.sublist( "Test 2" );
//  Teuchos::ParameterList& userDefined = conv.sublist( "Test 3" );
  normF.set( "Test Type", "NormF" );
  normF.set( "Tolerance", tolF );
  normF.set( "Norm Type", "Two Norm" );
  normF.set( "Scale Type", "Unscaled" );
//  normWRMS.set( "Test Type", "NormWRMS" );
//  normWRMS.set( "Absolute Tolerance", 1.0e-8 );
//  normWRMS.set( "Relative Tolerance", 1.0e-5 );
//  normWRMS.set( "Tolerance", 1.0 );
//  normWRMS.set( "BDF Multiplier", 1.0 );
//  normWRMS.set( "Alpha", 1.0 );
//  normWRMS.set( "Beta", 0.5 );
  normUpdate.set( "Test Type", "NormUpdate" );
  normUpdate.set( "Norm Type", "One Norm" );
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
//  stagnation.set("Test Type", "Stagnation");
//  stagnation.set("Tolerance", 1.0);
//  stagnation.set("Consecutive Iterations", 5);
  maxiters.set( "Test Type", "MaxIters" );
  maxiters.set( "Maximum Iterations", maxI );
//  Teuchos::RCP<NOX::StatusTest::Generic>
  auto status_tests =
      NOX::StatusTest::buildStatusTests( stl, NOX::Utils() );

  return( status_tests );

} // end of create


} // end of namespace Pimpact
} // end of namespace NOX


#endif // end of #ifndef NOX_PIMPACT_STATUSTEST_HPP
