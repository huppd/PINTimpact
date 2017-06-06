#pragma once
#ifndef NOX_PIMPACT_STATUSTEST_HPP
#define NOX_PIMPACT_STATUSTEST_HPP


#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "NOX_StatusTest_Generic.H" // this order is necessary
#include "NOX_StatusTest_Factory.H"
#include "NOX_Utils.H"            // class data element




namespace NOX {
namespace Pimpact {

Teuchos::RCP< NOX::StatusTest::Generic >
createStatusTest( int maxI=10, double tolF=1.e-6, double tolUpdate=1.e-6 );

Teuchos::RCP<Teuchos::ParameterList>
createNOXSolverParameter(
  const std::string& solverName = "NonlinearCG",
  const std::string& lineSearchName = "NonlinearCG" );

} // end of namespace Pimpact
} // end of namespace NOX


#endif // end of #ifndef NOX_PIMPACT_STATUSTEST_HPP
