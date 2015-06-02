#pragma once
#ifndef NOX_PIMPACT_STATUSTEST_HPP
#define NOX_PIMPACT_STATUSTEST_HPP

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "NOX_Utils.H"            // class data element
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Factory.H"



namespace NOX{
	namespace Pimpact {

		Teuchos::RCP<NOX::StatusTest::Generic>
			createStatusTest( int maxI=10, double tolF=1.e-6, double tolUpdate=1.e-6 );

		Teuchos::RCP<Teuchos::ParameterList>
			createNOXSolverParameter(
				const std::string& solverName = "NonlinearCG",
				const std::string& lineSearchName = "NonlinearCG" ); 
	}
}


#endif // end of #ifndef NOX_PIMPACT_STATUSTEST_HPP
