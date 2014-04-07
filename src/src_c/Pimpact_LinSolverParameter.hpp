#pragma once
#ifndef PIMPACT_LINSOLVERPARAMETER_HPP
#define PIMPACT_LINSOLVERPARAMETER_HPP

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "BelosTypes.hpp"



namespace Pimpact {


Teuchos::RCP< Teuchos::ParameterList > createLinSolverParameter(
    const std::string& solver_name="GMRES",
    double tol=1.e-6,
    int outfreq=100 ) {

  auto parameter_ = Teuchos::parameterList(solver_name);
  int verbo;
//  if(verbosity)
    verbo =
				Belos::Errors +
				Belos::Warnings +
				Belos::IterationDetails +
				Belos::OrthoDetails +
				Belos::FinalSummary +
				Belos::TimingDetails +
				Belos::StatusTestDetails +
				Belos::Debug ;
//  else
//    verbo =
//				Belos::Errors +
//				Belos::Warnings +
//				Belos::FinalSummary +
//				Belos::TimingDetails +
//				Belos::Debug ;

	int style = Belos::Brief;

	if( solver_name=="CG" ) {
	  parameter_->set( "Maximum Iterations",        500    );
		parameter_->set( "Verbosity",									verbo	 );
		parameter_->set( "Output Style",							style	 );
		parameter_->set( "Convergence Tolerance",			tol		 );
	}
	else if( solver_name=="GMRES" ) {
		parameter_->set( "Num Blocks",								300	  );
    parameter_->set( "Maximum Iterations",        1000  );
		parameter_->set( "Maximum Restarts",					20	  );

		parameter_->set( "Convergence Tolerance",			tol           );
		parameter_->set( "Implicit Residual Scaling",	"Norm of RHS" );
		parameter_->set( "Explicit Residual Scaling",	"Norm of RHS" );
		parameter_->set( "Deflation Quorum",					-1		        );
		parameter_->set( "Orthogonalization",					"DGKS"        );

		parameter_->set( "Verbosity",									verbo	 );
		parameter_->set( "Output Frequency",					outfreq);
		parameter_->set( "Output Style",							style	 );
		parameter_->set( "Timer Label",								"Belos");
//	parameter_->set( "Show Maximum Residual Norm Only", int(0)		 );
	}
	else if( solver_name=="GCRODR" ) {
		parameter_->set( "Num Blocks",								25		 );
    parameter_->set( "Maximum Iterations",        500    );
		parameter_->set( "Maximum Restarts",					20		 );

		parameter_->set( "Num Recycled Blocks",				5		   );

		parameter_->set( "Convergence Tolerance",			tol    );
//			parameter_->set( "Orthogonalization",					"DGKS" );
//			parameter_->set( "Implicit Residual Scaling",	"None" );
//			parameter_->set( "Explicit Residual Scaling",	"None" );

		parameter_->set( "Output Frequency",					outfreq);
		parameter_->set( "Verbosity",									verbo	 );
		parameter_->set( "Output Style",							style	 );
		parameter_->set( "Timer Label",								"Belos");
	}
	else if( solver_name=="LSQR" ) {
		parameter_->set( "Maximum Iterations",				1000	 );
		parameter_->set( "Block Size",								1			 );
		parameter_->set( "Condition Limit",						1.e20	 );
		parameter_->set( "Term Iter Max",							20		 );
		parameter_->set( "Rel RHS Err",								1.e-16 );
		parameter_->set( "Rel Mat Err",								tol    );
		parameter_->set( "Orthogonalization",					"DGKS" );
		parameter_->set( "Verbosity",									verbo	 );
		parameter_->set( "Output Frequency",					outfreq);
		parameter_->set( "Output Style",							0    	 );
		parameter_->set( "Lambda",										1.e-12 );
		parameter_->set( "Timer Label",								"Belos");
	}
	else if( solver_name=="TFQMR" ) {
//			parameter_->set( "Maximum Iterations",				1000	 );
//			parameter_->set( "Convergence Tolerance",			tol    );
//			parameter_->set( "Verbosity",									verbo	 );
//			parameter_->set( "Output Style",							style	 );
//			parameter_->set( "Output Frequency",					outfreq);
//			parameter_->set( "Timer Label",								"Belos");
	}
	else {
		std::cout << "!!!Warning!!! solver_name:\t" << solver_name << "\tnot kown in Pimpact!";
	}

	return( parameter_ );

} // end of createLinSolverParamter




} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_LINSOLVERPARAMETER_HPP
