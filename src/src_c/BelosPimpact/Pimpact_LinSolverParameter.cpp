#include "Pimpact_LinSolverParameter.hpp"




Teuchos::RCP< Teuchos::ParameterList > Pimpact::createLinSolverParameter(
		const std::string& solver_name,
		double tol,
		int outfreq,
		const Teuchos::RCP<std::ostream>& outStream,
		int maxIters,
	 	const std::string& label ) {

  auto parameter_ = Teuchos::parameterList(solver_name);
  int verbo;
  //  if(verbosity)
  verbo =
      Belos::Errors +
      Belos::Warnings +
      Belos::IterationDetails +
      Belos::OrthoDetails +
      Belos::FinalSummary +
//      Belos::TimingDetails +
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
  }
  else if( solver_name=="GMRES" ) {
    parameter_->set( "Num Blocks",         10	  );
    parameter_->set( "Maximum Restarts",   maxIters/10	);

    parameter_->set( "Implicit Residual Scaling",	"Norm of RHS" );
    parameter_->set( "Explicit Residual Scaling",	"Norm of RHS" );
    parameter_->set( "Deflation Quorum",					-1		        );
    parameter_->set( "Orthogonalization",					"DGKS"        );

  }
  else if( solver_name=="Block GMRES" ) {
    parameter_->set( "Num Blocks",         10	  );
    parameter_->set( "Maximum Restarts",   maxIters/10	);

    parameter_->set( "Block Size", 2);
		parameter_->set( "Adaptive Block Size", true );		

    parameter_->set( "Flexible Gmres", true );

    parameter_->set( "Implicit Residual Scaling",	"Norm of RHS" );
    parameter_->set( "Explicit Residual Scaling",	"Norm of RHS" );
    parameter_->set( "Orthogonalization",					"DGKS"        );
  }
  else if( solver_name=="GCRODR" ) {
    parameter_->set( "Num Blocks",         20	  );
    parameter_->set( "Maximum Restarts",   50	);

    parameter_->set( "Num Recycled Blocks",				2		 );
  }
  else if( solver_name=="LSQR" ) {
    //		parameter_->set( "Block Size",								1			 );
    //		parameter_->set( "Condition Limit",						1.e20	 );
    //		parameter_->set( "Term Iter Max",							20		 );
    //		parameter_->set( "Rel RHS Err",								1.e-16 );
    //		parameter_->set( "Rel Mat Err",								tol    );
    //		parameter_->set( "Orthogonalization",					"DGKS" );
    //		parameter_->set( "Lambda",										1.e-12 );
  }
  else if( solver_name=="TFQMR" ) {
  }
  else if( solver_name=="Fixed Point" ) {
  }
  else if( solver_name=="Pseudo Block TFQMR" ) {
    parameter_->set( "Block Size", 2);
  }
  else {
    std::cout << "!!!Warning!!! solver_name:\t" << solver_name << "\tnot known in Pimpact!";
  }
	parameter_->set( "Maximum Iterations",	  maxIters    );
	parameter_->set( "Convergence Tolerance",	tol         );
	parameter_->set( "Verbosity",			        verbo	      );
	parameter_->set( "Timer Label",	          label       );
	parameter_->set( "Output Frequency",      outfreq     );
  parameter_->set( "Output Stream",         outStream   );
	parameter_->set( "Output Style",		      style	      );

  return( parameter_ );



} // end of Pimpact::createLinSolverParamter
