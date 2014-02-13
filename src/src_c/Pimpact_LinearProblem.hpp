#pragma once
#ifndef PIMPACT_LINEARPROBLEM_HPP
#define PIMPACT_LINEARPROBLEM_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
//#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"


//#include"Pimpact_MultiField.hpp"
//#include"Pimpact_VectorField.hpp"
//#include"Pimpact_ScalarField.hpp"
//#include"Pimpact_OperatorMV.hpp"
//#include"Pimpact_Operator.hpp"
#include"BelosPimpactAdapter.hpp"



namespace Pimpact {



template< class Scalar, class MV, class OP>
class LinearProblem {

	 Teuchos::RCP< Belos::SolverManager<Scalar, MV, OP> >	solver_;
	 Teuchos::RCP< Belos::LinearProblem<Scalar, MV, OP> > problem_;

public:

	 /// constructor
	 LinearProblem(
			 Teuchos::RCP< Belos::SolverManager<Scalar, MV, OP> >	solver,
			 Teuchos::RCP< Belos::LinearProblem<Scalar, MV, OP> > problem ):
				 solver_(solver),problem_(problem) {};


  /// \name base methods
  //@{

	 /// \brief applys Operator of Linear problem.
	 ///
	 ///( not wokring properly yet).
	 void apply( const Teuchos::RCP<const MV>& x, const Teuchos::RCP<MV> & y ) {
//	 void apply( const MV& x, MV> & y ) {
	   problem_->applyOp( *x, *y );
	 }

	 Belos::ReturnType solve( const Teuchos::RCP<MV>& x, const Teuchos::RCP<const MV>& rhs) {
		 problem_->setProblem(x,rhs);
		 return( solver_->solve() );
	 }

  //@}
  /// \name getter methods
  //@{

	 const Teuchos::RCP<const Belos::SolverManager<Scalar, MV, OP> > getSolver() const {
		 return( solver_ );
	 }
	 const Teuchos::RCP<const Belos::LinearProblem<Scalar, MV, OP> > getProblem() const {
		 return( problem_ );
	 }


  //@}
  /// \name setter methods
  //@{

  /// \brief Set the parameters to use when solving the linear problem.
  ///
  /// \param params [in/out] List of parameters to use when solving
  ///   the linear problem.  This list will be modified as necessary
  ///   to include default parameters that need not be provided.  If
  ///   params is null, then this method uses default parameters.
  ///
  /// \note The ParameterList returned by \c getValidParameters() has
  ///   all the parameters that the solver understands, possibly
  ///   including human-readable documentation and validators.
  void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params ) {
		 solver_->setParameters( params );
	 }


	/// \brief Set left preconditioner (\c LP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setLeftPrec(const Teuchos::RCP<const OP> &LP) {
    problem_->setLeftPrec( LP );
  }


  /// \brief Set right preconditioner (\c RP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setRightPrec(const Teuchos::RCP<const OP> &RP) {
    problem_->setRightPrec( RP );
  }


  //@}

}; // end of class LinearProblem


template<class Scalar, class MV, class OP>
Teuchos::RCP< LinearProblem<Scalar,MV,OP> > createLinearProblem(
		const Teuchos::RCP<const OP>& A,
		const Teuchos::RCP<MV>& x,
		const Teuchos::RCP<const MV>& b,
		Teuchos::RCP<Teuchos::ParameterList> param,
		const std::string& solvername="GMRES" ) {

	 Belos::SolverFactory<Scalar, MV, OP> factory;

	 Teuchos::RCP<Belos::SolverManager<Scalar, MV, OP> > solver =
	   factory.create( solvername, param );

	 Teuchos::RCP<Belos::LinearProblem<Scalar, MV, OP> > problem =
	    Teuchos::rcp( new Belos::LinearProblem<Scalar, MV, OP> (A, x, b) );

	 solver->setProblem(problem);
//	 problem->setProblem(x,b);

	 return( Teuchos::rcp( new LinearProblem<Scalar,MV,OP>(solver,problem) ) );

}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_LINEARPROBLEM_HPP
