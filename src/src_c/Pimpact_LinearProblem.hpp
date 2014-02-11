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
//#include"BelosPimpactAdapter.hpp"

namespace Pimpact {

template< class Scalar, class MV, class OP>
class LinearProblem {
	 Teuchos::RCP< Belos::SolverManager<Scalar, MV, OP> >	solver_;
	 Teuchos::RCP< Belos::LinearProblem<Scalar, MV, OP> > problem_;
public:
	 LinearProblem(
			 Teuchos::RCP< Belos::SolverManager<Scalar, MV, OP> >	solver,
			 Teuchos::RCP< Belos::LinearProblem<Scalar, MV, OP> > problem ):
				 solver_(solver),problem_(problem) {};

	 void apply( const Teuchos::RCP<const MV>& x, const Teuchos::RCP<MV> & y ) {
	   problem_->applyOp( *x, *y );
	 }

	 Belos::ReturnType solve( const Teuchos::RCP<MV>& x, const Teuchos::RCP<const MV>& rhs) {
		 problem_->setProblem(x,rhs);
		 return( solver_->solve() );
	 }

	 const Teuchos::RCP<const Belos::SolverManager<Scalar, MV, OP> > getSolver() const {
		 return( solver_ );
	 }
	 const Teuchos::RCP<const Belos::LinearProblem<Scalar, MV, OP> > getProblem() const {
		 return( problem_ );
	 }

	 void setParameters( const Teuchos::RCP< Teuchos::ParameterList > &params ) {
		 solver_->setParameters( params );
	 }
};


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
