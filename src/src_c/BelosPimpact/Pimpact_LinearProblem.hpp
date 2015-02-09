#pragma once
#ifndef PIMPACT_LINEARPROBLEM_HPP
#define PIMPACT_LINEARPROBLEM_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "BelosSolverFactory.hpp"

#include "BelosPimpactAdapter.hpp"
#include "Pimpact_OperatorBase.hpp"



namespace Pimpact {


/// \tparam MultiFieldT has to be of type \c Pimpact::MultiField
template< class MultiFieldT >
class LinearProblem {

  typedef MultiFieldT MF;
  typedef typename MF::Scalar Scalar;
  typedef OperatorBase<MF> Op;

  Teuchos::RCP< Belos::SolverManager<Scalar, MF, Op> >	solver_;
  Teuchos::RCP< Belos::LinearProblem<Scalar, MF, Op> > problem_;

public:

  LinearProblem():solver_(Teuchos::null),problem_(Teuchos::null) {};
  /// constructor
  LinearProblem(
      Teuchos::RCP< Belos::SolverManager<Scalar, MF, Op> >	solver,
      Teuchos::RCP< Belos::LinearProblem<Scalar, MF, Op> > problem ):
        solver_(solver),problem_(problem) {};


  /// \name base methods
  //@{

  /// \brief applys Operator of Linear problem.
  ///
  ///( not wokring properly yet).
  void apply( const Teuchos::RCP<const MF>& x, const Teuchos::RCP<MF> & y ) {
    problem_->applyOp( *x, *y );
  }

  Belos::ReturnType solve( const Teuchos::RCP<MF>& x, const Teuchos::RCP<const MF>& rhs) {
    problem_->setProblem( x, rhs );
    return( solver_->solve() );
  }

  //@}
  /// \name getter methods
  //@{

  const Teuchos::RCP<const Belos::SolverManager<Scalar, MF, Op> > getSolver() const {
    return( solver_ );
  }
  const Teuchos::RCP<const Belos::LinearProblem<Scalar, MF, Op> > getProblem() const {
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
  void setLeftPrec(const Teuchos::RCP<const Op> &LP) {
    problem_->setLeftPrec( LP );
  }


  /// \brief Set right preconditioner (\c RP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setRightPrec(const Teuchos::RCP<const Op> &RP) {
    problem_->setRightPrec( RP );
  }


  //@}

}; // end of class LinearProblem


/// \relates LinearProblem
template<class MF>
Teuchos::RCP< LinearProblem<MF> > createLinearProblem(
    const Teuchos::RCP<const OperatorBase<MF> >& A,
    const Teuchos::RCP<MF>& x,
    const Teuchos::RCP<const MF>& b,
    Teuchos::RCP<Teuchos::ParameterList> param,
    const std::string& solvername="GMRES" ) {
  typedef typename MF::Scalar S;
  typedef OperatorBase<MF> Op;

  Belos::SolverFactory<S,MF,Op> factory;

  Teuchos::RCP<Belos::SolverManager<S,MF,Op> > solver =
      factory.create( solvername, param );

  Teuchos::RCP<Belos::LinearProblem<S,MF,Op> > problem =
      Teuchos::rcp( new Belos::LinearProblem<S,MF,Op> (A, x, b) );

  solver->setProblem(problem);

  return( Teuchos::rcp( new LinearProblem<MF>(solver,problem) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_LINEARPROBLEM_HPP