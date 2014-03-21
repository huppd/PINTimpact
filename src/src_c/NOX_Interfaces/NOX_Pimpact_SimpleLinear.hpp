#pragma once
#ifndef NOX_PIMPACT_SIMPLELINEAR_HPP
#define NOX_PIMPACT_SIMPLELINEAR_HPP

#include "NOX_Common.H"
#include "Teuchos_RCP.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_ModeField.hpp"
//#include "Pimpact_CompoundField.hpp"

#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_LinearProblem.hpp"




namespace NOX {
namespace Pimpact {

/**
 * \brief Provides a set of interfaces for users to provide information about the nonlinear problem to NOX.
 * Contains interfaces for the user to supply
 * (1) the evaluation of the nonlinear equations,
 * (2) the Jacobian, and
 * (3) any preconditioning if required.
 */

/// \brief Supplies NOX with the set nonlinear equations.
///
/// This is the minimum required information to solve a nonlinear
/// problem using the NOX::Epetra objects for the linear algebra
/// implementation.  Used by NOX::Epetra::Group to provide a link
/// to the external code for residual fills.
/// \deprecated
class SimpleLinear {

public:

  typedef double S;
  typedef int O;
  typedef ::Pimpact::VectorField<S,O> VF;
//  typedef ::Pimpact::MultiField<VF> BVF;
  typedef ::Pimpact::MultiField<VF> Field;

//  typedef ::Pimpact::OperatorMV< ::Pimpact::Nonlinear<S,O> >  Op;
//  typedef ::Pimpact::OperatorMV< ::Pimpact::Div_DtLinv_Grad<S,O> >  Schur;
  typedef ::Pimpact::OperatorBase<Field> Op;

  typedef ::Pimpact::LinearProblem<Field> JOp;

//  typedef BVF Field;
  typedef NOX::Pimpact::Vector<Field> Vector;

protected:

  Teuchos::RCP<Field> fu_;
  Teuchos::RCP<Op>   op_;
  Teuchos::RCP<JOp>    jop_;

public:

  /// Constructor
  SimpleLinear():
    fu_(Teuchos::null),
    op_(Teuchos::null),
    jop_(Teuchos::null) {};

  SimpleLinear(
      Teuchos::RCP<Field> fu,
      Teuchos::RCP<Op> op,
      Teuchos::RCP<JOp> lp ):
        fu_( fu ),
        op_(op),
        jop_(lp) {};

  /// Destructor
  ~SimpleLinear() {};

  /// Compute the function, F, given the specified input vector x. Returns true if computation was successful.
  NOX::Abstract::Group::ReturnType computeF(const Field& x, Field& f ) {
//    Teuchos::rcp_dynamic_cast<
//      ::Pimpact::OperatorPimpl<Field,
//      ::Pimpact::Nonlinear<typename Field::Scalar,typename Field::Ordinal> > >(
//          op_, true )->getOperatorPtr()->setU(Teuchos::null);

//    Teuchos::rcp_dynamic_cast< ::Pimpact::Nonlinear<double,int> >( nonlinear_, true )->setU( x.GetConstVec(0).clone()) ;
    op_->apply( x, f );
    f.add( 1., f, 1., *fu_ );
    return( NOX::Abstract::Group::Ok );
  }


  /// Compute the Jacobian Operator, given the specified input vector x. Returns true if computation was successful.
  NOX::Abstract::Group::ReturnType computeJacobian( const Field& x ) {
//    auto opJ = jop_->getProblem()->getOperator();
//    auto opJ2 = Teuchos::rcp_const_cast<BOp>(opJ);
//    auto opJ3 = Teuchos::rcp_dynamic_cast<
//        ::Pimpact::OperatorPimpl<Field,::Pimpact::NonlinearJacobian<typename Field::Scalar,typename Field::Ordinal> > >( opJ2 );
//    opJ3->getOperatorPtr()->setU( x.GetConstVec(0).clone() );
//    Teuchos::rcp_dynamic_cast<
//      ::Pimpact::OperatorPimpl<Field,
//      ::Pimpact::NonlinearJacobian<typename Field::Scalar,typename Field::Ordinal> > >(
//          Teuchos::rcp_const_cast< ::Pimpact::OperatorBase<Field> >
//          lp_->getProblem()->getOperator()
//          , true )->getOperatorPtr()->setU( x.GetConstVec(0).clone() );
//    Teuchos::rcp_dynamic_cast< ::Pimpact::Nonlinear<double,int> >( nonlinear_, true )->setU( x.GetConstVec(0).clone()) ;
    return( NOX::Abstract::Group::Ok );
  }

  NOX::Abstract::Group::ReturnType applyJacobian( const Field& x, Field& y, Belos::ETrans type=Belos::NOTRANS ) {
//    return( NOX::Abstract::Group::NotDefined );
    jop_->getProblem()->getOperator()->apply(x,y);
    y.scale(-1.);

    return( NOX::Abstract::Group::Ok );
//    return( computeF( x, y ) );
  }


  NOX::Abstract::Group::ReturnType applyJacobianInverse( Teuchos::ParameterList &params, const Field& x, Field& y ) {
    jop_->solve( Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x) );

    return( NOX::Abstract::Group::Ok );
  }


  NOX::Abstract::Group::ReturnType applyPreconditioner( const Field& x, Field& y ) {

    jop_->solve( Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x) );
    y.scale( -1. );

    return( NOX::Abstract::Group::Ok );
  }

}; // end of class SimpleNonlinear


Teuchos::RCP<SimpleLinear> createSimpleLinear(
    Teuchos::RCP<SimpleLinear::Field> fu,
    Teuchos::RCP<SimpleLinear::Op> op,
    Teuchos::RCP<SimpleLinear::JOp> lp ) {
  return(
      Teuchos::rcp( new SimpleLinear(fu,op,lp) )
      );
}


} // end of namespace Pimpact
} // end of namespace NOX


#endif // end of #ifndef NOX_PIMPACT_SIMPLELINEAR_HPP
