#pragma once
#ifndef NOX_PIMPACT_INTERFACE_HPP
#define NOX_PIMPACT_INTERFACE_HPP

#include "NOX_Common.H"
#include "NOX_Abstract_Group.H"
#include "Teuchos_RCP.hpp"


#include "Pimpact_MultiField.hpp"

#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_LinearProblem.hpp"



namespace NOX {
namespace Pimpact {

/// \brief Provides a set of interfaces for users to provide information about the nonlinear problem to NOX.
/// Contains interfaces for the user to supply
/// (1) the evaluation of the nonlinear equations,
/// (2) the Jacobian, and
/// (3) any preconditioning if required.

/// \brief Supplies NOX with the set nonlinear equations.
///
/// This is the minimum required information to solve a nonlinear
/// problem using the NOX::Epetra objects for the linear algebra
/// implementation.  Used by NOX::Epetra::Group to provide a link
/// to the external code for residual fills.
/// \tparam Field hast to be of type \c Pimpact::MultiField.
template<class F>
class Interface {

public:

  typedef F Field;
  typedef typename Field::Scalar S;
  typedef typename Field::Ordinal O;
//  typedef ::Pimpact::VectorField<S,O> VF;
////  typedef ::Pimpact::MultiField<VF> BVF;
//  typedef ::Pimpact::MultiField<VF> Field;
  typedef ::Pimpact::OperatorBase<Field> Op;
//  typedef ::Pimpact::LinearProblem<Field> JOp;
  typedef Op JOp;

//  typedef BVF Field;
  typedef NOX::Pimpact::Vector<Field> Vector;

protected:

  Teuchos::RCP<Field> fu_;
  Teuchos::RCP<Op>    op_;
  Teuchos::RCP<JOp>  jop_;
//  Teuchos::RCP<Op>  jop_;

public:

  /// Constructor
  Interface():
    fu_(Teuchos::null),
    op_(Teuchos::null),
    jop_(Teuchos::null) {};

  Interface(
      Teuchos::RCP<Field> fu,
      Teuchos::RCP<Op> op,
      Teuchos::RCP<JOp> lp ):
        fu_( fu ),
        op_(op),
        jop_(lp) {};

  /// Destructor
  ~Interface() {};


  /// Compute the function, F, given the specified input vector x.
  NOX::Abstract::Group::ReturnType computeF(const Field& x, Field& f ) {
    op_->apply( x, f );
    f.add( 1., f, -1., *fu_ );
    return( NOX::Abstract::Group::Ok );
  }


  /// \brief Compute the Jacobian Operator, given the specified input vector x.
  NOX::Abstract::Group::ReturnType computeJacobian( const Field& x ) {

//    auto prob = jop_->getProblem();
//    auto opJ = Teuchos::rcp_const_cast<Op>( prob->getOperator() );
//    opJ->assignField( x );
//    if( prob->isLeftPrec() ) {
//      auto opPrec = Teuchos::rcp_const_cast<Op>( prob->getLeftPrec() );
////      if( !opPrec.is_null() )
//        opPrec->assignField( x );
//    }
    jop_->assignField( x );

    return( NOX::Abstract::Group::Ok );
  }


  NOX::Abstract::Group::ReturnType applyJacobian( const Field& x, Field& y, Belos::ETrans type=Belos::NOTRANS ) {
//    if( Belos::NOTRANS==type ) {
//      jop_->getProblem()->getOperator()->apply(x,y);
//      return( NOX::Abstract::Group::Ok );
//    }
    return( NOX::Abstract::Group::NotDefined );
  }


  NOX::Abstract::Group::ReturnType applyJacobianInverse( Teuchos::ParameterList &params, const Field& x, Field& y ) {
//    jop_->solve( Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x) );
    jop_->apply( x, y );

    return( NOX::Abstract::Group::Ok );
  }


  NOX::Abstract::Group::ReturnType applyPreconditioner( const Field& x, Field& y ) {

//    jop_->solve( Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x) );
//    y.scale( -1. );

    return( NOX::Abstract::Group::Ok );
  }

}; // end of class Interface



/// \relates Interface
template<class Field>
Teuchos::RCP< Interface<Field> > createInterface(
    Teuchos::RCP<Field> fu,
    Teuchos::RCP<typename Interface<Field>::Op> op,
    Teuchos::RCP<typename Interface<Field>::JOp> lp ) {
  return(
      Teuchos::rcp( new Interface<Field>(fu,op,lp) )
      );
}


} // end of namespace Pimpact
} // end of namespace NOX


#endif // end of #ifndef NOX_PIMPACT_INTERFACE_HPP
