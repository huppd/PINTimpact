#pragma once
#ifndef NOX_PIMPACT_SIMPLENONLINEAR_HPP
#define NOX_PIMPACT_SIMPLENONLINEAR_HPP

#include "NOX_Common.H"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_ModeField.hpp"
//#include "Pimpact_CompoundField.hpp"

#include "Pimpact_MultiField.hpp"

#include "Pimpact_Nonlinear.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorMV.hpp"
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
class SimpleNonlinear {

public:

  typedef double S;
  typedef int O;
  typedef ::Pimpact::VectorField<S,O> VF;
//  typedef ::Pimpact::MultiField<VF> BVF;
  typedef ::Pimpact::MultiField<VF> Field;

  typedef ::Pimpact::OperatorMV< ::Pimpact::Nonlinear<S,O> >  Op;
//  typedef ::Pimpact::OperatorMV< ::Pimpact::Div_DtLinv_Grad<S,O> >  Schur;
  typedef ::Pimpact::OperatorBase<Field> BOp;

//  typedef ::Pimpact::LinearProblem<S,BVF,BOP> LP_DTL;

//  typedef BVF Field;
  typedef NOX::Pimpact::Vector<Field> Vector;

protected:

  Teuchos::RCP<Field> fu_;
  Teuchos::RCP<BOp> nonlinear_;

public:

  /// Constructor
  SimpleNonlinear():
    fu_(Teuchos::null),
    nonlinear_(Teuchos::null) {};

  SimpleNonlinear(
      Teuchos::RCP<Field> fu,
      Teuchos::RCP<BOp> nonlinear ):
        fu_( fu ),
        nonlinear_(nonlinear) {};

  /// Destructor
  ~SimpleNonlinear() {};

  /// Compute the function, F, given the specified input vector x. Returns true if computation was successful.
  bool computeF(const Field& x, Field& f ) {
    nonlinear_->apply( x, f );
    f.add( 1., f, 1., *fu_ );
    return( true );
  }

  /// Compute the Jacobian Operator, given the specified input vector x. Returns true if computation was successful.
  bool computeJacobian( const Field& x ) { return( true ); }

  bool applyJacobian( const Field& x, Field& y, Belos::ETrans type=Belos::NOTRANS ) {
    return( false );
//    return( computeF( x, y ) );
  }


  bool applyJacobianInverse( /* param, */ const Field& x, Field& y ) {
    return( false );
  }

}; // end of class SimpleNonlinear


Teuchos::RCP<SimpleNonlinear> createSimpleNonlinear(
    Teuchos::RCP<SimpleNonlinear::Field> fu,
    Teuchos::RCP<SimpleNonlinear::BOp> op ) {
  return(
      Teuchos::rcp( new SimpleNonlinear(fu,op) )
      );
}


} // end of namespace Pimpact
} // end of namespace NOX


#endif // end of #ifndef NOX_PIMPACT_SIMPLENONLINEAR_HPP
