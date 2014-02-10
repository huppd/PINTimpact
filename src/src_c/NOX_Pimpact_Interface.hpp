#pragma once
#ifndef NOX_PIMPACT_INTERFACE_HPP
#define NOX_PIMPACT_INTERFACE_HPP

#include "NOX_Common.H"

// Forward declarations
//class Epetra_Vector;

namespace NOX {
namespace Pimpact {

/**
 * \brief Provides a set of interfaces for users to provide information about the nonlinear problem to NOX.
 * Contains interfaces for the user to supply (1) the evaluation of the nonlinear equations, (2) the Jacobian, and (3) any preconditioning if required.
 */

/**
* \brief Supplies NOX with the set nonlinear equations.
* This is the minimum required information to solve a nonlinear
* problem using the NOX::Epetra objects for the linear algebra
*     implementation.  Used by NOX::Epetra::Group to provide a link
*     to the external code for residual fills.
*/
class Interface_Stokes {

public:

  typedef ... Field;
  typedef ... Operator;
  typedef ... LinearProblem;

  /// Constructor
  Interface() {};

  /// Destructor
  ~Required() {};

  /// Compute the function, F, given the specified input vector x. Returns true if computation was successful.
  bool computeF(const Field& x, Field& F ) ;

  /// Compute the Jacobian Operator, given the specified input vector x. Returns true if computation was successful.
  bool computeJacobian(const Field& x, Operator& Jay ) ;

}; // end of class Interface

} // end of namespace Pimpact
} // end of namespace NOX

#endif // end of #ifndef NOX_EPETRA_INTERFACE_REQUIRED_HPP

