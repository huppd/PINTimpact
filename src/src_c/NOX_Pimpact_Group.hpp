#pragma once
#ifndef NOX_PIMPACT_GROUP_HPP
#define NOX_PIMPACT_GROUP_HPP

#include "NOX_Abstract_Group.H"   // base class
#include "NOX_Pimpact_Vector.hpp" // class data element
#include "NOX_Utils.H"          // class data element
#include "NOX_Common.H"         // class data element (std::string)
//#include "NOX_Epetra_LinearSystem.H"  // class data element
#include "NOX_SharedObjectTemplate.H"  // class data element
#include "Teuchos_RCP.hpp"  // class data element

#include "BelosTypes.hpp"

//#include "Teuchos_ParameterList.hpp"
//#include "NOX_Utils.H"
//#include "NOX_Epetra_Interface_Required.H"
//#include "Epetra_Vector.h"
//#include "Epetra_Operator.h"
//#include "AztecOO_ConditionNumber.h"

// Forward declares
//namespace NOX {
//  namespace Epetra {
//    class Scaling;
//    namespace Interface {
//      class Required;
//    }
//  }
//  namespace Parameter {
//    class List;
//  }
//}
//class Epetra_Vector;
//class Epetra_Operator;
//class Epetra_RowMatrix;
//class AztecOO;
//class AztecOOConditionNumber;
//class Ifpack_IlukGraph;
//class Ifpack_CrsRiluk;

namespace NOX {
namespace PIMPACT {

/**
 * \brief Concrete implementation of NOX::Abstract::Group for Pimpact.
 *
 *  This group is set up to use the linear algebra services provided
 *  through Pimpact with Belos for the linear %solver.
 */
template<class Interface>
class Group : public virtual NOX::Abstract::Group {

public:
  typedef Interface::Field   Field;
  typedef NOX::Pimpact::Vector<Field>  Vector;
  typedef Interface::Operator Operator;
  typedef Interface::LinearProblem LinearProblem;

protected:

  /// Printing Utilities object
  const NOX::Utils utils;

  /** @name Vectors */
  //@{
  /// Solution vector pointer.
  Teuchos::RCP<Vector> xVectorPtr;
  /// Solution vector.
  Vector& xVector;
  /// Right-hand-side vector pointer (function evaluation).
  Teuchos::RCP<Vector> RHSVectorPtr;
  /// Right-hand-side vector (function evaluation).
  Vector& RHSVector;
  /// Gradient vector pointer (steepest descent vector).
  Teuchos::RCP<Vector> gradVectorPtr;
  /// Gradient vector (steepest descent vector).
  Vector& gradVector;
  /// Newton direction vector pointer.
  Teuchos::RCP<Vector> NewtonVectorPtr;
  /// Newton direction vector.
  Vector& NewtonVector;
  /// An extra temporary vector, only allocated if needed.
  mutable Teuchos::RCP<Vector> tmpVectorPtr;

  //@}
  /**
   * \name IsValid flags
   *
   * True if the current solution is up-to-date with respect to the
   * currect xVector.
   */
  //@{

  bool isValidRHS;
  bool isValidJacobian;
  bool isValidGrad;
  bool isValidNewton;
  bool isValidNormNewtonSolveResidual;
  mutable bool isValidPreconditioner;
  mutable bool isValidSolverJacOp;
  bool isValidConditionNumber;

  //@}
  /// 2-Norm of the Newton solve residual: ||Js+f||
  double normNewtonSolveResidual;

  /// condition number of Jacobian
  double conditionNumber;

//  /// Pointer to the condition number object.
//  Teuchos::RCP<AztecOOConditionNumber> azConditionNumberPtr;

  /** @name Shared Operators */
  //@{
  /// Pointer to shared Interface
  Teuchos::RCP<
    NOX::SharedObject< Interface, NOX::Pimpact::Group>
    > sharedInterfacePtr;

  /// Reference to shared Interface
  NOX::SharedObject< Interface, NOX::Pimpact::Group>&
  sharedInterface;

  //@}

  // Internal flag to disable linear resid computation. False unless set.
  bool linearResidCompDisabled;


public:

  /**
   * \brief Constructor with NO linear system (VERY LIMITED).
   *
   * \warning: If this constructor is used, then methods that require
   * a Jacobian or preconditioning will not be available.  You will be
   * limited to simple algorithms like nonlinear-CG with no
   * preconditioning.
  */
  Group( Teuchos::ParameterList& printingParams,
      const Teuchos::RCP<Interface>& i,
      const Vector& initialGuess):
        utils(printParams),
        xVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(x.clone(DeepCopy))),
        xVector(*xVectorPtr),
        RHSVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(x.clone(ShapeCopy))),
        RHSVector(*RHSVectorPtr),
        gradVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(x.clone(ShapeCopy))),
        gradVector(*gradVectorPtr),
        NewtonVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(x.clone(ShapeCopy))),
        NewtonVector(*NewtonVectorPtr),
        normNewtonSolveResidual(0),
        conditionNumber(0.0),
        sharedInterfacePtr(i),
        sharedInterface(*sharedInterfacePtr),
        linearResidCompDisabled(false) {
    // Set all isValid flags to false
    resetIsValid();
  }


  /// Standard Constructor.
//  Group(Teuchos::ParameterList& printingParams,
//      const Teuchos::RCP<Interface>& i,
//      const Vector& initialGuess,
//      const Teuchos::RCP<Operator>& op,
//      const Teuchos::RCP<Operator>& linSys ):
//    utils( printParams ),
//    xVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(x.clone(DeepCopy))),
//    xVector(*xVectorPtr),
//    RHSVectorPtr(rcp_dynamic_cast<Vector>(x.clone(ShapeCopy))),
//    RHSVector(*RHSVectorPtr),
//    gradVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(x.clone(ShapeCopy))),
//    gradVector(*gradVectorPtr),
//    NewtonVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(x.clone(ShapeCopy))),
//    NewtonVector(*NewtonVectorPtr),
//    normNewtonSolveResidual(0),
//    conditionNumber(0.0),
//    sharedOperatorPtr(Teuchos::rcp(new NOX::SharedObject<Operator, NOX::Pimpact::Group>(op))),
//    sharedOperator(*sharedOperatorPtr),
//    sharedLinearSystemPtr(Teuchos::rcp(new NOX::SharedObject<Operator, NOX::Pimpact::Group>(linSys))),
//    sharedLinearSystem(*sharedLinearSystemPtr),
//    linearResidCompDisabled(false),
//    userInterfacePtr(i) {
//    // Set all isValid flags to false
//    resetIsValid();
//  };


  /**
   * \brief Copy constructor. If type is DeepCopy, takes ownership of valid
   * shared linear system.
   */
  Group(const NOX::Pimpact::Group& source, NOX::CopyType type = NOX::DeepCopy):
      utils(source.utils),
      xVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(source.xVector.clone(type))),
      xVector(*xVectorPtr),
      RHSVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(source.RHSVector.clone(type))),
      RHSVector(*RHSVectorPtr),
      gradVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(source.gradVector.clone(type))),
      gradVector(*gradVectorPtr),
      NewtonVectorPtr(Teuchos::rcp_dynamic_cast<Vector>(source.NewtonVector.clone(type))),
      NewtonVector(*NewtonVectorPtr),
      sharedInterfacePtr(source.sharedInterfacePtr),
      sharedInterface(*sharedInterface),
      linearResidCompDisabled(source.linearResidCompDisabled) {

    switch (type) {

    case DeepCopy:

      isValidRHS = source.isValidRHS;
      isValidJacobian = source.isValidJacobian;
      isValidGrad = source.isValidGrad;
      isValidNewton = source.isValidNewton;
      isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;
      isValidConditionNumber = source.isValidConditionNumber;
      normNewtonSolveResidual = source.normNewtonSolveResidual;
      conditionNumber = source.conditionNumber;
      isValidPreconditioner = source.isValidPreconditioner;
      isValidSolverJacOp = source.isValidSolverJacOp;

      // New copy takes ownership of the shared Jacobian for DeepCopy
      if (isValidJacobian)
        sharedLinearSystem.getObject(this);

      break;

    case ShapeCopy:
      resetIsValid();
      break;

    default:
      std::cerr << "ERROR: Invalid ConstructorType for group copy constructor." << std::endl;
      throw "NOX Error";
    }
  }


  /// Destructor.
  virtual ~Group() {}


  virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source) {
    return( operator=(dynamic_cast<const Group&> (source)) );
  }


  /// See operator=(const NOX::Abstract::Group&);
  virtual NOX::Abstract::Group& operator=(const NOX::Pimpact::Group<Interface>& source) {
    // Copy the xVector
    xVector = source.xVector;

    // Copy reference to sharedJacobian
    sharedInterfacePtr = source.sharedInterfacePtr;

    // Update the isValidVectors
    isValidRHS = source.isValidRHS;
    isValidGrad = source.isValidGrad;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;
    isValidPreconditioner = source.isValidPreconditioner;
    isValidSolverJacOp = source.isValidSolverJacOp;
    isValidConditionNumber = source.isValidConditionNumber;

    // Only copy vectors that are valid
    if (isValidRHS) {
      RHSVector = source.RHSVector;
    }

    if (isValidGrad)
      gradVector = source.gradVector;

    if (isValidNewton)
      NewtonVector = source.NewtonVector;

    if (isValidNormNewtonSolveResidual)
      normNewtonSolveResidual = source.normNewtonSolveResidual;

    // If valid, this takes ownership of the shared Jacobian
    if (isValidJacobian)
      sharedInterface.getObject(this);

    if (isValidConditionNumber)
      conditionNumber = source.conditionNumber;

    linearResidCompDisabled = source.linearResidCompDisabled;

    return( *this );
  }


  /** @name "Compute" functions. */
  //@{

  virtual void setX(const Vector& y) {
    if (isPreconditioner()) {
//      sharedLinearSystem.getObject(this)->destroyPreconditioner();
    }
    resetIsValid();
    xVector = y;
    return;
  }
  virtual void setX(const NOX::Abstract::Vector& y) {
    setX(dynamic_cast<const Vector&> (y));
    return;
  }


  virtual void computeX( const Group& grp,
      const Vector& d,
      double step ) {
    if( isPreconditioner() )
//      sharedLinearSystem.getObject(this)->destroyPreconditioner();
    resetIsValid();
    xVector.update(1.0, grp.xVector, step, d);
    return;
  }
  virtual void computeX(const NOX::Abstract::Group& grp,
      const NOX::Abstract::Vector& d,
      double step);{
    // Cast to appropriate type, then call the "native" computeX
    const Group& pimpgrp = dynamic_cast<const Group&> (grp);
    const Vector& pimpd =  dynamic_cast<const Vector&> (d);
    computeX( pimpgrp, pimpd, step );
    return;
  }


  virtual NOX::Abstract::Group::ReturnType computeF() {
    if( isF() )
      return( Abstract::Group::Ok );

    bool status = false;

    status = userInterfacePtr->computeF( xVector.getField(),
        RHSVector.getField() );

    if (status == false) {
      std::cout << "ERROR: Pimpact::Group::computeF() - fill failed!!!"
          << std::endl;
      throw "NOX Error: Fill Failed";
    }

    isValidRHS = true;

    return( Abstract::Group::Ok );
  }


  virtual NOX::Abstract::Group::ReturnType computeJacobian() {
    // Skip if the Jacobian is already valid
    if (isJacobian())
      return( Abstract::Group::Ok );

    // Fill the Jacobian
    bool status = false;

//    status = sharedLinearSystem.getObject(this)->
//        computeJacobian(xVector);
    status = sharedInterfacePtr->computeJacobian( xVector.getField() );

    if (status == false) {
      std::cout << "ERROR: NOX::Pimpact::Group::computeJacobian() - fill failed!!!"
          << std::endl;
      throw "NOX Error: Fill Failed";
    }

    // Update status of Jacobian wrt solution vector
    isValidJacobian = true;

    return( Abstract::Group::Ok );
  }


  virtual NOX::Abstract::Group::ReturnType computeGradient() {
    if( isGradient() )
      return( Abstract::Group::Ok );

    if (!isF()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::computeGradient() - RHS is out of date wrt X!" << std::endl;
      throw "NOX Error";
    }

    if (!isJacobian()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::computeGradient() - Jacobian is out of date wrt X!" << std::endl;
      throw "NOX Error";
    }

    // Compute grad = Jacobian^T * RHS.
    sharedInterface.getObject(this)->applyJacobian( RHSVector,
							       gradVector, Belos::TRANS );

    // Update state
    isValidGrad = true;

    // Return result
    return( Abstract::Group::Ok );
  }


  virtual NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params) {
    if (isNewton())
      return( Abstract::Group::Ok );

    if (!isF()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::computeNewton() - invalid RHS" << std::endl;
      throw "NOX Error";
    }

    if (!isJacobian()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::computeNewton() - invalid Jacobian" << std::endl;
      throw "NOX Error";
    }

    Abstract::Group::ReturnType status;

    // Zero out the Newton Vector
    NewtonVector.init( 0.0 );

    // Create Epetra problem for the linear solve
    status = applyJacobianInverse(params, RHSVector, NewtonVector);

    // Scale soln by -1
    NewtonVector.scale(-1.0);

    // Update state EVEN IF LINEAR SOLVE FAILED
    // We still may want to use the vector even it it just missed it's
    isValidNewton = true;

    // Compute the 2-norm of the linear solve residual ||Js+f||
    // Can be disabled, but then disallows inexact Newton methods
    if (!linearResidCompDisabled)
      computeNormNewtonSolveResidual();

    // Return solution
    return( status );
  }


  //@}
  /** \name Jacobian operations.
   *
   * Operations using the Jacobian matrix. These may not be defined in
   * matrix-free scenarios. */
  //@{

  virtual NOX::Abstract::Group::ReturnType
  applyJacobian(const Vector& input, Vector& result) const {
    // Check validity of the Jacobian
    if (!isJacobian())
      return( Abstract::Group::BadDependency );

    // Apply the Jacobian
    Belos::ReturnType status = sharedInterface.getObject()->applyJacobian( result, input );

    switch(status) {
    case Belos::Converged: return( Abstract::Group::Ok );
    case Belos::Unconverged: return( Abstract::Group::Failed );
    }
  }
  virtual NOX::Abstract::Group::ReturnType
  applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const {
    const Vector& pimpinput =  dynamic_cast<const Vector&> (input);
    Vector& pimpresult = dynamic_cast<Vector&> (result);

    return( applyJacobian( pimpinput, pimpresult) );
  }


  virtual NOX::Abstract::Group::ReturnType
  applyJacobianTranspose(const Vector& input, Vector& result) const {
    // Check validity of the Jacobian
    if( !isJacobian() )
      return( Abstract::Group::BadDependency );

    bool status = sharedInterface.getObject()->applyJacobian( input, result, Belos::Trans );

    return( status == true ? Abstract::Group::Ok : Abstract::Group::Failed );
  }
  virtual NOX::Abstract::Group::ReturnType
  applyJacobianTranspose( const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result ) const {
    const Vector& pimpinput = dynamic_cast<const Vector&> (input);
    Vector& pimpresult = dynamic_cast<Vector&> (result);
    return( applyJacobianTranspose( pimpinput, pimpresult ) );
  }


  /**
   * \brief Applies the inverse of the Jacobian matrix to the given
   * input vector and puts the answer in result.
   *
   * Computes
   * \f[ v = J^{-1} u, \f]
   * where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is
   * the result vector.
   * The "Tolerance" parameter specifies that the solution should be such that
   * \f[ \frac{\| J v - u \|_2}{\max \{ 1, \|u\|_2\} } < \mbox{Tolerance} \f]
   *
   * \return
   * <ul>
   * <li> NOX::Abstract::Group::NotDefined - Returned by default implementation
   * in NOX::Abstract::Group
   * <li> NOX::Abstract::Group::BadDependency - If \f$J\f$ has not been computed     <li> NOX::Abstract::Group::NotConverged - If the linear solve fails to satisfy the "Tolerance"
   * specified in \c params
   * <li> NOX::Abstract::Group::Failed - If the computation fails
   * <li> NOX::Abstract::Group::Ok - Otherwise
   * </ul>
   * The parameter "Tolerance" may be added/modified in the list of
   * parameters - this is the ideal solution tolerance for an iterative
   * linear solve.
   *
   * The parameter "Reuse Preconditioner" is a boolean that tells the group to turn off control of preconditioner recalculation.  This is a dangerous flag but can really speed the computations if the user knows what they are doing.  Toggling this flag is left to the user (ideally it should be done through a status test).  Defaults to false.
   */
  virtual NOX::Abstract::Group::ReturnType
  applyJacobianInverse(Teuchos::ParameterList &params, const Vector &input, Vector &result) const {
    if (!isJacobian())
      return( Abstract::Group::BadDependency );

    if (!isValidSolverJacOp) {
//      sharedLinearSystem.getObject(this)->setJacobianOperatorForSolve(sharedLinearSystem.getObject(this)->getJacobianOperator());
      isValidSolverJacOp = true;
    }

    // Compute the preconditioner
//    NOX::Epetra::LinearSystem::PreconditionerReusePolicyType precPolicy =
//        sharedLinearSystem.getObject(this)->getPreconditionerPolicy();

//    if (!isPreconditioner()) {
//      if (precPolicy == NOX::Epetra::LinearSystem::PRPT_REBUILD) {
//        sharedLinearSystem.getObject(this)->destroyPreconditioner();
//        sharedLinearSystem.getObject(this)->
//            createPreconditioner(xVector, p, false);
//        isValidPreconditioner = true;
//      }
//      else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_RECOMPUTE) {
//        sharedLinearSystem.getObject(this)->recomputePreconditioner(xVector, p);
//        isValidPreconditioner = true;
//      }
//      else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_REUSE) {
//        // Do Nothing!!!
//      }
//    }

    bool status = sharedInterface.getObject(this)->applyJacobianInverse(params, input, result);

    return( status == true ? Abstract::Group::Ok : Abstract::Group::NotConverged );
  }
  virtual NOX::Abstract::Group::ReturnType
  applyJacobianInverse(Teuchos::ParameterList &params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const {
    const Vector& pimpInput = dynamic_cast<const Vector&>(input);
    Vector& pimpResult = dynamic_cast<Vector&>(result);
    return( applyJacobianInverse(params, pimpInput, pimpResult) );
  }


  virtual NOX::Abstract::Group::ReturnType
  applyRightPreconditioning(bool useTranspose,
          Teuchos::ParameterList& params,
          const Vector& input,
          Vector& result) const {
    bool success = false;

//    if( !isPreconditioner() ) {
//      sharedLinearSystem.getObject(this)->destroyPreconditioner();
//      sharedLinearSystem.getObject(this)->
//          createPreconditioner(xVector, linearSolverParams, false);
//      isValidPreconditioner = true;
//    }
//
//    success = sharedLinearSystem.getObject()->
//        applyRightPreconditioning(useTranspose, linearSolverParams, input, result);

    if( success == true )
      return( Abstract::Group::Ok );
    else
      return( Abstract::Group::Failed );
  }
  virtual NOX::Abstract::Group::ReturnType
  applyRightPreconditioning(bool useTranspose,
          Teuchos::ParameterList& params,
          const NOX::Abstract::Vector& input,
          NOX::Abstract::Vector& result) const {
    const Vector& pimpInput = dynamic_cast<const Vector&>(input);
    Vector& pimpResult = dynamic_cast<Vector&>(result);

    return( applyRightPreconditioning(useTranspose, params, pimpInput, pimpResult) );
  }


  //@}
  /**
   * \name "Is" functions
   *
   * Checks to see if various objects have been computed. Returns true
   * if the corresponding "compute" function has been called since the
   * last update to the solution vector (via instantiation or
   * computeX).
   */
  //@{

  virtual bool isF() const { return( isValidRHS ); }
  virtual bool isJacobian() const { return( ((sharedInterface.isOwner(this)) && (isValidJacobian)) ); }
  virtual bool isGradient() const { return( isValidGrad ); }
  virtual bool isNewton() const { return( isValidNewton ); }

  /**
   * \brief Returns true if the value of the Norm of the linear model
   * for a full Newton step ||Js + f|| is valid with respect to the
   * current solution vector.
   */
  virtual bool isNormNewtonSolveResidual() const {
    return( isValidNormNewtonSolveResidual );
  }

  /**
   * \brief Returns true if an explicitly constructed preconditioner
   * exists (i.e. one that is computed and saved for further use in
   * multiple calls to applyRightPreconditioner).
   */
  virtual bool isPreconditioner() const {
    return ((sharedInterface.isOwner(this)) && (isValidPreconditioner) &&
        (sharedInterface.getObject(this)->isPreconditionerConstructed()));
  }


  /**
   * \brief Returns true if the condition number has been computed.
   */
  virtual bool isConditionNumber() const {
    return( isValidConditionNumber );
  }


  //@}
  /**
   * \name "Get" functions
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose.
   */
  //@{

  virtual const NOX::Abstract::Vector& getX() const {
    return( xVector );
  }


  virtual const NOX::Abstract::Vector& getF() const {
    if (!isF()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::getF() - invalid RHS" << std::endl;
      throw "NOX Error";
    }

    return( RHSVector );
  }


  virtual double getNormF() const {
    if (!isF()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::getNormF() - invalid RHS" << std::endl;
      throw "NOX Error";
    }

    return( RHSVectorPtr->norm() );
  }


  virtual const NOX::Abstract::Vector& getGradient() const {
    if (!isGradient()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::getGradient() - invalid gradient" << std::endl;
      throw "NOX Error";
    }

    return( gradVector );
  }


  virtual const NOX::Abstract::Vector& getNewton() const {
    if (!isNewton()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::getNewton() - invalid Newton vector" << std::endl;
      throw "NOX Error";
    }

    return( *NewtonVectorPtr );
  }


  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const {return( xVectorPtr );};
  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const {return( RHSVectorPtr );};
  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const {return( gradVectorPtr );};
  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const {return( NewtonVectorPtr );};


  /**
   * \brief Returns the 2-norm of the residual of the linear model used in the
   * Newton solve computation, ||Js+f||.  This does not account for
   * line search adjustments to the step length!
   */
  virtual NOX::Abstract::Group::ReturnType getNormLastLinearSolveResidual(double & residual) const {
    // Make sure value is not already calculated
    if (isValidNormNewtonSolveResidual) {
      residual = normNewtonSolveResidual;
      return( NOX::Abstract::Group::Ok );
    }

    // Otherwise give warning since a Newton direction has not been calculated
    // wrt this solution group
    if( utils.isPrintType(Utils::Warning) ) {
      std::cout << "ERROR: NOX::Epetra::Group::getNormLastLinearSolveResidual() - "
          << "Group has not performed a Newton solve corresponding to this "
          << "solution vector, or disableLinearSolveResidual(true) was set!" << std::endl;
    }
    return( NOX::Abstract::Group::BadDependency );
  }


  //@}

  virtual Teuchos::RCP<NOX::Abstract::Group>
  clone(CopyType type = DeepCopy) const {
    Teuchos::RCP<NOX::Abstract::Group> newgrp =
        Teuchos::rcp(new NOX::Pimpact::Group(*this, type));
    return( newgrp );
  }

//  /// Return the userInterface.
//  virtual Teuchos::RCP<Interface>
//  getRequiredInterface(); {
//    return userInterfacePtr;
//  }

  /// Return the Linear System.
  virtual Teuchos::RCP<const Interface> getInterface() const {
    return( sharedInterface.getObject(this) );
  }


  /// Return the Linear System.
  virtual Teuchos::RCP<Interface> getLinearSystem() {
    return( sharedInterface.getObject() );
  }


  /**
   *  \brief Computes the condition number of the Jacobian matrix.
   *
   *  Uses GMRES to estimate the condtion number.  The tolerance and
   *  maxIters are used to control the GMRES iterations.  Typically the
   *  solves do not have to be tight to get the estimate.
   */
//  virtual NOX::Abstract::Group::ReturnType
//  computeJacobianConditionNumber(int maxIters, double tolerance,
//         int krylovSubspaceSize=100,
//         bool printOutput=false) {
//    if (!isConditionNumber()) {
//      if (!isJacobian()) {
//        std::cerr << "ERROR: NOX::Epetra::Group::computeJacobianConditionNumber()"
//            << " - Jacobian is invalid wrt the solution." << std::endl;
//        throw "NOX Error";
//      }
//
//      if (Teuchos::is_null(azConditionNumberPtr))
//        azConditionNumberPtr = Teuchos::rcp(new AztecOOConditionNumber);
//
//        azConditionNumberPtr->
//        initialize(*(sharedLinearSystem.getObject()->getJacobianOperator()),
//            AztecOOConditionNumber::GMRES_, krylovSubspaceSize,
//            printOutput);
//
//        azConditionNumberPtr->computeConditionNumber(maxIters, tolerance);
//
//        conditionNumber = azConditionNumberPtr->getConditionNumber();
//
//        isValidConditionNumber = true;
//    }
//    return( NOX::Abstract::Group::Ok );
//  }


//  /// Returns the condition number of the Jacobian matrix.
//  virtual double getJacobianConditionNumber() const {
//    if( !isConditionNumber() ) {
//      std::cerr << "ERROR: NOX::Epetra::Group::getJacobianConditionNumber()"
//          << " - condition number has not yet been computed!" << std::endl;
//      throw "NOX Error";
//    }
//    return( conditionNumber );
//  }


  /**
   * \brief Sets option to disable linear resid computation. If disabled,
   * this saves on a MatVec per Newton but disallows inexact Newton methods
   */
  virtual void disableLinearResidualComputation(const bool disableChoice) {
    linearResidCompDisabled = disableChoice;
  }

protected:
  /// resets the isValid flags to false
  virtual void resetIsValid() {
    isValidRHS = false;
    isValidJacobian = false;
    isValidGrad = false;
    isValidNewton = false;
    isValidNormNewtonSolveResidual = false;
    isValidPreconditioner = false;
    isValidSolverJacOp = false;
    isValidConditionNumber = false;
    return;
  }


  /**
   * \brief Computes the 2-norm of the residual of the linear model used in
   * the Newton solve computation, ||Js+f||.
   */
  virtual bool computeNormNewtonSolveResidual() {
    // Make sure value is not already calculated
    if (isValidNormNewtonSolveResidual)
      return( true );

    // Make sure NewtonVector and RHSVector are valid
    // We could return false, but for now we will throw errors
    if (!isValidRHS) {
      std::cerr << "ERROR: NOX::Pimpact::Group::computeNormNewtonSolveResidual() - invalid RHS"
          << std::endl;
      throw "NOX Error";
    }
    if (!isValidNewton) {
      std::cerr << "ERROR: NOX::Pimpact::Group::computeNormNewtonSolveResidual() - invalid "
          << "Newton direction" << std::endl;
      throw "NOX Error";
    }

    // Allocate the tmpVectorPtr if not already done (deleted in ~Group)
    if (Teuchos::is_null(tmpVectorPtr)) {
      tmpVectorPtr = RHSVector.clone();
//          Teuchos::rcp(new Epetra_Vector(RHSVector.getEpetraVector()));
    }
    Vector tmpNoxVector(*tmpVectorPtr, ShapeCopy);

    sharedInterface.getObject()->applyJacobian(NewtonVector, tmpNoxVector);
    tmpNoxVector.update(1.0, RHSVector, 1.0);
    normNewtonSolveResidual = tmpNoxVector.norm();

    isValidNormNewtonSolveResidual = true;

    return( true );
  }



}; // end of class Group

} // end of namespace PIMPACT
} // end of namespace NOX

#endif // end of #ifndef NOX_PIMPACT_GROUP_HPP
