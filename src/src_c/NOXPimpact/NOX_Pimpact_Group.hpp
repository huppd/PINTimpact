#pragma once
#ifndef NOX_PIMPACT_GROUP_HPP
#define NOX_PIMPACT_GROUP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "BelosTypes.hpp"

#include "NOX_Abstract_Group.H"   // base class
#include "NOX_Common.H"           // class data element (std::string)
#include "NOX_SharedObjectTemplate.H"  // class data element
#include "NOX_Utils.H"            // class data element

#include "NOX_Pimpact_Vector.hpp" // class data element


namespace NOX {
namespace Pimpact  {

/// \brief Concrete implementation of NOX::Abstract::Group for Pimpact.
///
/// This group is set up to use the linear algebra services provided
/// through Pimpact with Belos for the linear %solver.
template<class InterfaceT>
class Group : public virtual NOX::Abstract::Group {

public:

  using FieldT = typename InterfaceT::FieldT;
  using VectorT = NOX::Pimpact::Vector<FieldT>;

protected:

  /// Printing Utilities object
  const NOX::Utils utils_;

  /// @name Vectors
  //@{
  /// Solution vector pointer.
  Teuchos::RCP<VectorT> xVectorPtr_;
  /// Solution vector.
  VectorT& xVector_;
  /// Right-hand-side vector pointer (function evaluation).
  Teuchos::RCP<VectorT> RHSVectorPtr_;
  /// Right-hand-side vector (function evaluation).
  VectorT& RHSVector;
  /// Newton direction vector pointer.
  Teuchos::RCP<VectorT> NewtonVectorPtr_;
  /// Newton direction vector.
  VectorT& NewtonVector;

  //@}
  /// \name IsValid flags
  ///
  /// True if the current solution is up-to-date with respect to the
  /// currect xVector_.
  //@{

  bool isValidRHS_;
  bool isValidJacobian_;
  bool isValidNewton;
  bool isValidNormNewtonSolveResidual;
  mutable bool isValidPreconditioner;
  mutable bool isValidSolverJacOp;

  //@}
  /// 2-Norm of the Newton solve residual: ||Js+f||
  double normNewtonSolveResidual;

  /// 2-Norm of the residual F 
  double normF_;


  /// @name Shared Operators
  //@{
  /// Pointer to shared Interface
  Teuchos::RCP< NOX::SharedObject< InterfaceT, NOX::Pimpact::Group<InterfaceT> > >
    sharedInterfacePtr_;

  /// Reference to shared Interface
  NOX::SharedObject< InterfaceT, NOX::Pimpact::Group<InterfaceT> >& sharedInterface_;

  //@}

  // Internal flag to disable linear resid computation. False unless set.
  bool linearResidCompDisabled;

  Teuchos::RCP<Teuchos::Time> noxCompF_;
  Teuchos::RCP<Teuchos::Time> noxUpdateX_;
  Teuchos::RCP<Teuchos::Time> noxCompdx_;
  Teuchos::RCP<Teuchos::Time> noxAssign_;

public:

  /// \brief Constructor with NO linear system.
  ///
  /// \warning: If this constructor is used, then methods that require
  /// a Jacobian or preconditioning will not be available.  You will be
  /// limited to simple algorithms like nonlinear-CG with no
  /// preconditioning.
  Group(
    Teuchos::ParameterList& printingParams,
    const Teuchos::RCP<InterfaceT>& i,
    const VectorT& x):
    utils_(printingParams),
    xVectorPtr_(Teuchos::rcp_dynamic_cast<VectorT>(x.clone(DeepCopy))),
    xVector_(*xVectorPtr_),
    RHSVectorPtr_(Teuchos::rcp_dynamic_cast<VectorT>(x.clone(ShapeCopy))),
    RHSVector(*RHSVectorPtr_),
    NewtonVectorPtr_(Teuchos::rcp_dynamic_cast<VectorT>(x.clone(ShapeCopy))),
    NewtonVector(*NewtonVectorPtr_),
    normNewtonSolveResidual(0),
    normF_(0.0),
    sharedInterfacePtr_(Teuchos::rcp(new NOX::SharedObject<InterfaceT, NOX::Pimpact::Group<InterfaceT> >(i))),
    sharedInterface_(*sharedInterfacePtr_),
    linearResidCompDisabled(false),
    noxCompF_( Teuchos::TimeMonitor::getNewCounter("NOX: compute F") ),
    noxUpdateX_( Teuchos::TimeMonitor::getNewCounter("NOX: update X") ),
    noxCompdx_( Teuchos::TimeMonitor::getNewCounter("NOX: solve dx") ),
    noxAssign_( Teuchos::TimeMonitor::getNewCounter("NOX: assign DF") ) {

    // Set all isValid flags to false
    resetIsValid();
  }


  /// \brief Copy constructor. If type is DeepCopy, takes ownership of valid
  /// shared linear system.
  Group( const NOX::Pimpact::Group<InterfaceT>& source, NOX::CopyType type = NOX::DeepCopy ):
    utils_(source.utils_),
    xVectorPtr_(Teuchos::rcp_dynamic_cast<VectorT>(source.xVector_.clone(type))),
    xVector_(*xVectorPtr_),
    RHSVectorPtr_(Teuchos::rcp_dynamic_cast<VectorT>(source.RHSVector.clone(type))),
    RHSVector(*RHSVectorPtr_),
    NewtonVectorPtr_(Teuchos::rcp_dynamic_cast<VectorT>(source.NewtonVector.clone(type))),
    NewtonVector(*NewtonVectorPtr_),
    sharedInterfacePtr_( source.sharedInterfacePtr_ ),
    sharedInterface_(*sharedInterfacePtr_),
    linearResidCompDisabled(source.linearResidCompDisabled),
    noxCompF_( Teuchos::TimeMonitor::getNewCounter("NOX: compute F") ),
    noxUpdateX_( Teuchos::TimeMonitor::getNewCounter("NOX: update X") ),
    noxCompdx_( Teuchos::TimeMonitor::getNewCounter("NOX: solve dx") ),
    noxAssign_( Teuchos::TimeMonitor::getNewCounter("NOX: assign DF") ) {

      switch (type) {
        case DeepCopy:
          isValidRHS_ = source.isValidRHS_;
          isValidJacobian_ = source.isValidJacobian_;
          isValidNewton = source.isValidNewton;
          isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;
          normNewtonSolveResidual = source.normNewtonSolveResidual;
          normF_ = source.normF_;
          isValidPreconditioner = source.isValidPreconditioner;
          isValidSolverJacOp = source.isValidSolverJacOp;

          // New copy takes ownership of the shared Jacobian for DeepCopy
          if( isValidJacobian_ )
            sharedInterface_.getObject(this);
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


  virtual NOX::Abstract::Group& operator=( const NOX::Abstract::Group& source ) {
    return operator=(dynamic_cast<const Group&> (source));
  }


  /// See operator=(const NOX::Abstract::Group&);
  virtual NOX::Abstract::Group& operator=( const NOX::Pimpact::Group<InterfaceT>& source ) {
    // Copy the xVector_
    xVector_ = source.xVector_;

    // Copy reference to sharedJacobian
    sharedInterfacePtr_ = source.sharedInterfacePtr_;

    // Update the isValidVectors
    isValidRHS_ = source.isValidRHS_;
    isValidNewton = source.isValidNewton;
    isValidJacobian_ = source.isValidJacobian_;
    isValidNormNewtonSolveResidual = source.isValidNormNewtonSolveResidual;
    isValidPreconditioner = source.isValidPreconditioner;
    isValidSolverJacOp = source.isValidSolverJacOp;

    // Only copy vectors that are valid
    if (isValidRHS_)
      RHSVector = source.RHSVector;

    if (isValidNewton)
      NewtonVector = source.NewtonVector;

    if (isValidNormNewtonSolveResidual)
      normNewtonSolveResidual = source.normNewtonSolveResidual;

    // If valid, this takes ownership of the shared Jacobian
    if (isValidJacobian_) sharedInterface_.getObject(this);

    normF_ = source.normF_;

    linearResidCompDisabled = source.linearResidCompDisabled;

    noxCompF_   = source.noxCompF_  ;
    noxUpdateX_ = source.noxUpdateX_;
    noxCompdx_  = source.noxCompdx_ ;
    noxAssign_  = source.noxAssign_ ;

    return *this;
  }


  /// \name "Compute" functions.
  //@{

  virtual void setX( const VectorT& y ) {

    resetIsValid();
    xVector_ = y;
  }

  virtual void setX( const NOX::Abstract::Vector& y ) {

    setX( dynamic_cast<const VectorT&> (y) );
  }


  virtual void computeX( const Group& grp, const VectorT& d, double step ) {

    resetIsValid();
    xVector_.update(1.0, grp.xVector_, step, d);
  }


  virtual void computeX( const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d,
      double step ) {

    Teuchos::TimeMonitor bla(*noxUpdateX_);

    // Cast to appropriate type, then call the "native" computeX
    const Group& pimpgrp = dynamic_cast<const Group&>( grp );
    const VectorT& pimpd =  dynamic_cast<const VectorT&>( d );

    computeX( pimpgrp, pimpd, step );
  }


  virtual NOX::Abstract::Group::ReturnType computeF() {

    Teuchos::TimeMonitor bla(*noxCompF_);

    if( isF() ) return Abstract::Group::Ok;

    if( !isJacobian() ) computeJacobian();

    NOX::Abstract::Group::ReturnType status;

    status = sharedInterfacePtr_->getObject(this)->computeF( xVector_.getField(), RHSVector.getField() );

    if (status != Abstract::Group::Ok ) {
      std::cout << "ERROR: Pimpact::Group::computeF() - fill failed!!!"
                << std::endl;
      throw "NOX Error: Fill Failed";
    }

    normF_ = RHSVectorPtr_->norm();

    isValidRHS_ = true;

    return status;
  }


  virtual NOX::Abstract::Group::ReturnType computeJacobian() {

    Teuchos::TimeMonitor bla(*noxAssign_);

    // Skip if the Jacobian is already valid
    if( isJacobian() ) return Abstract::Group::Ok;

    // Fill the Jacobian
    NOX::Abstract::Group::ReturnType status;
    status = sharedInterfacePtr_->getObject(this)->computeJacobian( xVector_.getField() );

    if (status != NOX::Abstract::Group::Ok) {
      std::cout << "ERROR: NOX::Pimpact::Group::computeJacobian() - fill failed!!!"
                << std::endl;
      throw "NOX Error: Fill Failed";
    }
    // Update status of Jacobian wrt solution vector
    isValidJacobian_ = true;

    return status;
  }


  virtual NOX::Abstract::Group::ReturnType computeGradient() {

    return Abstract::Group::NotDefined;
  }


  virtual NOX::Abstract::Group::ReturnType computeNewton( Teuchos::ParameterList& params ) {

    Teuchos::TimeMonitor bla( *noxCompdx_ );

    if (isNewton())
      return Abstract::Group::Ok;

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
    status = applyJacobianInverse( params, RHSVector, NewtonVector );

    // Scale soln by -1
    NewtonVector.scale(-1.0);
    //
    // Update state EVEN IF LINEAR SOLVE FAILED
    // We still may want to use the vector even it it just missed it's
    isValidNewton = true;

    // Compute the 2-norm of the linear solve residual ||Js+f||
    // Can be disabled, but then disallows inexact Newton methods
    if (!linearResidCompDisabled)
      computeNormNewtonSolveResidual();

    // Return solution
    return status;
  }


  //@}
  /// \name Jacobian operations.
  ///
  /// Operations using the Jacobian matrix. These may not be defined in
  /// matrix-free scenarios.
  //@{

  virtual NOX::Abstract::Group::ReturnType
  applyJacobian( const VectorT& input, VectorT& result ) const {
    //    // Check validity of the Jacobian
    if (!isJacobian())
      return Abstract::Group::BadDependency;

    // Apply the Jacobian
    //    NOX::Abstract::Group::ReturnType status = sharedInterface_.getObject()->applyJacobian( input.getConstField(), result.getField() );
    NOX::Abstract::Group::ReturnType status = sharedInterface_.getObject(this)->applyJacobian( input.getConstField(), result.getField() );

    return status;
  }

  virtual NOX::Abstract::Group::ReturnType
  applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const {
    const VectorT& pimpinput =  dynamic_cast<const VectorT&> (input);
    VectorT& pimpresult = dynamic_cast<VectorT&> (result);

    return applyJacobian( pimpinput, pimpresult);
  }


  //  virtual NOX::Abstract::Group::ReturnType
  //  applyJacobianTranspose(const Vector& input, Vector& result) const {
  //    // Check validity of the Jacobian
  //    if( !isJacobian() )
  //      return Abstract::Group::BadDependency;
  //
  //    bool status = sharedInterface_.getObject()->applyJacobian( input, result, Belos::Trans );
  //
  //    return status == true ? Abstract::Group::Ok : Abstract::Group::Failed;
  //  }
  virtual NOX::Abstract::Group::ReturnType
  applyJacobianTranspose( const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result ) const {
    return Abstract::Group::NotDefined;
  }



  /// \brief Applies the inverse of the Jacobian matrix to the given
  /// input vector and puts the answer in result.
  ///
  /// Computes
  /// \f[ v = J^{-1} u, \f]
  /// where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is
  /// the result vector.
  /// The "Tolerance" parameter specifies that the solution should be such that
  /// \f[ \frac{\| J v - u \|_2}{\max \{ 1, \|u\|_2\} } < \mbox{Tolerance} \f]
  ///
  /// \return
  /// <ul>
  /// <li> NOX::Abstract::Group::NotDefined - Returned by default implementation
  /// in NOX::Abstract::Group
  /// <li> NOX::Abstract::Group::BadDependency - If \f$J\f$ has not been computed
  /// <li> NOX::Abstract::Group::NotConverged - If the linear solve fails to satisfy the "Tolerance"
  /// specified in \c params
  /// <li> NOX::Abstract::Group::Failed - If the computation fails
  /// <li> NOX::Abstract::Group::Ok - Otherwise
  /// </ul>
  /// The parameter "Tolerance" may be added/modified in the list of
  /// parameters - this is the ideal solution tolerance for an iterative
  /// linear solve.
  /// The parameter "Reuse Preconditioner" is a boolean that tells the group to turn off control of preconditioner recalculation.  This is a dangerous flag but can really speed the computations if the user knows what they are doing.  Toggling this flag is left to the user (ideally it should be done through a status test).  Defaults to false.
  virtual NOX::Abstract::Group::ReturnType
  applyJacobianInverse( Teuchos::ParameterList &params, const VectorT& input, VectorT& result ) const {
    if (!isJacobian())
      return Abstract::Group::BadDependency;

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
    //            createPreconditioner(xVector_, p, false);
    //        isValidPreconditioner = true;
    //      }
    //      else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_RECOMPUTE) {
    //        sharedLinearSystem.getObject(this)->recomputePreconditioner(xVector_, p);
    //        isValidPreconditioner = true;
    //      }
    //      else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_REUSE) {
    //        // Do Nothing!!!
    //      }
    //    }

    NOX::Abstract::Group::ReturnType status = sharedInterface_.getObject(this)->applyJacobianInverse(params, input.getConstField(), result.getField() );

    return status;
  }

  virtual NOX::Abstract::Group::ReturnType
  applyJacobianInverse( Teuchos::ParameterList& params, const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result ) const {
    const VectorT& pimpInput = dynamic_cast<const VectorT&>(input);
    VectorT& pimpResult = dynamic_cast<VectorT&>(result);
    return applyJacobianInverse(params, pimpInput, pimpResult);
  }


  virtual NOX::Abstract::Group::ReturnType
  applyRightPreconditioning(
    bool useTranspose,
    Teuchos::ParameterList& params,
    const VectorT& input,
    VectorT& result ) const {

    //    params.print();
    //    std::cout << "Call applyRightPrecon.... !!!!!!\n";

    return sharedInterfacePtr_->getObject(this)->applyPreconditioner( xVector_.getField(), RHSVector.getField() );
  }

  virtual NOX::Abstract::Group::ReturnType
  applyRightPreconditioning(bool useTranspose,
                            Teuchos::ParameterList& params,
                            const NOX::Abstract::Vector& input,
                            NOX::Abstract::Vector& result) const {
    //    return Abstract::Group::NotDefined;
    const VectorT& pimpInput = dynamic_cast<const VectorT&>(input);
    VectorT& pimpResult = dynamic_cast<VectorT&>(result);

    return applyRightPreconditioning(useTranspose, params, pimpInput, pimpResult);
  }

  //@}
  /// \name "Is" functions
  ///
  /// Checks to see if various objects have been computed. Returns true
  /// if the corresponding "compute" function has been called since the
  /// last update to the solution vector (via instantiation or
  /// computeX).
  //@{

  virtual bool isF() const {
    return isValidRHS_;
  }
  virtual bool isJacobian() const {
    return ((sharedInterface_.isOwner(this)) && (isValidJacobian_));
  }
  virtual bool isGradient() const {
    return false;
  }
  virtual bool isNewton() const {
    return isValidNewton;
  }


  /// \brief Returns true if the value of the Norm of the linear model.
  ///
  /// for a full Newton step ||Js + f|| is valid with respect to the
  /// current solution vector.
  virtual bool isNormNewtonSolveResidual() const {
    return isValidNormNewtonSolveResidual;
  }

  /// \brief Returns true if an explicitly constructed preconditioner.
  ///
  /// exists (i.e. one that is computed and saved for further use in
  /// multiple calls to applyRightPreconditioner).
  virtual bool isPreconditioner() const {
    return false;
    //    return ((sharedInterface_.isOwner(this)) && (isValidPreconditioner) &&
    //        (sharedInterface_.getObject(this)->isPreconditionerConstructed()));
  }


  //@}
  /// \name "Get" functions
  ///
  /// Note that these function do not check whether or not the vectors
  /// are valid. Must use the "Is" functions for that purpose.
  //@{

  virtual const NOX::Abstract::Vector& getX() const {
    return xVector_;
  }


  virtual const NOX::Abstract::Vector& getF() const {
    if (!isF()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::getF() - invalid RHS" << std::endl;
      throw "NOX Error";
    }

    return RHSVector;
  }


  /// return pre computed 2-Norm of RHS F
  virtual double getNormF() const {
    if (!isF()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::getNormF() - invalid RHS" << std::endl;
      throw "NOX Error";
    }

    return normF_;
  }


  virtual const NOX::Abstract::Vector& getGradient() const {
    std::cerr << "ERROR: NOX::Pimpact::Group::getGradient() - gradient not implemented" << std::endl;
    throw "NOX Error";
  }


  virtual const NOX::Abstract::Vector& getNewton() const {
    if (!isNewton()) {
      std::cerr << "ERROR: NOX::Pimpact::Group::getNewton() - invalid Newton vector" << std::endl;
      throw "NOX Error";
    }

    return *NewtonVectorPtr_;
  }


  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const {
    return xVectorPtr_;
  };
  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const {
    return RHSVectorPtr_;
  };
  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const {
    std::cerr << "ERROR: NOX::Pimpact::Group::getGradient() - gradient not implemented" << std::endl;
    throw "NOX Error";
  };

  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const {
    return NewtonVectorPtr_;
  };


  /// \brief Returns the 2-norm of the residual of the linear model used in the.
  ///
  /// Newton solve computation, ||Js+f||.  This does not account for
  /// line search adjustments to the step length!
  virtual NOX::Abstract::Group::ReturnType getNormLastLinearSolveResidual(double & residual) const {
    // Make sure value is not already calculated
    if (isValidNormNewtonSolveResidual) {
      residual = normNewtonSolveResidual;
      return NOX::Abstract::Group::Ok;
    }

    // Otherwise give warning since a Newton direction has not been calculated
    // wrt this solution group
    if( utils_.isPrintType(Utils::Warning) ) {
      std::cout << "ERROR: NOX::Epetra::Group::getNormLastLinearSolveResidual() - "
                << "Group has not performed a Newton solve corresponding to this "
                << "solution vector, or disableLinearSolveResidual(true) was set!" << std::endl;
    }
    return NOX::Abstract::Group::BadDependency;
  }


  //@}

  virtual Teuchos::RCP<NOX::Abstract::Group>
  clone(CopyType type = DeepCopy) const {
    Teuchos::RCP<NOX::Abstract::Group> newgrp =
      Teuchos::rcp( new NOX::Pimpact::Group<InterfaceT>( *this, type ) );

    return newgrp;
  }

  /// Return the Interface.
  virtual Teuchos::RCP<InterfaceT> getInterface() const {
    return sharedInterface_.getObject(this);
  }


  /// \brief Computes the condition number of the Jacobian matrix.
  ///
  /// Uses GMRES to estimate the condtion number.  The tolerance and
  /// maxIters are used to control the GMRES iterations.  Typically the
  /// solves do not have to be tight to get the estimate.
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
  //    }
  //    return NOX::Abstract::Group::Ok;
  //  }


  //  /// Returns the condition number of the Jacobian matrix.
  //  virtual double getJacobianConditionNumber() const {
  //    if( !isConditionNumber() ) {
  //      std::cerr << "ERROR: NOX::Epetra::Group::getJacobianConditionNumber()"
  //          << " - condition number has not yet been computed!" << std::endl;
  //      throw "NOX Error";
  //    }
  //    return normF_;
  //  }


  /// \brief Sets option to disable linear resid computation.
  ///
  /// If disabled, this saves on a MatVec per Newton but disallows inexact
  /// Newton methods.
  virtual void disableLinearResidualComputation(const bool disableChoice) {
    linearResidCompDisabled = disableChoice;
  }

protected:
  /// resets the isValid flags to false
  virtual void resetIsValid() {
    isValidRHS_ = false;
    isValidJacobian_ = false;
    isValidNewton = false;
    isValidNormNewtonSolveResidual = false;
    isValidPreconditioner = false;
    isValidSolverJacOp = false;
    return;
  }


  /// \brief Computes the 2-norm of the residual of the linear model used in
  /// the Newton solve computation, ||Js+f||.
  virtual bool computeNormNewtonSolveResidual() {
    //    // Make sure value is not already calculated
    //    if (isValidNormNewtonSolveResidual)
    //      return true;
    //
    //    // Make sure NewtonVector and RHSVector are valid
    //    // We could return false, but for now we will throw errors
    //    if (!isValidRHS_) {
    //      std::cerr << "ERROR: NOX::Pimpact::Group::computeNormNewtonSolveResidual() - invalid RHS"
    //          << std::endl;
    //      throw "NOX Error";
    //    }
    //    if (!isValidNewton) {
    //      std::cerr << "ERROR: NOX::Pimpact::Group::computeNormNewtonSolveResidual() - invalid "
    //          << "Newton direction" << std::endl;
    //      throw "NOX Error";
    //    }
    //
    //    sharedInterface_.getObject()->applyJacobian(NewtonVector, tmpNoxVector);
    //    tmpNoxVector.update(1.0, RHSVector, 1.0);
    //    normNewtonSolveResidual = tmpNoxVector.norm();
    //
    //    isValidNormNewtonSolveResidual = true;

    return true;
  }

}; // end of class Group



/// \relates Group
template< class InterfaceT >
Teuchos::RCP<NOX::Pimpact::Group<InterfaceT> > createGroup(
  const Teuchos::RCP<Teuchos::ParameterList>& list,
  const Teuchos::RCP<InterfaceT>& i,
  const Teuchos::RCP<typename NOX::Pimpact::Group<InterfaceT>::VectorT>& x ) {

  return Teuchos::rcp( new NOX::Pimpact::Group<InterfaceT>( *list, i, *x ) );
}


} // end of namespace Pimpact
} // end of namespace NOX


#endif // end of #ifndef NOX_PIMPACT_GROUP_HPP
