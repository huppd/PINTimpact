#pragma once
#ifndef NOX_PIMPACT_STATUSTEST_HPP
#define NOX_PIMPACT_STATUSTEST_HPP


#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_Vector.H"        // for NormType
#include "NOX_Common.H"
#include "NOX_Solver_Generic.H"
#include "NOX_StatusTest_Generic.H"    // base class
#include "NOX_Utils.H"                  // for std::cerr

#include "NOX_Pimpact_Group.hpp"



namespace NOX {
namespace Pimpact {


template<class InterfaceT>
class RefinementTest : public NOX::StatusTest::Generic {

private:

  /// %Status
  NOX::StatusTest::StatusType status_;

  //! Tolerance required for convergence.
  double tolerance_;

  /// \brief Norm of F to be compared to trueTolerance
  double normF_;

  /// \brief Norm of refined F
  double normRF_;

  //! Ostream used to print errors
  NOX::Utils utils_;

  //! Ostream to print out norms
  Teuchos::RCP<std::ostream> out_;

public:

  /// Constructor
  RefinementTest( double tolerance=1., const Teuchos::RCP<std::ostream>& out=Teuchos::null,
      const NOX::Utils* u = NULL ):
    tolerance_(tolerance),
    normF_(0.0),
    normRF_(0.0),
    out_(out) {

      if (u != NULL)
        utils_ = *u;
    }


  /// Destructor.
  virtual ~RefinementTest() {};

  /// derived
  virtual NOX::StatusTest::StatusType
  checkStatus(const NOX::Solver::Generic& problem, NOX::StatusTest::CheckType checkType) {

    int nf = 0;

    if( checkType == NOX::StatusTest::None ) {
      normF_ = 0.0;
      normRF_ = 0.0;
      status_ = NOX::StatusTest::Unevaluated;
    }
    else {
      normF_ = computeFNorm( problem.getSolutionGroup() );
      normRF_ = computeRFNorm( problem.getSolutionGroup() );

    Teuchos::RCP<const NOX::Abstract::Vector> x =
      problem.getSolutionGroup().getXPtr();

    nf = Teuchos::rcp_dynamic_cast<const NOX::Pimpact::Vector<typename
      InterfaceT::FieldT> >(x)->getConstFieldPtr()->space()->nGlo(3) + 1;


      if( normRF_==0. )
        status_ = NOX::StatusTest::Unconverged;
      else
        status_ = ( normF_ < nf*tolerance_*normRF_ ) ?
          NOX::StatusTest::Converged :
          NOX::StatusTest::Unconverged;
    }

    if( !out_.is_null() )
      (*out_) << problem.getNumIterations() << nf << normF_ << "\t" << normRF_ << "\n";

    return status_;
  }

  /// derived
  virtual NOX::StatusTest::StatusType getStatus() const {
    return status_;
  }

  virtual std::ostream& print(std::ostream& stream, int indent = 0) const {

    for (int j = 0; j < indent; j ++)
      stream << ' ';
    stream << status_;
    stream << "F-Norm/FR-Norm = " << Utils::sciformat(normF_,3);
    stream << " < " << Utils::sciformat(tolerance_*normRF_, 3);
    stream << "\n";

    for (int j = 0; j < indent; j ++)
      stream << ' ';
    stream << std::setw(13) << " ";

    stream << std::endl;

    return stream;
  }

  /// @name Reset Functions
  ///   Used to change the tolerances in the status test after construction.
  ///@{

  /// Resets the user specified absolute or relative tolerance.
  virtual void reset(double tolerance) {

    tolerance_ = tolerance;
  }

  //@}

  /// @name Accessor Functions
  /// Used to query current values of variables in the status test.
  ///@{

  ///// Returns the value of the F-norm computed in the last call to checkStatus.
  //virtual double getRefinementTest() const {
    //return normF_;
  //}

  ///@}

private:

  /// \brief Calculate the norm of F for the given group according to the scaling type,
  /// norm type, and tolerance type.
  ///
  /// \note Returns -1.0 if F(x) has not been calculated for the given grp (i.e.,
  /// grp.isF() is false).
  double computeFNorm(const NOX::Abstract::Group& grp) {

    if (!grp.isF())
      return -1.0;

    return grp.getNormF();
  }

  /// \brief Calculate the norm of the refinement F for the given group according to the
  /// scaling type, / norm type, and tolerance type.
  ///
  /// \note Returns -1.0 if F(x) has not been calculated for the given grp (i.e.,
  /// grp.isF() is false).
  double computeRFNorm(const NOX::Abstract::Group& grp) {

    if (!grp.isJacobian()) const_cast<NOX::Abstract::Group&>(grp).computeJacobian();

    Teuchos::RCP<InterfaceT> sharedInterface =
      dynamic_cast<const NOX::Pimpact::Group<InterfaceT>& >(grp).getInterface();

    Teuchos::RCP<const typename InterfaceT::OperatorT> op = sharedInterface->getOperatorPtr();

    Teuchos::RCP<const NOX::Abstract::Vector> x = grp.getXPtr();

    if( Teuchos::rcp_dynamic_cast<const NOX::Pimpact::Vector<typename InterfaceT::FieldT>
      >(x)->getConstFieldPtr()->space()->nGlo(3)>0 )
      return op->getOperatorPtr()->getOpV2V()->compRefRes( Teuchos::rcp_dynamic_cast<const
          NOX::Pimpact::Vector<typename InterfaceT::FieldT>
          >(x)->getConstFieldPtr()->getField(0).getVField() );
    else
      return 0.;

  }
};


} // namespace Pimpact
} // namespace NOX


#endif // end of #ifndef NOX_PIMPACT_STATUSTEST_HPP
