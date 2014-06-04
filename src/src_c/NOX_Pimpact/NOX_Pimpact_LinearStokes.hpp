#pragma once
#define NOX_PIMPACT_LINEARSTOKES_HPP
#ifndef NOX_PIMPACT_LINEARSTOKES_HPP

#include "NOX_Common.H"
#include "NOX_Abstract_Group.H"   // base class

#include "Pimpact_VectorField.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_CompoundField.hpp"

#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_LinearProblem.hpp"

// Forward declarations
//class Epetra_Vector;



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
class LinearStokes {

public:

  typedef double S;
  typedef int O;
  typedef ::Pimpact::VectorField<S,O> VF;
  typedef ::Pimpact::ScalarField<S,O> SF;
  typedef ::Pimpact::ModeField<VF> MVF;
  typedef ::Pimpact::ModeField<SF> MSF;
  typedef ::Pimpact::MultiField<MVF> BVF;
  typedef ::Pimpact::MultiField<MSF> BSF;
  typedef ::Pimpact::CompoundField< BVF, BSF> Field;

  typedef ::Pimpact::MultiOpWrap< ::Pimpact::DtL<S,O> >  DTL;
  typedef ::Pimpact::MultiOpWrap< ::Pimpact::DivDtLinvGrad<S,O> >  Schur;
  typedef ::Pimpact::Div<S,O>   DD;
  typedef ::Pimpact::Grad<S,O>  GG;
  typedef ::Pimpact::MultiOpWrap< ::Pimpact::ModeOpWrap<DD> > D;
  typedef ::Pimpact::MultiOpWrap< ::Pimpact::ModeOpWrap<GG> > G;

  typedef ::Pimpact::LinearProblem<BVF> LP_DTL;
  typedef ::Pimpact::LinearProblem<BSF> LP_Schur;

  typedef NOX::Pimpact::Vector<Field> Vector;

protected:

  Teuchos::RCP<BVF> fu_;
  Teuchos::RCP<BSF> fp_;
  Teuchos::RCP<LP_DTL> lp_DTL_;
  Teuchos::RCP<D> div_;
  Teuchos::RCP<G> grad_;
  Teuchos::RCP<LP_Schur> lp_Schur_;

public:

  /// Constructor
  LinearStokes():
    fu_(Teuchos::null),fp_(Teuchos::null),
    lp_DTL_(Teuchos::null),div_(Teuchos::null),
    grad_(Teuchos::null),lp_Schur_(Teuchos::null) {};

  LinearStokes(
      Teuchos::RCP<BVF> fu,
      //      Teuchos::RCP<BVF> tempu,
      Teuchos::RCP<BSF> fp,
      //      Teuchos::RCP<BSF> tempp,
      Teuchos::RCP<LP_DTL > lp_DTL,
      Teuchos::RCP<LP_Schur > lp_Schur ):
        fu_( fu ),
        fp_( fp ),
        lp_DTL_( lp_DTL),
        div_( ::Pimpact::createMultiModeOpWrap<DD>() ),
        grad_( ::Pimpact::createMultiModeOpWrap<GG>() ),
        lp_Schur_( lp_Schur ) {};

  /// Destructor
  ~LinearStokes() {};

  /// Compute the function, F, given the specified input vector x. Returns true if computation was successful.
  NOX::Abstract::Group::ReturnType computeF(const Field& x, Field& f ) {
    auto xv = x.getConstVField();
    auto xs = x.getConstSField();
    auto yv = f.getVField();
    auto ys = f.getSField();

    //    yv->random();
    //    ys->random();
    auto tempv = xv->clone();
    //
    //    lp_DTL_->apply( xv, tempv );
    //    lp_DTL_->apply( *xv, *tempv );
    lp_DTL_->getProblem()->getOperator()->apply( *xv, *tempv );

    grad_->apply( *xs, *yv );

    yv->add(  1., *tempv, 1., *yv );
    yv->add( -1., *fu_,   1., *yv );

    div_->apply( *xv, *ys );
    ys->add( -1., *fp_,   1., *ys );


    return( NOX::Abstract::Group::Ok );
  }

  /// Compute the Jacobian Operator, given the specified input vector x. Returns true if computation was successful.
  NOX::Abstract::Group::ReturnType computeJacobian( const Field& x ) {
    return( NOX::Abstract::Group::Ok );
  }

  NOX::Abstract::Group::ReturnType applyJacobian( const Field& x, Field& y, Belos::ETrans type=Belos::NOTRANS ) {
    return( computeF( x, y ) );
  }


  NOX::Abstract::Group::ReturnType applyJacobianInverse( Teuchos::ParameterList &params, const Field& x, Field& y ) {
    params.print();

    auto xv = x.getConstVField();
    auto xs = x.getConstSField();
    auto yv = y.getVField();
    auto ys = y.getSField();

    auto tempv = xv->clone();
    auto temps = xs->clone();

    // solve stationary stokes
    lp_DTL_->solve( tempv, fu_ );
    //  tempv->write( 2000 );

    div_->apply( *tempv, *temps );
    //temps->write( 2002 );
    //
    //    solverParams = Pimpact::createLinSolverParameter( solver_name_1, 1.e-6*l1*l2/n1/n2 );
    //    solverParams->get()->set( "Output Stream", outLap2 );
    //    solverParams->get()->set ("Verbosity", int( Belos::Errors) );
    //    lap_problem->setParameters( solverParams->get() );
    lp_Schur_->solve( ys, temps );

    grad_->apply( *ys, *tempv );
    //tempv->write(2006);
    //
    tempv->add( -1., *tempv, 1., *fu_ );
    //
    ////    solverParams->get()->set ("Verbosity",  Belos::Errors + Belos::Warnings + Belos::IterationDetails +
    ////        Belos::OrthoDetails + Belos::FinalSummary + Belos::TimingDetails +
    ////        Belos::StatusTestDetails + Belos::Debug );
    ////    solverParams->get()->set( "Output Stream", outLap2 );
    ////    lap_problem->setParameters( solverParams->get() );
    //
    lp_DTL_->solve( yv, tempv );
    return( NOX::Abstract::Group::Ok );
  }

  NOX::Abstract::Group::ReturnType applyPreconditioner( const Field& x, Field& y ) {
    return( NOX::Abstract::Group::NotDefined );
  }

}; // end of class LinearStokes

Teuchos::RCP<LinearStokes> createLinearStokes(
    Teuchos::RCP<LinearStokes::BVF> fu,
    Teuchos::RCP<LinearStokes::BSF> fp,
    Teuchos::RCP<LinearStokes::LP_DTL> lp_DTL,
    Teuchos::RCP<LinearStokes::LP_Schur> lp_Schur ) {
  return(
      Teuchos::rcp( new LinearStokes(fu,fp,lp_DTL,lp_Schur) )
  );
}

} // end of namespace Pimpact
} // end of namespace NOX

#endif // end of #ifndef NOX_PIMPACT_LINEARSTOKES_HPP

