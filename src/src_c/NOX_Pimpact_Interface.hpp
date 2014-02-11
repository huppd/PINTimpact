#pragma once
#ifndef NOX_PIMPACT_INTERFACE_HPP
#define NOX_PIMPACT_INTERFACE_HPP

#include "NOX_Common.H"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_CompoundField.hpp"

#include "Pimpact_MultiField.hpp"

#include "Pimpact_LinearProblem.hpp"

// Forward declarations
//class Epetra_Vector;

namespace NOX {
namespace PIMPACT {

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

  typedef double S;
  typedef int O;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::MultiField<VF> BVF;
  typedef Pimpact::MultiField<SF> BSF;

  typedef Pimpact::OperatorMV< Pimpact::DtL<Scalar,Ordinal> >  DTL;
  typedef Pimpact::OperatorMV< Pimpact::Div_DtLinv_Grad<Scalar,Ordinal> >  Schur;
  typedef Pimpact::OperatorMV< Pimpact::Div<Scalar,Ordinal> >  D;
  typedef Pimpact::OperatorMV< Pimpact::Grad<Scalar,Ordinal> > G;


public:
  typedef Pimpact::CompoundField< BVF, BSF> Field;

protected:

  Teuchos::RCP<BVF> f_;
  Teuchos::RCP<Pimpact::LinearProblem<S,BVF,DTL> > dTL_;
  Teuchos::RCP<D> div_;
  Teuchos::RCP<G> grad_;
  Teuchos::RCP<Pimpact::LinearProblem<S,BSF,Schur> > schur_;

public:

  /// Constructor
  Interface():
    f_(Teuchos::null), dTL_(Teuchos::null),div_(Teuchos::null),
    grad_(Teuchos::null),schur_(Teuchos::null) {};

  Interface(  Teuchos::RCP<BVF> f,
      Teuchos::RCP<Pimpact::LinearProblem<S,BVF,DTL> > dTL,
      Teuchos::RCP<Pimpact::LinearProblem<S,BVF,Schur> > schur ):
        f_( f->clone() ),
        div_( Pimpact::createOperatorMV<Pimpact::Div<S,O> >() ),
        dTL_( dTL),
        grad_( Pimpact::createOperatorMV<Pimpact::Grad<S,O> >() ),
        shur_( shur ) {};

  /// Destructor
  ~Required() {};

  /// Compute the function, F, given the specified input vector x. Returns true if computation was successful.
  bool computeF(const Field& x, Field& F ) {
    auto xv = x.getVField();
    auto xs = x.getSField();
    auto yv = y.getVField();
    auto ys = y.getSField();

    auto tempv = xv.clone();

    dTL_->apply( *xv, *temvp );
    grad_->apply( *xs, *yv );

    tempv->add( -1., *tempv, 1., *f );
    yv->add(yv,tempv);

    div_->apply( *xv, *ys );


    return( true );
  }

  /// Compute the Jacobian Operator, given the specified input vector x. Returns true if computation was successful.
  bool computeJacobian( const Field& x ) { return( true ); }
  bool applyJacobian( const Field& x, Field& y, Belos::ETrans type=Belos::NOTRANS ) {

    return( true );
  }


  bool applyJacobianInverse( /* param, */ const Field& x, Field& y ) {

    auto xv = x.getVField();
    auto xs = x.getSField();
    auto yv = y.getVField();
    auto ys = y.getSField();

    auto tempv = xv->clone();
    auto temps = xs->clone();
    // solve stationary stokes
    dTL_->solve( *tempv, *f_ );
    //  tempv->write( 2000 );

    div_->apply( *xv, *temps );
    //temps->write( 2002 );

//    solverParams = Pimpact::createLinSolverParameter( solver_name_1, 1.e-6*l1*l2/n1/n2 );
//    solverParams->get()->set( "Output Stream", outLap2 );
//    solverParams->get()->set ("Verbosity", int( Belos::Errors) );
//    lap_problem->setParameters( solverParams->get() );
    schur_->solve( *ys, *temps );

    grad_->apply( *ys, *tempv );
    //tempv->write(2006);

    tempv->add( -1., *tempv, 1., *f );

//    solverParams->get()->set ("Verbosity",  Belos::Errors + Belos::Warnings + Belos::IterationDetails +
//        Belos::OrthoDetails + Belos::FinalSummary + Belos::TimingDetails +
//        Belos::StatusTestDetails + Belos::Debug );
//    solverParams->get()->set( "Output Stream", outLap2 );
//    lap_problem->setParameters( solverParams->get() );

    dTL_->solve( *ys, *tempv );
    return( true );
  }

}; // end of class Interface

} // end of namespace PIMPACT
} // end of namespace NOX

#endif // end of #ifndef NOX_EPETRA_INTERFACE_REQUIRED_HPP

