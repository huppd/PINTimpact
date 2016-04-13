#pragma once
#ifndef NOX_PIMPACT_INTERFACE_HPP
#define NOX_PIMPACT_INTERFACE_HPP


#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_OperatorBase.hpp"




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

  using Field = F;
  //using S = typename Field::SpaceT::Scalar;
  //using O = typename Field::SpaceT::Ordinal;
  using Op = ::Pimpact::OperatorBase<Field>;

protected:

  Teuchos::RCP<Field> fu_;
  Teuchos::RCP<Op>    op_;
  Teuchos::RCP<Op>   jop_;

public:

  /// Constructor
  Interface(
      Teuchos::RCP<Field> fu=Teuchos::null,
      Teuchos::RCP<Op>    op=Teuchos::null,
      Teuchos::RCP<Op>   jop=Teuchos::null ):
        fu_( fu ),
        op_(op),
        jop_(jop) {};

  /// Compute the function, F, given the specified input vector x.
  NOX::Abstract::Group::ReturnType computeF(const Field& x, Field& f ) {

    op_->assignField( x );
    op_->apply( x, f );
    f.add( 1., f, -1., *fu_ );
//		f.level();
    return( NOX::Abstract::Group::Ok );

  }


  /// \brief Compute the Jacobian Operator, given the specified input vector x.
  NOX::Abstract::Group::ReturnType computeJacobian( const Field& x ) {

	 jop_->assignField( x );
    return( NOX::Abstract::Group::Ok );

  }


  NOX::Abstract::Group::ReturnType applyJacobian( const Field& x, Field& y, Belos::ETrans type=Belos::NOTRANS ) {
    return( NOX::Abstract::Group::NotDefined );
  }


  NOX::Abstract::Group::ReturnType applyJacobianInverse( Teuchos::ParameterList &params, const Field& x, Field& y ) {

    jop_->apply( x, y );
    return( NOX::Abstract::Group::Ok );

  }


  NOX::Abstract::Group::ReturnType applyPreconditioner( const Field& x, Field& y ) {
    return( NOX::Abstract::Group::NotDefined );
  }

}; // end of class Interface



/// \relates Interface
template<class Field>
Teuchos::RCP< Interface<Field> > createInterface(
    Teuchos::RCP<Field> fu=Teuchos::null,
    Teuchos::RCP< ::Pimpact::OperatorBase<Field> >  op=Teuchos::null,
    Teuchos::RCP< ::Pimpact::OperatorBase<Field> > jop=Teuchos::null ) {
  return(
      Teuchos::rcp( new Interface<Field>(fu,op,jop) )
  );
}


} // end of namespace Pimpact
} // end of namespace NOX


#ifdef COMPILE_ETI
#include "Pimpact_Space.hpp"
#include "Pimpact_Fields.hpp"
extern template class NOX::Pimpact::Interface<
Pimpact::MultiField<Pimpact::VectorField<Pimpact::Space<double,int,3,2>
> > >; 
extern template class NOX::Pimpact::Interface<
		Pimpact::CompoundField<
			Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<Pimpact::Space<double,int,3,4> > > >,
			Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<Pimpact::Space<double,int,3,4> > > >
		>
	>; 
extern template class NOX::Pimpact::Interface<
		Pimpact::MultiField<
			Pimpact::CompoundField<
				Pimpact::MultiHarmonicField<Pimpact::VectorField<Pimpact::Space<double,int,3,4> > >,
				Pimpact::MultiHarmonicField<Pimpact::ScalarField<Pimpact::Space<double,int,3,4> > > 
			>
		>
	>; 
extern template class NOX::Pimpact::Interface<
		Pimpact::MultiField<
			Pimpact::CompoundField<
				Pimpact::TimeField<Pimpact::VectorField<Pimpact::Space<double,int,4,4> > >,
				Pimpact::TimeField<Pimpact::ScalarField<Pimpact::Space<double,int,4,4> > > 
			>
		>
	>; 
#endif

#endif // end of #ifndef NOX_PIMPACT_INTERFACE_HPP
