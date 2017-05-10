#pragma once
#ifndef NOX_PIMPACT_INTERFACE_HPP
#define NOX_PIMPACT_INTERFACE_HPP


#include <mpi.h>
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
/// \tparam FieldT hast to be of type \c Pimpact::MultiField.
template< class FT, class OpT=::Pimpact::OperatorBase<FT>, class IOpT=::Pimpact::OperatorBase<FT> >
class Interface {

public:

	using FieldT = FT;

protected:

	Teuchos::RCP<FieldT> fu_;
	Teuchos::RCP<FieldT> sol_;
	Teuchos::RCP<OpT>    op_;
	Teuchos::RCP<IOpT>   jopInv_;

	Teuchos::RCP<std::ostream> eStream;

public:

	/// Constructor
	Interface(
			Teuchos::RCP<FieldT> fu=Teuchos::null,
			Teuchos::RCP<OpT>    op=Teuchos::null,
			Teuchos::RCP<IOpT>   jop=Teuchos::null,
			Teuchos::RCP<FieldT> sol=Teuchos::null	):
		fu_( fu ),
		sol_( sol ),
		op_(op),
		jopInv_(jop) {

			int world_rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
			if( sol_!=Teuchos::null && world_rank==0 )
				eStream = Teuchos::rcp( new std::ofstream( "errorIter.txt" ) );
			else
				eStream = Teuchos::null;
		};

	/// Compute the function, F, given the specified input vector x.
	NOX::Abstract::Group::ReturnType computeF( const FieldT& x, FieldT& f ) {

		if( !sol_.is_null() ) {
			auto temp = x.getField(0).getVField().clone();
			double solNorm = sol_->getField(0).getVField().norm();
			temp->add( 1., sol_->getField(0).getVField(), -1., x.getField(0).getVField() );
			//temp->add( 1., *sol_, -1., x );
			double error = temp->norm()/solNorm;
			if( eStream!=Teuchos::null )
				*eStream << error << "\n";
		}

		op_->assignField( x );
		op_->apply( x, f );
		f.add( 1., f, -1., *fu_ );

		return( NOX::Abstract::Group::Ok );
	}


	/// \brief Compute the Jacobian Operator, given the specified input vector x.
	NOX::Abstract::Group::ReturnType computeJacobian( const FieldT& x ) {

		jopInv_->assignField( x );
		return( NOX::Abstract::Group::Ok );
	}


	NOX::Abstract::Group::ReturnType applyJacobian( const FieldT& x, FieldT& y, const Belos::ETrans& type=Belos::NOTRANS ) {
		return( NOX::Abstract::Group::NotDefined );
	}


	NOX::Abstract::Group::ReturnType applyJacobianInverse( Teuchos::ParameterList &params, const FieldT& x, FieldT& y ) {

		double atol = params.get<double>("Tolerance");
		if( atol>0. ) {
			double res = x.norm();

			Teuchos::RCP<Teuchos::ParameterList> para = Teuchos::parameterList("Linear Solver");
			para->set<double>( "Convergence Tolerance", atol*res );

			jopInv_->setParameter( para );
		}

		jopInv_->apply( x, y );
		return( NOX::Abstract::Group::Ok );
	}


	NOX::Abstract::Group::ReturnType applyPreconditioner( const FieldT& x, FieldT& y ) {
		return( NOX::Abstract::Group::NotDefined );
	}

}; // end of class Interface



/// \relates Interface
template<class FT, class OpT=::Pimpact::OperatorBase<FT>, class IOpT=::Pimpact::OperatorBase<FT> >
Teuchos::RCP< Interface<FT,OpT,IOpT> > createInterface(
    Teuchos::RCP<FT> fu=Teuchos::null,
    Teuchos::RCP<OpT>  op=Teuchos::null,
    Teuchos::RCP<IOpT> jop=Teuchos::null,
    Teuchos::RCP<FT> sol=Teuchos::null	) {

	return( Teuchos::rcp( new Interface<FT,OpT,IOpT>(fu,op,jop,sol) ) );
}


} // end of namespace Pimpact
} // end of namespace NOX


#endif // end of #ifndef NOX_PIMPACT_INTERFACE_HPP
