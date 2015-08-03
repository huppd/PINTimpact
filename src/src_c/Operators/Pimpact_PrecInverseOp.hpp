#pragma once
#ifndef PIMPACT_PRECINVERSEOP_HPP
#define PIMPACT_PRECINVERSEOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_LinearProblem.hpp"

#include "Pimpact_MultiField.hpp"
#include "Pimpact_LinSolverParameter.hpp"


namespace Pimpact{



/// \ingroup Operator
/// hides all the linear solver typeerasiure stuff
/// \tparam OperatorT has to be of type \c Pimpact::MultiWrapOp
template< class OT, template<class> class PT >
class PrecInverseOp {

public:

	typedef OT OperatorT;
	typedef PT<OT> PreconditionerT;

	typedef typename OperatorT::SpaceT SpaceT;


	typedef typename OperatorT::DomainFieldT DomainFieldT;
	typedef typename OperatorT::RangeFieldT  RangeFieldT;

	typedef typename Pimpact::MultiField<DomainFieldT> MF;

	typedef OperatorBase<MF> Op;

protected:

	Teuchos::RCP< LinearProblem<MF> > linprob_;

public:

	PrecInverseOp( const Teuchos::RCP<OperatorT>& op ) { 

		auto para = 
			createLinSolverParameter("GMRES",1.e-12,-1, Teuchos::rcp( &std::cout, false ), 200 );

		linprob_ = createLinearProblem<MF>(
				createMultiOperatorBase( op ),
				create<MF>( op->space() ),
				create<MF>( op->space() ),
				para,
				"GMRES" );

		linprob_->setLeftPrec(
				createMultiOperatorBase(
					create<PreconditionerT>(op) ) );

	}


	PrecInverseOp( const Teuchos::RCP<OperatorT>& op,
			Teuchos::RCP<Teuchos::ParameterList> pl ):
		linprob_( createLinearProblem<MF>(
					createMultiOperatorBase( op ),
					create<MF>( op->space() ),
					create<MF>( op->space() ),
					Teuchos::rcpFromRef( pl->sublist("Solver") ), 
					pl->get<std::string>("Solver name","GMRES") ) ) { 

			if( pl->get<bool>( "LeftPrec", true ) )
				linprob_->setLeftPrec(
						createMultiOperatorBase(
							create<PreconditionerT>(op, Teuchos::rcpFromRef( pl->sublist("Preconditioner") )) ) );
			else
				linprob_->setRightPrec(
						createMultiOperatorBase(
							create<PreconditionerT>(op, Teuchos::rcpFromRef( pl->sublist("Preconditioner") )) ) );

		}


	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		linprob_->solve(
				createMultiField( Teuchos::rcpFromRef(y) ),
				createMultiField( Teuchos::rcpFromRef( const_cast<DomainFieldT&>(x) ) ) );

	}


	void assignField( const DomainFieldT& mvr ) {

		auto mv =
			createMultiField(
					Teuchos::rcpFromRef<DomainFieldT>( const_cast<DomainFieldT&>(mvr) ) );

		auto prob = linprob_->getProblem();

		Teuchos::rcp_const_cast<Op>( prob->getOperator() )->assignField( *mv );

		if( prob->isLeftPrec() ) {
			auto opPrec = Teuchos::rcp_const_cast<Op>( prob->getLeftPrec() );
			opPrec->assignField( mv );
		}

		if( prob->isRightPrec() ) {
			auto opPrec = Teuchos::rcp_const_cast<Op>( prob->getRightPrec() );
			opPrec->assignField( mv );
		}

	};


	Teuchos::RCP<const SpaceT> space() const { return(linprob_->space()); };

	Teuchos::RCP< LinearProblem<MF> > getLinearProblem() { return(linprob_); }


	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		auto prob = linprob_->getProblem();

		Teuchos::rcp_const_cast<Op>( prob->getOperator() )->setParameter( para );

		if( prob->isLeftPrec() ) {
			auto opPrec = Teuchos::rcp_const_cast<Op>( prob->getLeftPrec() );
			opPrec->setParameter( para );
		}

		if( prob->isRightPrec() ) {
			auto opPrec = Teuchos::rcp_const_cast<Op>( prob->getRightPrec() );
			opPrec->setParameter( para );
		}
	}


	/// \brief Set left preconditioner (\c LP) of linear problem \f$AX = B\f$.
	///
	/// The operator is set by pointer; no copy of the operator is made.
	void setLeftPrec(const Teuchos::RCP<const OperatorBase<MF> > &LP) {
		linprob_->getProblem()->setLeftPrec( LP );
	}

	/// \brief Set right preconditioner (\c RP) of linear problem \f$AX = B\f$.
	///
	/// The operator is set by pointer; no copy of the operator is made.
	void setRightPrec(const Teuchos::RCP<const OperatorBase<MF> > &RP) {
		linprob_->getProblem()->setRightPrec( RP );
	}

	bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const {
		return( linprob_->getProblem()->getOperator()->getLabel() + std::string("^-1(prec) ")  );
	};

	void print( std::ostream& out=std::cout ) const {
		out << "PrecInverse:\n";
		linprob_->print( out );
	}

}; // end of class PrecInverseOp



///// \relates PrecInverseOp
//template< class OpT, template<class> class PT >
//Teuchos::RCP< PrecInverseOp<OpT,PT> >
//createPrecInverseOp(
//		const Teuchos::RCP<OpT>& op, Teuchos::RCP< Teuchos::ParameterList > pl=Teuchos::null ) {
//
//	if( pl.is_null() )
//		return( Teuchos::rcp( new PrecInverseOp<OpT,PT>( op ) ) );
//	else
//		return( Teuchos::rcp( new PrecInverseOp<OpT,PT>( op, pl ) ) );
//
//}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_PRECINVERSEOP_HPP
