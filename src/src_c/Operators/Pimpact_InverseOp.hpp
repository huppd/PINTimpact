#pragma once
#ifndef PIMPACT_INVERSEOP_HPP
#define PIMPACT_INVERSEOP_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosSolverManager.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_MultiOpWrap.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// hides all the linear solver typeerasiure stuff
/// 
/// \todo merge with Prec...
/// \ingroup Operator
template< class OpT >
class InverseOp {

public:

	using OperatorT = OpT;

	using SpaceT = typename OperatorT::SpaceT;

	using DomainFieldT = typename OperatorT::DomainFieldT;
	using RangeFieldT = typename OperatorT::RangeFieldT;

	using MF = typename Pimpact::MultiField<DomainFieldT>;

	using MOpT = OperatorBase<MF>;

protected:

	using ST = typename SpaceT::Scalar;

	bool level_;
	bool levelRHS_;
	bool nullspaceOrtho_;
	bool initZero_;
	bool debug_;

	ST relTol_;

	std::string solverName_;

	Teuchos::RCP< Teuchos::ParameterList > solverParameter_;
	Teuchos::RCP< Belos::LinearProblem<ST, MF, MOpT> > problem_;

	Teuchos::RCP<const RangeFieldT> nullspace_;

public:

	/// should be avoided( used in peri_navierSMGXML)
	InverseOp( const Teuchos::RCP<const SpaceT>& space ):
		level_(false),
		levelRHS_(false),
		nullspaceOrtho_(false),
		initZero_(false),
		debug_(false),
		relTol_(1.),
		solverName_( "GMRES" ),
		solverParameter_( createLinSolverParameter("GMRES",1.e-1,-1, Teuchos::rcp( new Teuchos::oblackholestream ), 200 ) ),
		problem_( Teuchos::rcp( new Belos::LinearProblem<ST,MF,MOpT> () )) {}


	InverseOp( const Teuchos::RCP<OpT>& op ):
		level_(false),
		levelRHS_(false),
		nullspaceOrtho_(false),
		initZero_(false),
		debug_(false),
		relTol_(1.),
		solverName_( "GMRES" ),
		solverParameter_( createLinSolverParameter( "GMRES", 1.e-1, -1, Teuchos::rcp( new Teuchos::oblackholestream ),10 ) ),
		problem_(
				Teuchos::rcp(
					new Belos::LinearProblem<ST,MF,MOpT>(
						createMultiOperatorBase( op ),
						Teuchos::rcp( new MF(op->space()) ),
						Teuchos::rcp( new MF(op->space()) ) ) ) ) {}


	//template<class IOperatorT>
	InverseOp( const Teuchos::RCP<OperatorT>& op,
			const Teuchos::RCP<Teuchos::ParameterList>& pl ):
		level_( pl->get<bool>( "level", false ) ),
		levelRHS_( pl->get<bool>( "level RHS", false ) ),
		nullspaceOrtho_( pl->get<bool>( "nullspace ortho", false ) ),
		initZero_( pl->get<bool>( "initZero", false ) ),
		debug_( pl->get<bool>( "debug", false ) ),
		relTol_( pl->get<ST>( "relative Tol", 1. ) ),
		solverName_( pl->get<std::string>("Solver name", "GMRES") ),
		solverParameter_( Teuchos::sublist( pl, "Solver" ) ),
		problem_(
				Teuchos::rcp( new Belos::LinearProblem<ST,MF,MOpT>(
						createMultiOperatorBase( op ),
						Teuchos::rcp( new MF(op->space()) ),
						Teuchos::rcp( new MF(op->space()) ) ) )) {}


	void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N  ) const {
		apply(
				*wrapMultiField( Teuchos::rcpFromRef<DomainFieldT>(
						const_cast<DomainFieldT&>(x) ) ),
				*wrapMultiField( Teuchos::rcpFromRef<RangeFieldT>(y) ) );
	}


	void apply( const MF& x, MF& y, const Add& add=Add::N  ) const {

		if( levelRHS_ ) { x.level(); }
		if( initZero_ ) { y.init( ); }
		if( nullspaceOrtho_ ) {
			for( int i=0; i<x.getNumberVecs(); ++i ) {
				ST bla = -nullspace_->dot( x.getField(i) );
				if( 0==space()->rankST() ) std::cout << getLabel()<< ": nullspace contributtion[" << i<< "]: " << std::abs(bla)  << "\n";
				if( std::abs( bla ) >= Teuchos::ScalarTraits<ST>::eps() )
					const_cast<MF&>(x).getField(i).add( 1., x.getField(i), bla, *nullspace_ );
			}
			if( 0==space()->rankST() ) std::cout << "\n";
		}

		problem_->setProblem( Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x) );
		Belos::SolverFactory<ST,MF,MOpT> factory;

		Teuchos::RCP< Belos::SolverManager<ST,MF,MOpT> > solver =
			factory.create(solverName_,solverParameter_);

		solver->setProblem( problem_ );
		Belos::ReturnType succes = solver->solve();

		if( debug_ && Belos::ReturnType::Unconverged==succes ) {
			x.write();
			assert( Belos::ReturnType::Converged==succes );
			//TEUCHOS_TEST_FOR_EXCEPT( false );
			//TEUCHOS_TEST_FOR_EXCEPT( true );
		}

		if( level_ ) y.level();
	}


	void assignField( const DomainFieldT& mvr ) {

		auto mv = wrapMultiField(
				Teuchos::rcpFromRef<DomainFieldT>( const_cast<DomainFieldT&>(mvr) ) );

		Teuchos::rcp_const_cast<MOpT>( problem_->getOperator() )->assignField( *mv );

		if( problem_->isLeftPrec() ) {
			auto opPrec = Teuchos::rcp_const_cast<MOpT>( problem_->getLeftPrec() );
			opPrec->assignField( *mv );
		}

		if( problem_->isRightPrec() ) {
			auto opPrec = Teuchos::rcp_const_cast<MOpT>( problem_->getRightPrec() );
			opPrec->assignField( *mv );
		}
	};


	void setNullspace( const Teuchos::RCP<const RangeFieldT>& nullspace ) {
		nullspaceOrtho_ = true;
		nullspace_ = nullspace;
	}


	constexpr const Teuchos::RCP<const SpaceT>& space() {
		return( problem_->getOperator()->space() );
	};


	constexpr Teuchos::RCP<const MOpT> getOperator() {
		return( problem_->getOperator() );
	};


	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {

		if( para->name()=="Linear Solver" )
			solverParameter_->set<ST>( "Convergence Tolerance",
					para->get<double>("Convergence Tolerance")*relTol_ );

		Teuchos::rcp_const_cast<MOpT>( problem_->getOperator() )->setParameter( para );

		if( problem_->isLeftPrec() )
			Teuchos::rcp_const_cast<MOpT>( problem_->getLeftPrec() )->setParameter( para );

		if( problem_->isRightPrec() )
			Teuchos::rcp_const_cast<MOpT>( problem_->getRightPrec() )->setParameter( para );
	}


  /// \brief Set left preconditioner (\c LP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setLeftPrec( const Teuchos::RCP<const MOpT>& LP ) {
    problem_->setLeftPrec( LP );
  }

  /// \brief Set right preconditioner (\c RP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setRightPrec( const Teuchos::RCP<const MOpT>& RP ) {
    problem_->setRightPrec( RP );
  }

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( problem_->getOperator()->getLabel() + std::string("^-1 ")  ); };

  void print( std::ostream& out=std::cout ) const {
		out << "Inverse:\n";
		problem_->getOperator()->print( out );
  }


}; // end of class InverseOp



/// \relates InverseOp
template<class OpT>
Teuchos::RCP< InverseOp<OpT> >
createInverseOp( const Teuchos::RCP<OpT>& op ) {
	return( Teuchos::rcp( new InverseOp<OpT>( op ) ) );
}


/// \relates InverseOp
template<class OpT>
Teuchos::RCP< InverseOp<OpT> >
createInverseOp(
		const Teuchos::RCP<OpT>& op,
		const Teuchos::RCP<Teuchos::ParameterList>& pl ) {

	return( Teuchos::rcp( new InverseOp<OpT>( op, pl ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSEOP_HPP
