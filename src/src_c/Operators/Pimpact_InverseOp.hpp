#pragma once
#ifndef PIMPACT_INVERSEOP_HPP
#define PIMPACT_INVERSEOP_HPP


#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_MultiOpWrap.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{



/// \ingroup Operator
/// hides all the linear solver typeerasiure stuff
/// \tparam OperatorT has to be of type \c Pimpact::MultiWrapOp
template< class MOT >
class InverseOp {

public:

  using OperatorT = MOT;

  using SpaceT = typename OperatorT::SpaceT;

  using MF = typename OperatorT::DomainFieldT;

  using Op = OperatorBase<MF>;

  using DomainFieldT = typename OperatorT::DomainFieldT;
  using RangeFieldT = typename OperatorT::RangeFieldT;

protected:

	using ScalarT = typename SpaceT::Scalar;
	bool level_;
	bool levelRHS_;
	bool nullspaceOrtho_;
	bool initZero_;
  Teuchos::RCP< LinearProblem<MF> > linprob_;
	Teuchos::RCP<const RangeFieldT> nullspace_;

public:

	InverseOp( const Teuchos::RCP<const SpaceT>& space ):
		level_(false),
		levelRHS_(false),
		nullspaceOrtho_(false),
		initZero_(false) {

			auto para = 
				createLinSolverParameter("GMRES",1.e-1,-1, Teuchos::rcp( new Teuchos::oblackholestream ), 200 );

			linprob_ = createLinearProblem<MF>(
					createOperatorBase( Teuchos::rcp( new OperatorT(space) ) ),
					create<MF>(space),
					create<MF>(space),
					para,
					"GMRES" );
	}

	template<class IOperatorT>
	InverseOp( const Teuchos::RCP<IOperatorT>& op ):
		level_(false),
		levelRHS_(false),
		nullspaceOrtho_(false),
		initZero_(false) {

			auto para = 
				createLinSolverParameter("GMRES",1.e-1,-1, Teuchos::rcp( new Teuchos::oblackholestream ), 10 );
			linprob_ = createLinearProblem<MF>(
					createOperatorBase( create<OperatorT>(op) ),
					create<MF>( op->space() ),
					create<MF>( op->space() ),
					para,
					"GMRES" );
		}

 template<class IOperatorT>
 InverseOp( const Teuchos::RCP<IOperatorT>& op,
		 const Teuchos::RCP<Teuchos::ParameterList>& pl ):
	 level_( pl->get<bool>( "level", false ) ),
	 levelRHS_( pl->get<bool>( "level RHS", false ) ),
	 nullspaceOrtho_( pl->get<bool>( "nullspace ortho", false ) ),
	 initZero_( pl->get<bool>( "initZero", false ) ),
	 linprob_( createLinearProblem<MF>(
				 createOperatorBase( create<OperatorT>(op) ),
				 create<MF>( op->space() ),
				 create<MF>( op->space() ),
				 Teuchos::rcpFromRef( pl->sublist("Solver") ), 
				 pl->get<std::string>("Solver name","GMRES") ) ) { }


 /// \todo clean up: rm clone make nullspace normalized.
 void apply( const MF& x, MF& y ) const {
	 if( levelRHS_ ) { x.level(); }
	 if( initZero_ ) { y.init( ); }
	 if( nullspaceOrtho_ ) {

		 //auto x_ = x.clone();
		 for( int i=0; i<x.getNumberVecs(); ++i ) {
			 ScalarT bla = -nullspace_->getFieldPtr(0)->dot( x.getField(i) );
			 std::cout << "bla: " << std::setprecision(10) << bla  << "\n";
			 //std::cout << "x_->getNumberVecs(): " << x_->getNumberVecs() << "\n";
			 //std::cout << "nullspace_->getNumberVecs(): " << nullspace_->getNumberVecs() << "\n";
			 const_cast<MF&>(x).getFieldPtr(i)->add( 1., x.getField(i), bla, nullspace_->getField(0) );
		 }

		 //linprob_->solve( Teuchos::rcpFromRef(y), x_ );
	 }
	 //else
    linprob_->solve( Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x) );
	 if( level_    ) { y.level(); }
 }


 void assignField( const DomainFieldT& mv ) {

	 auto prob = linprob_->getProblem();

	 Teuchos::rcp_const_cast<Op>( prob->getOperator() )->assignField( mv );

	 if( prob->isLeftPrec() ) {
		 auto opPrec = Teuchos::rcp_const_cast<Op>( prob->getLeftPrec() );
		 opPrec->assignField( mv );
	 }

	 if( prob->isRightPrec() ) {
		 auto opPrec = Teuchos::rcp_const_cast<Op>( prob->getRightPrec() );
		 opPrec->assignField( mv );
	 }
 }

 void setNullspace( const Teuchos::RCP<const RangeFieldT>& nullspace ) {
	 nullspaceOrtho_ = true;
	 nullspace_ = nullspace;
 }

 constexpr const Teuchos::RCP<const SpaceT>& space() const { return(linprob_->space()); };

 constexpr Teuchos::RCP<const LinearProblem<MF> > getLinearProblem() const {
	 return(linprob_);
 }

 constexpr Teuchos::RCP<const Op> getOperatorPtr() const {
	 return( getLinearProblem()->getOperatorPtr() );
 };

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
    linprob_->setLeftPrec( LP );
  }

  /// \brief Set right preconditioner (\c RP) of linear problem \f$AX = B\f$.
  ///
  /// The operator is set by pointer; no copy of the operator is made.
  void setRightPrec(const Teuchos::RCP<const OperatorBase<MF> > &RP) {
    linprob_->setRightPrec( RP );
  }

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( linprob_->getProblem()->getOperator()->getLabel() + std::string("^-1 ")  ); };

  void print( std::ostream& out=std::cout ) const {
		out << "Inverse:\n";
    linprob_->print( out );
  }


}; // end of class InverseOp



/// \relates InverseOp
template< class OpT>
Teuchos::RCP< InverseOp< MultiOpWrap<OpT> > >
createInverseOp( const Teuchos::RCP<OpT>& op ) {
	return( Teuchos::rcp( new InverseOp<MultiOpWrap<OpT> >( op ) ) );
}


/// \relates InverseOp
template< class OpT>
Teuchos::RCP< InverseOp< MultiOpWrap<OpT> > >
createInverseOp(
		const Teuchos::RCP<OpT>& op,
		const Teuchos::RCP<Teuchos::ParameterList>& pl ) {

	return( Teuchos::rcp( new InverseOp<MultiOpWrap<OpT> >( op, pl ) ) );
}

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSEOP_HPP
