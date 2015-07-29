#pragma once
#ifndef PIMPACT_INVERSEOP_HPP
#define PIMPACT_INVERSEOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_LinearProblem.hpp"

#include "Pimpact_LinSolverParameter.hpp"


namespace Pimpact{



/// \ingroup Operator
/// hides all the linear solver typeerasiure stuff
/// \tparam OperatorT has to be of type \c Pimpact::MultiWrapOp
template< class OT >
class InverseOp {

public:

  typedef OT OperatorT;
  typedef typename OperatorT::SpaceT SpaceT;

  typedef typename OperatorT::DomainFieldT MF;

  typedef OperatorBase<MF>   Op;

  typedef typename OperatorT::DomainFieldT DomainFieldT;
  typedef typename OperatorT::RangeFieldT  RangeFieldT;

protected:

  Teuchos::RCP< LinearProblem<MF> > linprob_;

public:

	template<class IOperatorT>
	InverseOp( const Teuchos::RCP<IOperatorT>& op ) {

		auto para = 
			createLinSolverParameter("GMRES",1.e-12,-1, Teuchos::rcp( &std::cout, false ), 200 );
		linprob_ = createLinearProblem<MF>(
				createOperatorBase( create<OperatorT>(op) ),
				create<MF>( op->space() ),
				create<MF>( op->space() ),
				para,
				"GMRES" );
	}

 template<class IOperatorT>
 InverseOp( const Teuchos::RCP<IOperatorT>& op,
		 Teuchos::RCP<Teuchos::ParameterList> pl ):
	 linprob_( createLinearProblem<MF>(
				 createOperatorBase( create<OperatorT>(op) ),
				 create<MF>( op->space() ),
				 create<MF>( op->space() ),
				 Teuchos::rcpFromRef( pl->sublist("Solver") ), 
				 pl->get<std::string>("Solver name","GMRES") ) ) { }


  void apply( const MF& x, MF& y ) const {
    linprob_->solve( Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x) );
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

	const std::string getLabel() const { return( linprob_->getProblem()->getOperator()->getLabel() + std::string("^-1 ")  ); };

  void print( std::ostream& out=std::cout ) const {
		out << "Inverse:\n";
    linprob_->print( out );
  }

}; // end of class InverseOp



/// \relates InverseOp
template< class OpT>
Teuchos::RCP< InverseOp<OpT> >
createInverseOp(
    const Teuchos::RCP<OpT>& op, Teuchos::RCP< Teuchos::ParameterList > pl=Teuchos::null ) {

	if( pl.is_null() )
    return( Teuchos::rcp( new InverseOp<OpT>( op ) ) );
	else
    return( Teuchos::rcp( new InverseOp<OpT>( op, pl ) ) );

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_INVERSEOP_HPP
