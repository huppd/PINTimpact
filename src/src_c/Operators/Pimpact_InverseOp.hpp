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

  typedef typename OperatorT::DomainFieldT  MF;

  typedef OperatorBase<MF>     Op;

  typedef typename OperatorT::DomainFieldT DomainFieldT;
  typedef typename OperatorT::RangeFieldT  RangeFieldT;

protected:

  Teuchos::RCP< LinearProblem<MF> > linprob_;

public:

  InverseOp( const Teuchos::RCP<OperatorT>& op ) {
    linprob_ = createLinearProblem<MF>(
        createOperatorBase(op),
        create<MF>( op->getSpace() ),
        create<MF>( op->getSpace() ),
        createLinSolverParameter("GMRES",1.e-12,-1),
        "GMRES" );
  }

  template<class IOperatorT>
  InverseOp( const Teuchos::RCP<IOperatorT>& op ) {
    linprob_ = createLinearProblem<MF>(
        createOperatorBase( create<OperatorT>(op) ),
        create<MF>( op->getSpace() ),
        create<MF>( op->getSpace() ),
//				createLinSolverParameter("GMRES",1.e-12,-1),
			 Teuchos::parameterList(),
        "GMRES" );
  }


  void apply( const MF& x, MF& y, Belos::ETrans trans=Belos::NOTRANS ) const {
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
	

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		Teuchos::rcp_const_cast<Op>( linprob_->getProblem()->getOperator() )->setParameter( para );
	}


  bool hasApplyTranspose() const { return( false ); }


}; // end of class InverseOp



/// \relates InverseOp
template< class OpT>
Teuchos::RCP< InverseOp<OpT> >
createInverseOp(
    const Teuchos::RCP<OpT>& op ) {

    return( Teuchos::rcp( new InverseOp<OpT>( op ) ) );

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_INVERSEOP_HPP
