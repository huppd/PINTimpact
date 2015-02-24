#pragma once
#ifndef PIMPACT_INVERSEOPERATOR_HPP
#define PIMPACT_INVERSEOPERATOR_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_LinearProblem.hpp"




namespace Pimpact{



/// \ingroup Operator
/// \tparam FieldT has to be of type \c Pimpact::MultiField
/// \todo merge with InverseOp
/// \deprecated
template< class FieldT >
class InverseOperator {

  typedef FieldT               MF;
  typedef typename MF::Scalar  S;
  typedef OperatorBase<MF>     Op;

public:

  typedef FieldT  DomainFieldT;
  typedef FieldT  RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP< LinearProblem<MF> > linprob_;

public:

  InverseOperator( const Teuchos::RCP< LinearProblem<MF> >& linprob=Teuchos::null ):
    linprob_(linprob) {};

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

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

}; // end of class InverseOperator



/// \relates InverseOperator
template< class MF>
Teuchos::RCP< InverseOperator<MF> >
createInverseOperator(
    const Teuchos::RCP<LinearProblem<MF> >& linprob ) {

    return( Teuchos::rcp( new InverseOperator<MF>( linprob ) ) );

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_INVERSEOPERATOR_HPP
