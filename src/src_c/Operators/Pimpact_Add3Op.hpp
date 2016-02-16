#pragma once
#ifndef PIMPACT_ADD3OP_HPP
#define PIMPACT_ADD3OP_HPP


#include "Teuchos_RCP.hpp"

#include "BelosTypes.hpp"




namespace Pimpact {



/// \brief combines three operators.
///
/// Both operators are applied and the result is added.
/// the \c DomainFieldT  \c RangeFieldT has to be equal for both \c OP1, \c OP2, \c OP3.
/// \ingroup Operator
template< class OP1, class OP2, class OP3 >
class Add3Op {

public:

  using DomainFieldT = typename OP1::DomainFieldT;
  using RangeFieldT = typename OP1::RangeFieldT;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  Teuchos::RCP<RangeFieldT> temp_;

  Teuchos::RCP<OP1> op1_;
  Teuchos::RCP<OP2> op2_;
  Teuchos::RCP<OP3> op3_;


public:

  Add3Op(
      const Teuchos::RCP<OP1>& op1,
      const Teuchos::RCP<OP2>& op2=Teuchos::null,
      const Teuchos::RCP<OP3>& op3=Teuchos::null ):
        temp_(op1->space()),
        op1_(op1), op2_(op2), op3_(op3) {};

  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {
    if( !op1_.is_null() ) {
      op1_->apply( x, y  );
      if( !op2_.is_null() ) {
        op2_->apply( x, *temp_ );
        y.add( 1., y, 1., *temp_ );
      }
      if( !op3_.is_null() ) {
        op3_->apply( x, *temp_ );
        y.add( 1., y, 1., *temp_ );
      }
    }
    else if( !op2_.is_null() ) {
      op2_->apply( x, y );
      if( !op3_.is_null() ) {
        op3_->apply( x, *temp_ );
        y.add( 1., y, 1., *temp_ );
      }
    }
  }

  void assignField( const DomainFieldT& field ) {
    if( !op1_.is_null() )
      op1_->assignField( field );
    if( !op2_.is_null() )
      op2_->assignField( field );
    if( !op3_.is_null() )
      op3_->assignField( field );
  };

	Teuchos::RCP<const SpaceT> space() const { return(op1_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {
    if( !op1_.is_null() ) op1_->setParameter( para );
    if( !op2_.is_null() ) op2_->setParameter( para );
    if( !op3_.is_null() ) op3_->setParameter( para );
	}

  bool hasApplyTranspose() const { return( op1_->hasApplyTranspose() && op2_->hasApplyTranspose() && op3_->hasApplyTranspose() ); }

	const std::string getLabel() const { return( op1_->getLabel() + std::string(" + ") + op2_->getLabel() + std::string(" + ") + op3_->getLabel() ); };

}; // end of class Add3Op



/// \relates Add3Op
template<class OP1, class OP2, class OP3=OP1 >
Teuchos::RCP< Add3Op<OP1, OP2, OP3> > createAdd3Op(
    const Teuchos::RCP<OP1>& op1,
    const Teuchos::RCP<OP2>& op2=Teuchos::null,
    const Teuchos::RCP<OP3>& op3=Teuchos::null ) {

  return( Teuchos::rcp( new Add3Op<OP1,OP2,OP3>( op1, op2, op3 ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_ADD3OP_HPP
