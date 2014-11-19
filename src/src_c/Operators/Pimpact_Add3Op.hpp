#pragma once
#ifndef PIMPACT_ADD3OP_HPP
#define PIMPACT_ADD3OP_HPP


#include "Teuchos_RCP.hpp"


namespace Pimpact {




/// \brief combines three operators.
///
/// Both operators are applied and the result is added.
/// the \c DomainFieldT  \c RangeFieldT has to be equal for both \c OP1, \c OP2, \c OP3.
/// \ingroup Operator
template< class OP1, class OP2, class OP3 >
class Add3Op {

public:

  typedef typename OP1::DomainFieldT DomainFieldT;
  typedef typename OP1::RangeFieldT  RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP<RangeFieldT> temp_;

  Teuchos::RCP<OP1> op1_;
  Teuchos::RCP<OP2> op2_;
  Teuchos::RCP<OP3> op3_;


public:

  Add3Op(
      const Teuchos::RCP<DomainFieldT>& temp=Teuchos::null,
      const Teuchos::RCP<OP1>&           op1=Teuchos::null,
      const Teuchos::RCP<OP2>&           op2=Teuchos::null,
      const Teuchos::RCP<OP3>&           op3=Teuchos::null ):
        temp_(temp),
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

  bool hasApplyTranspose() const { return( op1_->hasApplyTranspose() && op2_->hasApplyTranspose() && op3_->hasApplyTranspose() ); }

}; // end of class Add3Op



/// \relates Add3Op
template<class OP1, class OP2, class OP3=OP1 >
Teuchos::RCP< Add3Op<OP1, OP2, OP3> > createAdd3Op(
    const Teuchos::RCP<typename OP1::DomainFieldT>& temp=Teuchos::null,
    const Teuchos::RCP<OP1>& op1=Teuchos::null,
    const Teuchos::RCP<OP2>& op2=Teuchos::null,
    const Teuchos::RCP<OP3>& op3=Teuchos::null
      ) {
  return( Teuchos::rcp( new Add3Op<OP1,OP2,OP3>( temp, op1, op2, op3) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_ADD3OP_HPP
