#pragma once
#ifndef PIMPACT_TRIPLECOMPOSITIONOP_HPP
#define PIMPACT_TRIPLECOMPOSITIONOP_HPP


#include "Teuchos_RCP.hpp"


namespace Pimpact {


/// \brief make the composition of two operators.
///
/// Both operators are applied sequentially.
/// the \c DomainFieldT \c OP1 has to equal to the \c RangeFieldT \c OP2.
/// \todo insert static_assert
/// \todo remove dirty hack with MultiField casting
/// \ingroup Operator
template< class OP1, class OP2, class OP3=OP2 >
class TripleCompositionOp {

public:

  typedef typename OP1::DomainFieldT DomainFieldT;
  typedef typename OP3::RangeFieldT  RangeFieldT;

private:

  Teuchos::RCP<typename OP1::RangeFieldT> temp1_; // has to be equal to OP2::DomainFieldT
  Teuchos::RCP<typename OP3::DomainFieldT> temp2_; // has to be euqal to OP3::DomainFieldT

  Teuchos::RCP<OP1> op1_;
  Teuchos::RCP<OP2> op2_;
  Teuchos::RCP<OP3> op3_;

public:

  TripleCompositionOp(
      const Teuchos::RCP<typename OP1::RangeFieldT>& temp1=Teuchos::null,
      const Teuchos::RCP<typename OP3::DomainFieldT>& temp2=Teuchos::null,
      const Teuchos::RCP<OP1>&          op1=Teuchos::null,
      const Teuchos::RCP<OP2>&          op2=Teuchos::null,
      const Teuchos::RCP<OP3>&          op3=Teuchos::null
      ):
        temp1_(temp1),temp2_(temp2),
        op1_(op1), op2_(op2), op3_(op3)
{};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    op1_->apply( x, *temp1_);
    op2_->apply( *Pimpact::createMultiField(temp1_), *Pimpact::createMultiField(temp2_) );
    op3_->apply( *temp2_, y );
  }

  void assignField( const RangeFieldT& field ) {
//    if( !op1_.is_null() )
//      op1_->assignField( field );
//    if( !op2_.is_null() )
//      op2_->assignField( field );
//    if( !op2_.is_null() )
//      op2_->assignField( field );
  };

  bool hasApplyTranspose() const { return( op1_->hasApplyTranspose() && op2_->hasApplyTranspose() && op3_->hasApplyTranspose() ); }

}; // end of class TripleCompositionOp



/// \relates TripleCompositionOp
template<class OP1, class OP2, class OP3>
Teuchos::RCP< TripleCompositionOp<OP1, OP2, OP3> > createTripleCompositionOp(
    const Teuchos::RCP<typename OP1::RangeFieldT>& temp1=Teuchos::null,
    const Teuchos::RCP<typename OP3::DomainFieldT>& temp2=Teuchos::null,
    const Teuchos::RCP<OP1>& op1=Teuchos::null,
    const Teuchos::RCP<OP2>& op2=Teuchos::null,
    const Teuchos::RCP<OP3>& op3=Teuchos::null
      ) {
  return( Teuchos::rcp( new TripleCompositionOp<OP1,OP2,OP3>( temp1, temp2, op1, op2, op3 ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRIPLECOMPOSITIONOP_HPP
