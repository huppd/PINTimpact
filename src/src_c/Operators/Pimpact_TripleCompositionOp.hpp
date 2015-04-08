#pragma once
#ifndef PIMPACT_TRIPLECOMPOSITIONOP_HPP
#define PIMPACT_TRIPLECOMPOSITIONOP_HPP


#include "Teuchos_RCP.hpp"


namespace Pimpact {


/// \brief make the composition of two operators.
///
/// \f[ op1 op2 op3 \mathbf{x} = \mathbf{y} \f]
/// Both operators are applied sequentially.
/// the \c DomainFieldT \c OP1 has to equal to the \c RangeFieldT \c OP2.
/// \todo insert static_assert
/// \ingroup Operator
template< class OP1, class OP2, class OP3=OP2 >
class TripleCompositionOp {

public:

  typedef typename OP3::DomainFieldT DomainFieldT;
  typedef typename OP1::RangeFieldT  RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP<typename OP2::RangeFieldT> temp1_; // has to be equal to OP2::DomainFieldT
  Teuchos::RCP<typename OP2::DomainFieldT> temp2_; // has to be euqal to OP3::DomainFieldT

  Teuchos::RCP<OP1> op1_;
  Teuchos::RCP<OP2> op2_;
  Teuchos::RCP<OP3> op3_;

public:

  TripleCompositionOp(
      const Teuchos::RCP<OP1>&          op1,
      const Teuchos::RCP<OP2>&          op2,
      const Teuchos::RCP<OP3>&          op3
      ):
        temp1_( create<typename OP2::DomainFieldT>(op1->space()) ),
        temp2_( create<typename OP2::RangeFieldT>(op1->space()) ),
        op1_(op1),
        op2_(op2),
        op3_(op3)
{};

  void apply(const DomainFieldT& x, RangeFieldT& y, Belos::ETrans trans=Belos::NOTRANS) const {

//		temp1_->init(0.);
//		x.level();

    op3_->apply( x, *temp1_);

    op2_->apply( *temp1_, *temp2_ );

//		y.init(0.);
//		temp2_->level();

    op1_->apply( *temp2_, y );

  }

  /// \note here nothing happens, because it is assumed to be done somewhere else
  void assignField( const RangeFieldT& field ) {
//    if( !op1_.is_null() )
//      op1_->assignField( field );
//    if( !op2_.is_null() )
//      op2_->assignField( field );
//    if( !op2_.is_null() )
//      op2_->assignField( field );
  };
//  void assignField( const DomainFieldT& field ) {
//    if( !op1_.is_null() )
//      op1_->assignField( field );
//    if( !op2_.is_null() )
//      op2_->assignField( field );
//    if( !op2_.is_null() )
//      op2_->assignField( field );
//  };

	Teuchos::RCP<const SpaceT> space() const { return(op1_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
	 if( !op1_.is_null() ) op1_->setParameter( para );
	 if( !op2_.is_null() ) op2_->setParameter( para );
	 if( !op3_.is_null() ) op3_->setParameter( para );
	}

  bool hasApplyTranspose() const { return( op1_->hasApplyTranspose() && op2_->hasApplyTranspose() && op3_->hasApplyTranspose() ); }

}; // end of class TripleCompositionOp



/// \relates TripleCompositionOp
template<class OP1, class OP2, class OP3>
Teuchos::RCP< TripleCompositionOp<OP1, OP2, OP3> > createTripleCompositionOp(
    const Teuchos::RCP<OP1>& op1,
    const Teuchos::RCP<OP2>& op2,
    const Teuchos::RCP<OP3>& op3
      ) {

  return( Teuchos::rcp(	new TripleCompositionOp<OP1,OP2,OP3>( op1, op2, op3 ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRIPLECOMPOSITIONOP_HPP
