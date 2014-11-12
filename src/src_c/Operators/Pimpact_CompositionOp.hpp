#pragma once
#ifndef PIMPACT_COMPOSITIONOP_HPP
#define PIMPACT_COMPOSITIONOP_HPP


#include "Teuchos_RCP.hpp"




namespace Pimpact {



/// \brief make the composition of two operators.
///
/// Both operators are applied sequentially.
/// the \c DomainFieldT \c OP2 has to equal to the \c RangeFieldT \c OP1.
/// \ingroup Operator
template< class OP1, class OP2 >
class CompositionOp {

public:

  typedef typename OP2::DomainFieldT DomainFieldT;
  typedef typename OP1::RangeFieldT  RangeFieldT;

  typedef typename OP2::RangeFieldT TempFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP<TempFieldT> temp_;

  Teuchos::RCP<OP1> op1_;
  Teuchos::RCP<OP2> op2_;

public:

  CompositionOp(
      const Teuchos::RCP<TempFieldT>& temp=Teuchos::null,
      const Teuchos::RCP<OP1>&           op1=Teuchos::null,
      const Teuchos::RCP<OP2>&           op2=Teuchos::null ):
        temp_(temp),
        op1_(op1), op2_(op2)//, op3_(op3)
{};

  void apply(
      const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {

    if( op1_.is_null() ) {
//      op2_->apply( x, y );
    }
    else {
      op2_->apply( x, *temp_ );
      op1_->apply( *temp_, y );
    }
  }

  void assignField( const DomainFieldT& field ) {
    if( !op1_.is_null() )
//      op1_->assignField( field );
    if( !op2_.is_null() )
      op2_->assignField( field );
  };

  bool hasApplyTranspose() const { return( op1_->hasApplyTranspose() && op2_->hasApplyTranspose() /*&& op3_->hasApplyTranspose()*/ ); }

}; // end of class CompositionOp



/// \relates CompositionOp
template<class OP1, class OP2>
Teuchos::RCP< CompositionOp<OP1, OP2> > createCompositionOp(
    const Teuchos::RCP<typename OP2::RangeFieldT>& temp=Teuchos::null,
    const Teuchos::RCP<OP1>& op1=Teuchos::null,
    const Teuchos::RCP<OP2>& op2=Teuchos::null
      ) {
  return( Teuchos::rcp( new CompositionOp<OP1,OP2>( temp, op1, op2 ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOSITIONOP_HPP
