#pragma once
#ifndef PIMPACT_COMPOUNDOP_HPP
#define PIMPACT_COMPOUNDOP_HPP

#include "Teuchos_RCP.hpp"

//#include "Pimpact_ScalarField.hpp"
//
//#include "Pimpact_Operator.hpp"
//#include "Pimpact_OperatorBase.hpp"


namespace Pimpact {


/// \brief combines two operators.
///
/// Both operators are applied and the result is added.
/// the \c DomainFieldT of \c OP2 has to be equal to the \c RangeFieldT of \c OP1.
/// both operators should have the same OpType(OpType is obsolete)
/// \ingroup Operator
template< class OP1, class OP2 >
class AddOp {
public:
//  typedef MultiField<ModeField<VectorField<Scalar,Ordinal> > > MVF;
//  typedef typename MVF::Scalar S;
//  typedef typename MVF::Ordinal O;
  typedef typename OP1::DomainFieldT DomainFieldT;
  typedef typename OP2::RangeFieldT  RangeFieldT;
//  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
//  typedef OperatorBase<MVF> OP;
  typedef NonModeOp OpType;
private:
  Teuchos::RCP<OP1> op1_;
  Teuchos::RCP<OP2> op2_;
  Teuchos::RCP<typename OP1::RangeFieldT> temp_;

public:
  AddOp():
        op1_(Teuchos::null),
        op2_(Teuchos::null),
        temp_(Teuchos::null)
       {};

  AddOp(
      const Teuchos::RCP<OP1>& op1,
      const Teuchos::RCP<OP2>& op2,
      const Teuchos::RCP<DomainFieldT> temp ):
        op1_(op1), op2_(op2),
        temp_(temp)
        {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    op1_->apply( x, *temp_ );
    op2_->apply( x, y );
    y.add( 1., *temp_, 1., y );
  }

  void assignField( const DomainFieldT& mv ) {
    op1_->assignField( mv );
    op2_->assignField( mv );
  };

  bool hasApplyTranspose() const { return( op1_->hasApplyTranspose() && op2_->hasApplyTranspose() ); }

}; // end of class AddOp



/// \relates AddOp
template<class OP1, class OP2 >
Teuchos::RCP< AddOp<OP1, OP2> > createAddOp(
    const Teuchos::RCP<OP1>& op1,
    const Teuchos::RCP<OP2>& op2,
    const Teuchos::RCP<typename OP1::DomainFieldT>& temp
      ) {
  return( Teuchos::rcp( new AddOp<OP1,OP2>( op1, op2, temp ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDOP_HPP
