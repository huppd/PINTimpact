#pragma once
#ifndef PIMPACT_COMPOUNDOP_HPP
#define PIMPACT_COMPOUNDOP_HPP

#include "Teuchos_RCP.hpp"

//#include "Pimpact_ScalarField.hpp"
//
//#include "Pimpact_Operator.hpp"
//#include "Pimpact_OperatorBase.hpp"


namespace Pimpact {


/// \brief calls for combining two operators, which are apllied sequential.
///
/// the \c DomainFieldT of \c OP2 has to be equal to the \c RangeFieldT of \c OP1.
/// both operators should have the same OpType
template< class OP1, class OP2 >
class CompoundOp {
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
  Teuchos::RCP<DomainFieldT> temp_;
//  Teuchos::RCP<Div<S,O>> div_;
//  Teuchos::RCP<Grad<S,O>> grad_;
//  Teuchos::RCP<LinearProblem<S, MVF, OP > > op_;

public:
  CompoundOp():
        op1_(Teuchos::null),
        op2_(Teuchos::null),
        temp_(Teuchos::null)
       {};

  CompoundOp(
      const Teuchos::RCP<OP1>& op1,
      const Teuchos::RCP<OP2>& op2,
      const Teuchos::RCP<DomainFieldT> temp ):
        op1_(op1), op2_(op2),
        temp_(temp)
        {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
//    temp_ = const_cast<Teuchos::RCP<DomainFieldT> >(x.clone()) ;
//    temp_ = x.clone() ;
//    temp_ = Teuchos::rcp( new DomainFieldT() );
    op1_->apply( x, *temp_ );
    op2_->apply( *temp_, y );
  }

  bool hasApplyTranspose() const { return( op1_->hasApplyTranspose() && op2_->hasApplyTranspose() ); }

}; // end of class CompoundOp



template<class OP1, class OP2 >
Teuchos::RCP< CompoundOp<OP1, OP2> > createCompoundOp(
    const Teuchos::RCP<OP1>& op1,
    const Teuchos::RCP<OP2>& op2,
    const Teuchos::RCP<typename OP1::DomainFieldT>& temp
      ) {
  return( Teuchos::rcp( new CompoundOp<OP1,OP2>( op1, op2, temp ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDOP_HPP
