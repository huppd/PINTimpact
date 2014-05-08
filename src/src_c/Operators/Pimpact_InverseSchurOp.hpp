#pragma once
#ifndef PIMPACT_INVERSESCHUROP_HPP
#define PIMPACT_INVERSESCHUROP_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"

#include "Pimpact_TripleComposition.hpp"
#include "Pimpact_LinearProblem.hpp"



namespace Pimpact {


/// \ingroup CompoundOperator
template<class OpV2Vinv,class OpS2V, class OpV2S>
class InverseSchurOp {

protected:
  typedef typename OpV2Vinv::DomainFieldT  VF;
  typedef typename OpS2V::DomainFieldT  SF;

  Teuchos::RCP<VF> tempv_;
  Teuchos::RCP<SF> temps_;

  Teuchos::RCP<OpV2Vinv> opV2Vinv_;
  Teuchos::RCP<OpS2V> opS2V_;
  Teuchos::RCP<OpV2S> opV2S_;

  Teuchos::RCP< LinearProblem< MultiField<SF> > > lp_;

public:

  typedef CompoundField<VF,SF>  DomainFieldT;
  typedef CompoundField<VF,SF>  RangeFieldT;

  InverseSchurOp(
      const Teuchos::RCP<VF> temp=Teuchos::null,
      const Teuchos::RCP<OpV2Vinv>& opV2V=Teuchos::null,
      const Teuchos::RCP<OpS2V>& opS2V=Teuchos::null,
      const Teuchos::RCP<OpV2S>& opV2S=Teuchos::null ):
        tempv_(temp),
        opV2Vinv_(opV2V),
        opS2V_(opS2V),
        opV2S_(opV2S) {

    if( opV2Vinv_.is_null() ) opV2Vinv_ = Teuchos::rcp( new OpV2Vinv() );
    if( opV2S_.is_null() ) opV2S_ = Teuchos::rcp( new OpV2S() );
    if( opS2V_.is_null() ) opS2V_ = Teuchos::rcp( new OpS2V() );

    auto opSchur =
        Pimpact::createMultiOperatorBase(
            Pimpact::createTripleCompositionOp(
                tempv_->getConstFieldPtr(0)->clone(),
                tempv_->getConstFieldPtr(0)->clone(),
                opS2V_,
                opV2Vinv_,
                opV2S_ ) );

    lp_ = Pimpact::createLinearProblem( opSchur, Teuchos::null, Teuchos::null,Teuchos::parameterList(),"GMRES" );

  };

  void apply(const DomainFieldT& x, RangeFieldT& y ) const {
    // ~ H^{-1} f_u
    opV2Vinv_->apply( x.getConstVField(), *tempv_ );
    // ~ D H^{-1} f_u
    opV2S_->apply( *tempv, *tempp_ );
    // f_p
    // ~ (D H^{-1} G)^{-1} p = f_p

    // ~grad
    opS2V_->apply( x.getConstSField(), *tempv_ );
    y.getVField().add( 1., y.getConstVField(), 1., *tempv_ );
    opV2S_->apply( x.getConstVField(), y.getSField() );
  }

  void assignField( const DomainFieldT& mv ) {
    opV2Vinv_->assignField( mv.getConstVField() );
//    opS2V_->assignField( mv.getConstVField() );
//    opV2S_->assignField( mv.getConstVField() );
  };

  bool hasApplyTranspose() const { return( false ); }

}; // end of class InverseSchurOp



/// \relates InverseSchurOp
template< class OpV2V, class OpS2V, class OpV2S >
Teuchos::RCP<InverseSchurOp<OpV2V,OpS2V,OpV2S> > createInverseSchurOp(
    const Teuchos::RCP<typename OpV2V::DomainFieldT>& temp,
    const Teuchos::RCP<OpV2V>& opV2V=Teuchos::null,
    const Teuchos::RCP<OpS2V>& opS2V=Teuchos::null,
    const Teuchos::RCP<OpV2S>& opV2S=Teuchos::null ) {

  return(
      Teuchos::rcp( new InverseSchurOp<OpV2V,OpS2V,OpV2S>(temp,opV2V,opS2V,opV2S) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSESCHUROP_HPP
