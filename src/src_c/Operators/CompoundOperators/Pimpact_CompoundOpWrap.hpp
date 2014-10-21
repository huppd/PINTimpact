#pragma once
#ifndef PIMPACT_COMPOUNDSOPWRAP_HPP
#define PIMPACT_COMPOUNDSOPWRAP_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"


namespace Pimpact {


/// \ingroup CompoundOperator
template<class OpV2V,class OpS2V, class OpV2S>
class CompoundOpWrap {

protected:

  typedef typename OpV2V::DomainFieldT  VF;
  typedef typename OpS2V::DomainFieldT  SF;

  Teuchos::RCP<VF> temp_;

  Teuchos::RCP<OpV2V> opV2V_;
  Teuchos::RCP<OpS2V> opS2V_;
  Teuchos::RCP<OpV2S> opV2S_;

public:

  typedef CompoundField<VF,SF>  DomainFieldT;
  typedef CompoundField<VF,SF>  RangeFieldT;

  CompoundOpWrap(
      const Teuchos::RCP<VF> temp/*=Teuchos::null*/,
      const Teuchos::RCP<OpV2V>& opV2V/*=Teuchos::null*/,
      const Teuchos::RCP<OpS2V>& opS2V/*=Teuchos::null*/,
      const Teuchos::RCP<OpV2S>& opV2S/*=Teuchos::null*/
      ):
        temp_(temp),
        opV2V_(opV2V),
        opS2V_(opS2V),
        opV2S_(opV2S) {
//    if( opV2V_.is_null() ) opV2V_ = Teuchos::rcp( new OpV2V() );
//    if( opV2S_.is_null() ) opV2S_ = Teuchos::rcp( new OpV2S() );
//    if( opS2V_.is_null() ) opS2V_ = Teuchos::rcp( new OpS2V() );
  };

  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS  ) const {
    // H-blockz
    opV2V_->apply( x.getConstVField(), y.getVField() );
    // ~grad
    opS2V_->apply( x.getConstSField(), *temp_ );
    y.getVField().add( 1., y.getConstVField(), 1., *temp_ );
    // ~div
    opV2S_->apply( x.getConstVField(), y.getSField() );
  }

  void assignField( const DomainFieldT& mv ) {
    opV2V_->assignField( mv.getConstVField() );
//    opS2V_->assignField( mv.getConstVField() );
//    opV2S_->assignField( mv.getConstVField() );
  };

  bool hasApplyTranspose() const { return( false ); }

}; // end of class CompoundOpWrap



/// \relates CompoundOpWrap
template< class OpV2V, class OpS2V, class OpV2S >
Teuchos::RCP<CompoundOpWrap<OpV2V,OpS2V,OpV2S> > createCompoundOpWrap(
    const Teuchos::RCP<typename OpV2V::DomainFieldT>& temp,
    const Teuchos::RCP<OpV2V>& opV2V=Teuchos::null,
    const Teuchos::RCP<OpS2V>& opS2V=Teuchos::null,
    const Teuchos::RCP<OpV2S>& opV2S=Teuchos::null ) {

  return(
      Teuchos::rcp( new CompoundOpWrap<OpV2V,OpS2V,OpV2S>(temp,opV2V,opS2V,opV2S) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDSOPWRAP_HPP
