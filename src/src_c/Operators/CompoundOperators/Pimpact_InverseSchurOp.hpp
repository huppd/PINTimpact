#pragma once
#ifndef PIMPACT_INVERSESCHUROP_HPP
#define PIMPACT_INVERSESCHUROP_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_TripleCompositionOp.hpp"
#include "Pimpact_LinearProblem.hpp"

#include "Pimpact_OperatorFactory.hpp"



namespace Pimpact {


/// \ingroup CompoundOperator
/// \relates TripleComposition
template<class OpV2Vinv,class OpS2V, class OpV2S>
class InverseSchurOp {

protected:

  typedef typename OpS2V::RangeFieldT  VF;
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
      const Teuchos::RCP<VF> tempv,
      const Teuchos::RCP<SF> temps,
      const Teuchos::RCP<OpV2Vinv>& opV2V,
      const Teuchos::RCP<OpS2V>& opS2V,
      const Teuchos::RCP<OpV2S>& opV2S,
      const Teuchos::RCP<Teuchos::ParameterList>& para=Teuchos::null):
        tempv_(tempv),
        temps_(temps),
        opV2Vinv_(opV2V),
        opS2V_(opS2V),
        opV2S_(opV2S) {



    auto opSchur =
        Pimpact::createOperatorBase< MultiField<SF> >(
            Pimpact::createTripleCompositionOp(
                createMultiField( tempv_->clone(Pimpact::ShallowCopy) ),
                createMultiField( tempv_->clone(Pimpact::ShallowCopy) ),
                createMultiOpWrap( opV2S_ ),
                opV2Vinv_,
                createMultiOpWrap( opS2V_ )
            )
        );

    lp_ = Pimpact::createLinearProblem< MultiField<SF> >( opSchur, Teuchos::null, Teuchos::null, para, "GMRES" );
//    lp_ = Pimpact::createLinearProblem< MultiField<SF> >( opSchur, Teuchos::null, Teuchos::null, Teuchos::null, "GmresPoly" );
//    lp_ = Pimpact::createLinearProblem< MultiField<SF> >( opSchur, Teuchos::null, Teuchos::null, para, "GCRODR" );
//    lp_ = Pimpact::createLinearProblem< MultiField<SF> >( opSchur, Teuchos::null, Teuchos::null, Teuchos::parameterList(), "GCRODR" );

  };

  void apply(const DomainFieldT& x, RangeFieldT& y ) const {
    //----full Schur complement
    // ~ H^{-1} f_u
    opV2Vinv_->apply( *createMultiField( Teuchos::rcp_const_cast<VF>(x.getConstVFieldPtr()) ), *createMultiField(tempv_) ); // should be correct
    // ~ D H^{-1} f_u
    opV2S_->apply( *tempv_, *temps_ );
    // ~ D H^{-1} f_u - f_p
    temps_->add( -1., x.getConstSField(), 1., *temps_ );
    // ~ (D H^{-1} G)^{-1} p = D H^{-1} f_u - f_p
    lp_->solve( createMultiField( y.getSFieldPtr()), createMultiField(temps_) );
    // ~ G p
    opS2V_->apply( y.getConstSField(), *tempv_ );
    // ~ f_u - G p
    tempv_->add( -1., *tempv_, 1., x.getConstVField() );
    // ~ H^{-1}(f_u -G p)
    opV2Vinv_->apply( *createMultiField(tempv_), *createMultiField( y.getVFieldPtr() ) );

//    //----Diag Schur
//    opV2Vinv_->apply( *createMultiField( Teuchos::rcp_const_cast<VF>(x.getConstVFieldPtr()) ), *createMultiField(y.getVFieldPtr()) ); // should be correct
////    y.getSFieldPtr()->add( 0., y.getConstSField(), 1., x.getConstSField() );
//    // ~ (D H^{-1} G)^{-1} p = D H^{-1} f_u - f_p
//    lp_->solve(
//        createMultiField( y.getSFieldPtr()),
//        createMultiField( Teuchos::rcp_const_cast<SF>(x.getConstSFieldPtr()) ) );
////    x.getSFieldPtr()->scale(-1);

//    ----Triangular Schur
//    temps_->add( -1., x.getConstSField(), 0., *temps_ );
//    lp_->solve( createMultiField( y.getSFieldPtr()), createMultiField(temps_) );
//    opS2V_->apply( y.getConstSField(), *tempv_ );
//    tempv_->add( -1., *tempv_, 1., x.getConstVField() );
//
//    opV2Vinv_->apply( *createMultiField(tempv_), *createMultiField( y.getVFieldPtr() ) );
////    y.getSFieldPtr()->add( 0., y.getConstSField(), 1., x.getConstSField() );
//    // ~ (D H^{-1} G)^{-1} p = D H^{-1} f_u - f_p

  }

  /// \todo fixme
  void assignField( const DomainFieldT& mv ) {
    opV2Vinv_->assignField( *createMultiField( Teuchos::rcp_const_cast<VF>(mv.getConstVFieldPtr()) ) );
//    opS2V_->assignField( mv.getConstVField() );
//    opV2S_->assignField( mv.getConstVField() );
  };

  bool hasApplyTranspose() const { return( false ); }

}; // end of class InverseSchurOp



/// \relates InverseSchurOp
template< class OpV2V, class OpS2V, class OpV2S >
Teuchos::RCP< InverseSchurOp<OpV2V,OpS2V,OpV2S> > createInverseSchurOp(
    const Teuchos::RCP<typename OpS2V::RangeFieldT>& tempv,
    const Teuchos::RCP<typename OpS2V::DomainFieldT>& temps,
    const Teuchos::RCP<OpV2V>& opV2V=Teuchos::null,
    const Teuchos::RCP<OpS2V>& opS2V=Teuchos::null,
    const Teuchos::RCP<OpV2S>& opV2S=Teuchos::null,
    const Teuchos::RCP<Teuchos::ParameterList>& para=Teuchos::null
    ) {

  return(
      Teuchos::rcp( new InverseSchurOp<OpV2V,OpS2V,OpV2S>(tempv,temps,opV2V,opS2V,opV2S,para) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSESCHUROP_HPP