#pragma once
#ifndef PIMPACT_INVERSETRIANGULAROP_HPP
#define PIMPACT_INVERSETRIANGULAROP_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_MultiField.hpp"

//#include "Pimpact_TripleCompositionOp.hpp"
//#include "Pimpact_LinearProblem.hpp"

//#include "Pimpact_OperatorFactory.hpp"



namespace Pimpact {


/// \ingroup CompoundOperator
template<class OpV2Vinv,class OpS2V, class OpS2Sinv>
class InverseTriangularOp {

protected:

  typedef typename OpS2V::RangeFieldT  VF;
  typedef typename OpS2V::DomainFieldT  SF;

public:

  typedef CompoundField<VF,SF>  DomainFieldT;
  typedef CompoundField<VF,SF>  RangeFieldT;

  Teuchos::RCP<VF> tempv_;
  Teuchos::RCP<SF> temps_;

  Teuchos::RCP<OpV2Vinv> opV2Vinv_;
  Teuchos::RCP<OpS2V>    opS2V_;
  Teuchos::RCP<OpS2Sinv> opS2Sinv_;



  InverseTriangularOp(
      const Teuchos::RCP<VF> tempv=Teuchos::null,
      const Teuchos::RCP<SF> temps=Teuchos::null,
      const Teuchos::RCP<OpV2Vinv>& opV2V=Teuchos::null,
      const Teuchos::RCP<OpS2V>& opS2V=Teuchos::null,
      const Teuchos::RCP<OpS2Sinv>& opS2Sinv=Teuchos::null ):
        tempv_(tempv),
        temps_(temps),
  		opV2Vinv_(opV2V),
        opS2V_(opS2V),
        opS2Sinv_(opS2Sinv) {

    if( opV2Vinv_.is_null() ) opV2Vinv_ = Teuchos::rcp( new OpV2Vinv() );
    if( opS2V_.is_null() )    opS2V_ = Teuchos::rcp( new OpS2V() );
    if( opS2Sinv_.is_null() )  opS2Sinv_ = Teuchos::rcp( new OpS2Sinv() );


  };

  void apply(const DomainFieldT& x, RangeFieldT& y ) const {

    //    ----Triangular Schur
    temps_->add( -1., x.getConstSField(), 0., *temps_ );
    opS2Sinv_->apply( *createMultiField(temps_) , *createMultiField( y.getSFieldPtr()));

    opS2V_->apply( y.getConstSField(), *tempv_ );

    tempv_->add( -1., *tempv_, 1., x.getConstVField() );

    opV2Vinv_->apply( *createMultiField(tempv_), *createMultiField( y.getVFieldPtr() ) );
//        y.getSFieldPtr()->add( 0., y.getConstSField(), 1., x.getConstSField() );
    // ~ (D H^{-1} G)^{-1} p = D H^{-1} f_u - f_p

  }

  /// \todo fixme
  void assignField( const DomainFieldT& mv ) {
    opV2Vinv_->assignField( *createMultiField( Teuchos::rcp_const_cast<VF>(mv.getConstVFieldPtr()) ) );
//    opS2V_->assignField( mv.getConstVField() );
//    opV2Sinv_->assignField( mv.getConstVField() );
  };

  bool hasApplyTranspose() const { return( false ); }

}; // end of class InverseTriangularOp



/// \relates InverseTriangularOp
template< class OpV2V, class OpS2V, class OpS2S >
Teuchos::RCP< InverseTriangularOp<OpV2V,OpS2V,OpS2S> >
createInverseTriangularOp(
    const Teuchos::RCP<typename OpS2V::RangeFieldT>& tempv,
    const Teuchos::RCP<typename OpS2V::DomainFieldT>& temps,
    const Teuchos::RCP<OpV2V>& opV2V,
    const Teuchos::RCP<OpS2V>& opS2V,
    const Teuchos::RCP<OpS2S>& opS2Sinv=Teuchos::null ) {

//  return Teuchos::null;
  return(
      Teuchos::rcp( new InverseTriangularOp<OpV2V,OpS2V,OpS2S>(tempv,temps,opV2V,opS2V,opS2Sinv) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSETRIANGULAROP_HPP
