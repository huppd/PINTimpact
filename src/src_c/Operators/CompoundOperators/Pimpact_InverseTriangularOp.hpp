#pragma once
#ifndef PIMPACT_INVERSETRIANGULAROP_HPP
#define PIMPACT_INVERSETRIANGULAROP_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_MultiField.hpp"




namespace Pimpact {



/// \ingroup CompoundOperator
///
/// \f[ \begin{bmatrix} opV2V^{-1} & opS2V \\ 0 & -opS2S^{-1} \end{bmatrix}^{-1} \mathbf{x} = \mathbf{y} \f]
template<class OpV2V,class OpS2V, class OpS2S>
class InverseTriangularOp {

  typedef typename OpS2V::RangeFieldT  VF;
  typedef typename OpS2V::DomainFieldT  SF;

public:

  typedef CompoundField<VF,SF>  DomainFieldT;
  typedef CompoundField<VF,SF>  RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP<VF> tempv_;
  Teuchos::RCP<SF> temps_;

  Teuchos::RCP<OpV2V> opV2V_;
  Teuchos::RCP<OpS2V> opS2V_;
  Teuchos::RCP<OpS2S> opS2S_;

public:

  /// \todo update
  InverseTriangularOp(
      const Teuchos::RCP<VF> tempv=Teuchos::null,
      const Teuchos::RCP<SF> temps=Teuchos::null,
      const Teuchos::RCP<OpV2V>& opV2V=Teuchos::null,
      const Teuchos::RCP<OpS2V>& opS2V=Teuchos::null,
      const Teuchos::RCP<OpS2S>& opS2S=Teuchos::null ):
        tempv_(tempv),
        temps_(temps),
        opV2V_(opV2V),
        opS2V_(opS2V),
        opS2S_(opS2S) {

    if( opV2V_.is_null() ) opV2V_ = Teuchos::rcp( new OpV2V() );
    if( opS2V_.is_null() )    opS2V_ = Teuchos::rcp( new OpS2V() );
    if( opS2S_.is_null() )  opS2S_ = Teuchos::rcp( new OpS2S() );


  };

  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    //    ----Triangular Schur
    temps_->add( -1., x.getConstSField(), 0., *temps_ );
    opS2S_->apply( *createMultiField(temps_) , *createMultiField( y.getSFieldPtr()));

    opS2V_->apply( y.getConstSField(), *tempv_ );

    tempv_->add( -1., *tempv_, 1., x.getConstVField() );

    opV2V_->apply( *createMultiField(tempv_), *createMultiField( y.getVFieldPtr() ) );
    //        y.getSFieldPtr()->add( 0., y.getConstSField(), 1., x.getConstSField() );
    // ~ (D H^{-1} G)^{-1} p = D H^{-1} f_u - f_p

  }

  /// \todo fixme
  void assignField( const DomainFieldT& mv ) {
    opV2V_->assignField( *createMultiField( Teuchos::rcp_const_cast<VF>(mv.getConstVFieldPtr()) ) );
    //    opS2V_->assignField( mv.getConstVField() );
    //    opV2S_->assignField( mv.getConstVField() );
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
    const Teuchos::RCP<OpS2S>& opS2S=Teuchos::null ) {

  //  return Teuchos::null;
  return(
      Teuchos::rcp( new InverseTriangularOp<OpV2V,OpS2V,OpS2S>(tempv,temps,opV2V,opS2V,opS2S) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSETRIANGULAROP_HPP
