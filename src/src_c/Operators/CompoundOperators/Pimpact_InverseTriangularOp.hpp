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

  InverseTriangularOp(
      const Teuchos::RCP<OpV2V>& opV2V,
      const Teuchos::RCP<OpS2V>& opS2V,
      const Teuchos::RCP<OpS2S>& opS2S ):
        tempv_( create<VF>(opV2V->space()) ),
        temps_( create<SF>(opV2V->space()) ),
        opV2V_(opV2V),
        opS2V_(opS2V),
        opS2S_(opS2S) {};

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

	Teuchos::RCP<const SpaceT> space() const { return(opV2V_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

}; // end of class InverseTriangularOp



/// \relates InverseTriangularOp
template< class OpV2V, class OpS2V, class OpS2S >
Teuchos::RCP< InverseTriangularOp<OpV2V,OpS2V,OpS2S> >
createInverseTriangularOp(
    const Teuchos::RCP<OpV2V>& opV2V,
    const Teuchos::RCP<OpS2V>& opS2V,
    const Teuchos::RCP<OpS2S>& opS2S ) {

  //  return Teuchos::null;
  return(
      Teuchos::rcp( new InverseTriangularOp<OpV2V,OpS2V,OpS2S>(opV2V,opS2V,opS2S) )
			);

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSETRIANGULAROP_HPP
