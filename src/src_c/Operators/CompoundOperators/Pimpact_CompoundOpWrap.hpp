#pragma once
#ifndef PIMPACT_COMPOUNDSOPWRAP_HPP
#define PIMPACT_COMPOUNDSOPWRAP_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"


namespace Pimpact {


/// \ingroup CompoundOperator
///
/// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class OpV2V,class OpS2V, class OpV2S>
class CompoundOpWrap {

public:

  typedef typename OpV2V::DomainFieldT  VF;
  typedef typename OpS2V::DomainFieldT  SF;

  typedef CompoundField<VF,SF>  DomainFieldT;
  typedef CompoundField<VF,SF>  RangeFieldT;

  typedef typename VF::SpaceT SpaceT;

protected:

  Teuchos::RCP<VF> temp_;

  Teuchos::RCP<OpV2V> opV2V_;
  Teuchos::RCP<OpS2V> opS2V_;
  Teuchos::RCP<OpV2S> opV2S_;

public:

  /// \todo constructor from space
  CompoundOpWrap(
      const Teuchos::RCP<OpV2V>& opV2V,
      const Teuchos::RCP<OpS2V>& opS2V,
      const Teuchos::RCP<OpV2S>& opV2S
      ):
        temp_( create<VF>(opV2V->space()) ),
        opV2V_(opV2V),
        opS2V_(opS2V),
        opV2S_(opV2S) {};

  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS  ) const {
    // H-blockz
    opV2V_->apply( x.getConstVField(), y.getVField() );
    // ~grad
    opS2V_->apply( x.getConstSField(), *temp_ );
    y.getVField().add( 1., y.getConstVField(), 1., *temp_ );
    // ~div
    opV2S_->apply( x.getConstVField(), y.getSField() );
//		y.getSField().level();
  }

  void assignField( const DomainFieldT& mv ) {
    opV2V_->assignField( mv.getConstVField() );
//    opS2V_->assignField( mv.getConstVField() );
//    opV2S_->assignField( mv.getConstVField() );
  };

	Teuchos::RCP<const SpaceT> space() const { return(opV2V_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "CompoundOpWrap( "+opV2V_->getLabel()+", "+opS2V_->getLabel()+", "+opV2S_->getLabel() +" ) " ); };

}; // end of class CompoundOpWrap



/// \relates CompoundOpWrap
template< class OpV2V, class OpS2V, class OpV2S >
Teuchos::RCP<CompoundOpWrap<OpV2V,OpS2V,OpV2S> > createCompoundOpWrap(
    const Teuchos::RCP<OpV2V>& opV2V=Teuchos::null,
    const Teuchos::RCP<OpS2V>& opS2V=Teuchos::null,
    const Teuchos::RCP<OpV2S>& opV2S=Teuchos::null ) {

  return(
      Teuchos::rcp( new CompoundOpWrap<OpV2V,OpS2V,OpV2S>(opV2V,opS2V,opV2S) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDSOPWRAP_HPP
