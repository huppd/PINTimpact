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

	typedef OpV2V OpV2VT;
	typedef OpS2V OpS2VT;
	typedef OpV2S OpV2ST;

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

	CompoundOpWrap( const Teuchos::RCP<const SpaceT>& space ):
		temp_(  create<VF   >(space) ),
		opV2V_( create<OpV2V>(space) ),
		opS2V_( create<OpS2V>(space) ),
		opV2S_( create<OpV2S>(space) ) {};

	CompoundOpWrap(
			const Teuchos::RCP<OpV2V>& opV2V,
			const Teuchos::RCP<OpS2V>& opS2V,
			const Teuchos::RCP<OpV2S>& opV2S ):
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

  };


	Teuchos::RCP<const SpaceT> space() const { return(opV2V_->space()); };


  Teuchos::RCP<OpV2V> getOpV2V() const { return( opV2V_ ); }
  Teuchos::RCP<OpS2V> getOpS2V() const { return( opS2V_ ); }
  Teuchos::RCP<OpV2S> getOpV2S() const { return( opV2S_ ); }

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "Compound( "+opV2V_->getLabel()+", "+opS2V_->getLabel()+", "+opV2S_->getLabel() +" )" ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		opV2V_->print( out );
		opS2V_->print( out );
		opV2S_->print( out );
  }

}; // end of class CompoundOpWrap



/// \relates CompoundOpWrap
template< class OpV2V, class OpS2V, class OpV2S >
Teuchos::RCP<CompoundOpWrap<OpV2V,OpS2V,OpV2S> > createCompoundOpWrap(
    const Teuchos::RCP<OpV2V>& opV2V,
    const Teuchos::RCP<OpS2V>& opS2V,
    const Teuchos::RCP<OpV2S>& opV2S ) {

  return(
      Teuchos::rcp( new CompoundOpWrap<OpV2V,OpS2V,OpV2S>(opV2V,opS2V,opV2S) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDSOPWRAP_HPP
