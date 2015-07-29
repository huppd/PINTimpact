#pragma once
#ifndef PIMPACT_COMPOUNDSMOOTHER_HPP
#define PIMPACT_COMPOUNDSMOOTHER_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_MultiField.hpp"




namespace Pimpact {



/// \ingroup CompoundOperator
///
/// \tparam CompoundOpT should be of compatible to CompoundOpWrap
/// \f[ \begin{bmatrix} opV2V^{-1} & opS2V \\ 0 & -opS2S^{-1} \end{bmatrix}^{-1} \mathbf{x} = \mathbf{y} \f]
template<class CompoundOpT, temlate<class> class vSmoother, template<class,class,class> class sSmoother >
class CompoundSmoother {

public:

	typedef typename CompoundOpT::OpV2V OpV2VT;
	typedef typename CompoundOpT::OpS2V OpS2VT;
	typedef typename CompoundOpT::OpV2S OpV2ST;

  typedef typename CompoundOpT::DomainFieldT DomainFieldT;
  typedef typename CompoundOpT::RangeFieldT  RangeFieldT;

  typedef typename CompoundOpT::SpaceT SpaceT;

protected:

  typedef typename OpS2V::RangeFieldT  VF;
  typedef typename OpS2V::DomainFieldT  SF;

	typedef vSmoother< OpV2V> OpVSmoother;
	typedef sSmoother< OpV2V, OpV2S, OpS2V >  OpSSmoother;

  Teuchos::RCP<OpV2VT> opV2V_;
  Teuchos::RCP<OpS2VT> opS2V_;
  Teuchos::RCP<OpV2ST> opV2S_;

  Teuchos::RCP<OpVSmoother> opVSmoother_;
  Teuchos::RCP<OpSSmoother> opSSmoother_;

  Teuchos::RCP<VF> tempv_;
  Teuchos::RCP<SF> temps_;

public:

	CompoundSmoother( const Teuchos::RCP< CompoundOpT >& op ):
		opV2V_( op->getOpV2V() ),
		opS2V_( op->getOpS2V() ),
		opV2S_( op->getOpV2S() ),
		tempv_( create<VF>(opV2V_->space()) ),
		temps_( create<SF>(opV2V_->space()) )
	{};

	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		// ~ (D H^{-1} G)^{-1} p = D H^{-1} f_u - f_p
		temps_->add( -1., x.getConstSField(), 0., *temps_ );
		opS2S_->apply( *temps_ ,  y.getSField() );

		opS2V_->apply( y.getConstSField(), *tempv_ );

		tempv_->add( -1., *tempv_, 1., x.getConstVField() );

		opV2V_->apply( *tempv_, y.getVField() );

	}


	void assignField( const DomainFieldT& mv ) {
		opV2V_->assignField( mv.getConstVField());
    //    opS2V_->assignField( mv.getConstVField() );
    //    opV2S_->assignField( mv.getConstVField() );
	};

	Teuchos::RCP<const SpaceT> space() const { return(opV2V_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "CompoundSmoother " ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
    opV2V_->print( out );
		opS2V_->print( out );
		opS2S_->print( out );
  }

}; // end of class CompoundSmoother



/// \relates CompoundSmoother
template< class OpV2V, class OpS2V, class OpS2S >
Teuchos::RCP< CompoundSmoother<OpV2V,OpS2V,OpS2S> >
createCompoundSmoother(
    const Teuchos::RCP<OpV2V>& opV2V,
    const Teuchos::RCP<OpS2V>& opS2V,
    const Teuchos::RCP<OpS2S>& opS2S ) {

  //  return Teuchos::null;
  return(
      Teuchos::rcp( new CompoundSmoother<OpV2V,OpS2V,OpS2S>(opV2V,opS2V,opS2S) )
			);

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDSMOOTHER_HPP
