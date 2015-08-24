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
//  Teuchos::RCP<VF> tempv2_;

  Teuchos::RCP<OpV2V> opV2V_;
  Teuchos::RCP<OpS2V> opS2V_;
  Teuchos::RCP<OpS2S> opS2S_;

public:

	InverseTriangularOp(
			const Teuchos::RCP<OpV2V>& opV2V,
			const Teuchos::RCP<OpS2V>& opS2V,
			const Teuchos::RCP<OpS2S>& opS2S ):
		tempv_( create<VF>(opV2V->space()) ),
//		tempv2_( create<VF>(opV2V->space()) ),
		opV2V_(opV2V),
		opS2V_(opS2V),
		opS2S_(opS2S) {};

	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		opS2S_->apply( x.getConstSField() ,  y.getSField() );
		y.getSField().scale( -1. );

		opS2V_->apply( y.getConstSField(), *tempv_ );

//		tempv2_->add( -1., *tempv_, 1., x.getConstVField() );
		tempv_->add( -1., *tempv_, 1., x.getConstVField() );

		opV2V_->apply( *tempv_, y.getVField() );
//		opV2V_->apply( x.getConstVField(), y.getVField() );

	}


	void assignField( const DomainFieldT& mv ) {
		opV2V_->assignField( mv.getConstVField());
    //    opS2V_->assignField( mv.getConstVField() );
    //    opV2S_->assignField( mv.getConstVField() );
	};

	Teuchos::RCP<const SpaceT> space() const { return(opV2V_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "InverseTriangularOp " ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
    opV2V_->print( out );
		opS2V_->print( out );
		opS2S_->print( out );
  }

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
