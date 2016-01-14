#pragma once
#ifndef PIMPACT_COMPOUNDSMOOTHER_HPP
#define PIMPACT_COMPOUNDSMOOTHER_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_TripleCompositionOp.hpp"




namespace Pimpact {



/// \ingroup CompoundOperator
///
/// \tparam CompoundOpT should be of compatible to CompoundOpWrap
/// \f[ \begin{bmatrix} opV2V^{-1} & opS2V \\ 0 & -opS2S^{-1} \end{bmatrix}^{-1} \mathbf{x} = \mathbf{y} \f]
template<class CompoundOpT, template<class> class vSmoother, class sSmoother>
class CompoundSmoother {

public:

	typedef typename CompoundOpT::OpV2VT OpV2VT;
	typedef typename CompoundOpT::OpS2VT OpS2VT;
	typedef typename CompoundOpT::OpV2ST OpV2ST;

  typedef typename CompoundOpT::DomainFieldT DomainFieldT;
  typedef typename CompoundOpT::RangeFieldT  RangeFieldT;

  typedef typename CompoundOpT::SpaceT SpaceT;

protected:

  typedef typename OpS2VT::RangeFieldT  VF;
  typedef typename OpS2VT::DomainFieldT  SF;

	typedef vSmoother< OpV2VT> OpVSmoother;
	typedef TripleCompositionOp< sSmoother,TripleCompositionOp<OpV2ST,OpV2VT,OpS2VT>,sSmoother> OpSSmootherT;

  Teuchos::RCP<OpS2VT> opS2V_;

  Teuchos::RCP<OpVSmoother> opVSmoother_;
  Teuchos::RCP<OpSSmootherT> opSSmoother_;

  Teuchos::RCP<VF> tempv_;
//  Teuchos::RCP<SF> temps_;

public:

	CompoundSmoother(
			const Teuchos::RCP< CompoundOpT >& op,
		  Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::null	):
		opS2V_( op->getOpS2V() ),
		opVSmoother_( Teuchos::rcp( new OpVSmoother( op->getOpV2V(), Teuchos::rcpFromRef(pl->sublist("VSmoother") ) ) ) ),
		tempv_( create<VF>( op->space() ) ) {

		auto bla = createTripleCompositionOp( op->getOpV2S(), op->getOpV2V(), op->getOpS2V() );
		auto blup = Teuchos::rcp( new sSmoother( space() ) );

		opSSmoother_ = createTripleCompositionOp( blup, bla, blup );

	};

	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		opSSmoother_->apply( x.getConstSField() ,  y.getSField() );
		y.getSField().scale( -1. );

		opS2V_->apply( y.getConstSField(), *tempv_ );

		tempv_->add( -1., *tempv_, 1., x.getConstVField() );

		opVSmoother_->apply( *tempv_, y.getVField() );

	}


	void assignField( const DomainFieldT& mv ) {
	};

	Teuchos::RCP<const SpaceT> space() const { return(opS2V_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "CompoundSmoother" ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
//    opV2V_->print( out );
//		opS2V_->print( out );
//		opS2S_->print( out );
  }

}; // end of class CompoundSmoother



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDSMOOTHER_HPP
