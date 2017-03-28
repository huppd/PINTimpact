#pragma once
#ifndef PIMPACT_MULTIHARMONICMULTIOPWRAP_HPP
#define PIMPACT_MULTIHARMONICMULTIOPWRAP_HPP


#include <Teuchos_RCP.hpp>

#include <BelosTypes.hpp>

#include "Pimpact_MultiField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief Operator wrapper.
/// \ingroup MultiHarmonicOperator
/// wraps a \c OpT and adds the functionality of handling \c MultiHarmonicField's which allows block methods
/// \notes coupled to InverseOp
template<class OpT>
class MultiHarmonicMultiOpWrap  {

public:

  using DomainFieldT = MultiHarmonicField<typename OpT::DomainFieldT>;
  using RangeFieldT = MultiHarmonicField<typename OpT::RangeFieldT>;

  using SpaceT = typename OpT::SpaceT;

protected:

  using Ordinal = typename SpaceT::Ordinal;

  using MDomainFT = MultiField<typename DomainFieldT::InnerFieldT>;
  using MRanageFT = MultiField<typename RangeFieldT::InnerFieldT>;

  Teuchos::RCP<OpT> op_;

public:

  MultiHarmonicMultiOpWrap( const Teuchos::RCP<OpT>& op ): op_(op) {};


	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		MDomainFT mx( space(), 0 );
		MRanageFT my( space(), 0 );

		for( Ordinal i=std::max(space()->si(F::U,3),1); i<=space()->ei(F::U,3); ++i ) {
			// making x 
			mx.push_back(
					Teuchos::rcpFromRef(
						const_cast<typename DomainFieldT::InnerFieldT&>(
							x.getCField(i) ) ) );
			mx.push_back(
					Teuchos::rcpFromRef(
						const_cast<typename DomainFieldT::InnerFieldT&>(
							x.getSField(i) ) ) );

			// making y
			my.push_back( Teuchos::rcpFromRef(y.getCField(i)) );
			my.push_back( Teuchos::rcpFromRef(y.getSField(i)) );
		}

		if( 0==space()->si(F::U,3) ) {
			mx.push_back(
					Teuchos::rcpFromRef(
						const_cast<typename DomainFieldT::InnerFieldT&>(
							x.get0Field() ) ) );
			my.push_back( Teuchos::rcpFromRef( y.get0Field() ) );
		}

		// applying MultiField operator
		op_->apply( mx, my );

		y.changed();
	};

	void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}

  Teuchos::RCP<OpT> getOperatorPtr() { return( op_ ); }

	const std::string getLabel() const { return( "MH(M)_"+op_->getLabel() ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		op_->print( out );
  }

}; // end of class MultiHarmonicMultiOpWrap


/// \relates MultiHarmonicMultiOpWrap
template<class MOT>
Teuchos::RCP< MultiHarmonicMultiOpWrap<MOT> >
createMultiHarmonicMultiOpWrap( const Teuchos::RCP<MOT>& op) {
  return( Teuchos::rcp( new MultiHarmonicMultiOpWrap<MOT>( op ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIHARMONICMULTIOPWRAP_HPP
