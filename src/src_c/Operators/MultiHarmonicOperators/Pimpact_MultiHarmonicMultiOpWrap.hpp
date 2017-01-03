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
/// wraps a \c MultiOperator and adds the functionality of handling \c MultiHarmonicField's which allows block methods
template<class MultiOperator>
class MultiHarmonicMultiOpWrap  {

public:

  using DomainFieldT = MultiHarmonicField<typename MultiOperator::DomainFieldT::InnerFieldT>;
  using RangeFieldT = MultiHarmonicField<typename MultiOperator::RangeFieldT::InnerFieldT>;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  using Ordinal = typename SpaceT::Ordinal;

  Teuchos::RCP<MultiOperator> op_;

public:

  MultiHarmonicMultiOpWrap( const Teuchos::RCP<MultiOperator>& op ): op_(op) {};

	void apply( const DomainFieldT& x, RangeFieldT& y, const Belos::ETrans&
			trans=Belos::NOTRANS) const {

		Teuchos::RCP<typename MultiOperator::DomainFieldT> mx =
			Teuchos::rcp( new typename MultiOperator::DomainFieldT( space(), 0 ) );

		Teuchos::RCP<typename MultiOperator::RangeFieldT>	 my = Teuchos::rcp( new
				typename MultiOperator::RangeFieldT ( space(), 0 ) );

		for( Ordinal i=std::max(space()->begin(F::U,3),1); i<=space()->end(F::U,3); ++i ) {
			// making x 
			mx->push_back(
					Teuchos::rcpFromRef(
						const_cast<typename MultiOperator::DomainFieldT::InnerFieldT&>(
							x.getConstCField(i) ) ) );
			mx->push_back(
					Teuchos::rcpFromRef(
						const_cast<typename MultiOperator::DomainFieldT::InnerFieldT&>(
							x.getConstSField(i) ) ) );

			// making y
			my->push_back( Teuchos::rcpFromRef(y.getCField(i)) );
			my->push_back( Teuchos::rcpFromRef(y.getSField(i)) );
		}

		if( 0==space()->begin(F::U,3) ) {
			mx->push_back(
					Teuchos::rcpFromRef(
						const_cast<typename MultiOperator::DomainFieldT::InnerFieldT&>(
							x.getConst0Field() ) ) );
			my->push_back( Teuchos::rcpFromRef( y.get0Field() ) );
		}

		// applying MultiField operator
		op_->apply( *mx, *my );

		y.changed();
	};

	void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}

  Teuchos::RCP<MultiOperator> getOperatorPtr() { return( op_ ); }

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
