#pragma once
#ifndef PIMPACT_MULTIHARMONICMULTIOPWRAP_HPP
#define PIMPACT_MULTIHARMONICMULTIOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_MultiField.hpp"



namespace Pimpact {


/// \brief Operator wrapper.
/// \ingroup MultiHarmonicOperator
/// wraps a \ref MultiOp and adds the functionality of handling \c MultiHarmonicField's which allows block methods
template<class MultiOperator>
class MultiHarmonicMultiOpWrap  {

public:

  typedef MultiHarmonicField<typename MultiOperator::DomainFieldT::FieldT> DomainFieldT;
  typedef MultiHarmonicField<typename MultiOperator::RangeFieldT::FieldT> RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP<MultiOperator> op_;

public:

  MultiHarmonicMultiOpWrap( const Teuchos::RCP<MultiOperator>& op ): op_(op) {};

  void apply( const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS) const {

		auto mx = Teuchos::rcp( new typename MultiOperator::DomainFieldT( space(), (int)0 ) );
		auto my = Teuchos::rcp( new typename MultiOperator::RangeFieldT ( space(), (int)0 ) );

		mx->push_back(
				Teuchos::rcp_const_cast<typename MultiOperator::DomainFieldT::FieldT>(
					x.getConst0FieldPtr()
					)
				);
		my->push_back( y.get0FieldPtr() );

    int m = x.getNumberModes();

    for( int i=0; i<m; ++i ) {
			mx->push_back(
				Teuchos::rcp_const_cast<typename MultiOperator::DomainFieldT::FieldT>(
						x.getConstCFieldPtr(i)
						)
					);
			mx->push_back(
					Teuchos::rcp_const_cast<typename MultiOperator::DomainFieldT::FieldT>(
						x.getConstSFieldPtr(i)
						)
					);
			my->push_back( y.getCFieldPtr(i) );

			my->push_back( y.getSFieldPtr(i) );
    }
//		std::cout << "x numberVecs: " << mx->getNumberVecs() << "\n";
//		std::cout << "y numberVecs: " << my->getNumberVecs() << "\n";
		op_->apply( *mx, *my );
  };

  void assignField( const DomainFieldT& mv ) {
//		op_->assignField( Pimpact::createMultiField(mv.getConst0FieldPtr()) );
  };

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}

  Teuchos::RCP<MultiOperator> getOperatorPtr() { return( op_ ); }

	const std::string getLabel() const { return( "MultiHarmonicMultiOpWrap " ); };

}; // end of class MultiHarmonicMultiOpWrap


/// \relates MultiHarmonicMultiOpWrap
template<class MOT>
Teuchos::RCP< MultiHarmonicMultiOpWrap<MOT> >
createMultiHarmonicMultiOpWrap( const Teuchos::RCP<MOT>& op) {
  return( Teuchos::rcp( new MultiHarmonicMultiOpWrap<MOT>( op ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIHARMONICMULTIOPWRAP_HPP
