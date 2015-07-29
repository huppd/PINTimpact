#pragma once
#ifndef PIMPACT_MULTIHARMONICMODEOPWRAP_HPP
#define PIMPACT_MULTIHARMONICMODEOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_MultiField.hpp"



namespace Pimpact {


/// \brief Operator wrapper.
/// \ingroup MultiHarmonicOperator
/// wraps a \ref BaseOperator "base operator" and adds the functionality of handling \c MultiHarmonicField's.
template<class OpT>
class MultiHarmonicOpWrap  {

public:

  typedef MultiHarmonicField<typename OpT::DomainFieldT> DomainFieldT;
  typedef MultiHarmonicField<typename OpT::RangeFieldT> RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP<OpT> op_;

public:

  MultiHarmonicOpWrap( const Teuchos::RCP<const SpaceT>& space ):
		op_( create<OpT>(space) ) {};

  MultiHarmonicOpWrap( const Teuchos::RCP<OpT>& op ): op_(op) {};

  void apply( const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS) const {

    op_->apply( x.getConst0Field(), y.get0Field() );

    int m = x.getNumberModes();

    for( int i=0; i<m; ++i ) {
      op_->apply( x.getConstCField(i), y.getCField(i) );
      op_->apply( x.getConstSField(i), y.getSField(i) );
    }
  };

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( mv.getConst0Field() );
  };

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}

  Teuchos::RCP<OpT> getOperatorPtr() { return( op_ ); }

	const std::string getLabel() const { return( "MH_"+op_->getLabel() ); };

  void print( std::ostream& out=std::cout ) const {
		out <<  getLabel() << ":\n";
    op_->print( out );
  }

}; // end of class MultiHarmonicOpWrap



/// \relates MultiHarmonicOpWrap
template<class OpT>
Teuchos::RCP< MultiHarmonicOpWrap<OpT> >
createMultiHarmonicOpWrap( const Teuchos::RCP<OpT>& op) {
  return( Teuchos::rcp( new MultiHarmonicOpWrap<OpT>( op ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIHARMONICMODEOPWRAP_HPP
