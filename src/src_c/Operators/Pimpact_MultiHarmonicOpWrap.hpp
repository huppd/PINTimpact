#pragma once
#ifndef PIMPACT_MULTIHARMONICOPWRAP_HPP
#define PIMPACT_MULTIHARMONICOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_MultiField.hpp"



namespace Pimpact {


/// \brief Operator wrapper.
/// \ingroup MultiHarmonicOperator
/// wraps a \ref BaseOperator "base operator" and adds the functionality of handling \c MultiHarmonicField's.
template<class Operator>
class MultiHarmonicOpWrap  {

  Teuchos::RCP<Operator> op_;

public:

  typedef MultiHarmonicField<typename Operator::DomainFieldT> DomainFieldT;
  typedef MultiHarmonicField<typename Operator::RangeFieldT> RangeFieldT;
  typedef typename Operator::OpType OpType;

  MultiHarmonicOpWrap():op_( Teuchos::rcp( new Operator() ) ) {};
  MultiHarmonicOpWrap( const Teuchos::RCP<Operator>& op ):op_(op) {};
  ~MultiHarmonicOpWrap() {op_=Teuchos::null;};


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
    op_->assignField( mv.getConstField(0) );
  };

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

  Teuchos::RCP<Operator> getOperatorPtr() { return( op_ ); }


}; // end of class MultiHarmonicOpWrap



/// \relates MultiHarmonicOpWrap
template<class Operator>
Teuchos::RCP< MultiHarmonicOpWrap<Operator> > createMultiHarmonicOpWrap( const Teuchos::RCP<Operator>& op=Teuchos::null) {
  if( Teuchos::is_null(op) )
    return( Teuchos::rcp( new MultiHarmonicOpWrap<Operator>( Teuchos::rcp( new Operator() ) ) ) );
  else
    return( Teuchos::rcp( new MultiHarmonicOpWrap<Operator>( op ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIHARMONICOPWRAP_HPP
