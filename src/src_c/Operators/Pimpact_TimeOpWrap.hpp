#pragma once
#ifndef PIMPACT_TIMEOPWRAP_HPP
#define PIMPACT_TIMEOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_TimeField.hpp"



namespace Pimpact {


/// \brief Operator wrapper for \c TimeField's.
/// \ingroup TimeOperator
/// wraps and \c Operator and adds the functionality of handling \c TimeField.
template<class Operator>
class TimeOpWrap  {

  Teuchos::RCP<Operator> op_;

public:

  typedef TimeField<typename Operator::DomainFieldT> DomainFieldT;
  typedef TimeField<typename Operator::RangeFieldT> RangeFieldT;

  TimeOpWrap():op_( Teuchos::rcp( new Operator() ) ) {};
  TimeOpWrap( const Teuchos::RCP<Operator>& op ):op_(op) {};
  ~TimeOpWrap() { op_=Teuchos::null; };


  /// \brief default apply
  void apply( const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS) const {

//    const_cast<DomainFieldT>(x).exchange();

    typename DomainFieldT::Iter j = y.beginI_;
    for( typename DomainFieldT::Iter i=x.beginI_; i<x.endI_; ++i )
      op_->apply( **i , **(j++) );
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( mv.getConstField(0) );
  };

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

  Teuchos::RCP<Operator> getOperatorPtr() { return( op_ ); }


}; // end of class TimeOpWrap



/// \relates TimeOpWrap
template<class Operator>
Teuchos::RCP< TimeOpWrap<Operator> > createTimeOpWrap( const Teuchos::RCP<Operator>& op=Teuchos::null) {
  if( Teuchos::is_null(op) )
    return( Teuchos::rcp( new TimeOpWrap<Operator>( Teuchos::rcp( new Operator() ) ) ) );
  else
    return( Teuchos::rcp( new TimeOpWrap<Operator>( op ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMEOPWRAP_HPP
