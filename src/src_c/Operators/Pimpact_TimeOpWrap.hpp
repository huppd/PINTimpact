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
template<class Operator, bool CNyes=false>
class TimeOpWrap  {

public:

  typedef TimeField<typename Operator::DomainFieldT> DomainFieldT;
  typedef TimeField<typename Operator::RangeFieldT> RangeFieldT;

protected:

  Teuchos::RCP<Operator> op_;
  Teuchos::RCP<typename Operator::RangeFieldT> temp_;

public:

  TimeOpWrap():op_( Teuchos::rcp( new Operator() ) ) {};

  TimeOpWrap(
      const Teuchos::RCP<Operator>& op,
      const Teuchos::RCP<typename Operator::RangeFieldT>& temp=Teuchos::null ):
        op_(op),temp_(temp) {};

//  ~TimeOpWrap() { op_=Teuchos::null; };


  /// \brief default apply
  void apply( const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS) const {

    if( true==CNyes ) {


      const_cast<DomainFieldT&>(x).exchange();
//      y.init(0.);

      typename DomainFieldT::Iter i= const_cast<DomainFieldT&>(x).mfs_.begin();
      typename RangeFieldT::Iter  j = y.mfs_.begin();

      for( ; i < const_cast<DomainFieldT&>(x).mfs_.end(); ++i ) {

        op_->apply( **(i), *temp_ );

        if( j>=y.beginI_ )
          (*j)->add( 1., **j, 0.5, *temp_ );
        ++j;
        if( j<y.endI_ )
          (*j)->add( 0., **j, 0.5, *temp_ );
      }

    }
    else{
      typename RangeFieldT::Iter j = y.beginI_;
      for( typename DomainFieldT::Iter i=x.beginI_; i<x.endI_; ++i )
        op_->apply( **i , **(j++) );
    }
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( mv.getConstField(0) );
  };

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

  Teuchos::RCP<Operator> getOperatorPtr() { return( op_ ); }


}; // end of class TimeOpWrap



/// \relates TimeOpWrap
template< class Operator, bool CNY=false >
Teuchos::RCP< TimeOpWrap<Operator,CNY> > createTimeOpWrap(
    const Teuchos::RCP<Operator>& op=Teuchos::null,
    const Teuchos::RCP<typename Operator::RangeFieldT>& temp=Teuchos::null ) {
  if( Teuchos::is_null(op) )
    return( Teuchos::rcp( new TimeOpWrap<Operator,CNY>( Teuchos::rcp(new Operator()), temp ) ) );
  else
    return( Teuchos::rcp( new TimeOpWrap<Operator,CNY>( op, temp ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMEOPWRAP_HPP
