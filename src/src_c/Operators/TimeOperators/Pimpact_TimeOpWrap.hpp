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

  typedef typename DomainFieldT::SpaceT SpaceT;

  typedef typename SpaceT::Ordinal Ordinal;

protected:

  Teuchos::RCP<Operator> op_;
  Teuchos::RCP<typename Operator::RangeFieldT> temp_;

public:

  TimeOpWrap(
      const Teuchos::RCP<Operator>& op,
      const Teuchos::RCP<typename Operator::RangeFieldT>& temp ):
        op_(op),temp_(temp) {};

  /// \brief default apply
  void apply( const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {

    if( true==CNyes ) {

      x.exchange();
//      y.init(0.);

//      typename DomainFieldT::Iter i= const_cast<DomainFieldT&>(x).mfs_.begin();
//      typename RangeFieldT::Iter  j = y.mfs_.begin();

//      for( ; i < const_cast<DomainFieldT&>(x).mfs_.end(); ++i ) {
      for( int i=0; i < op_->getSpace()->eInd(S,3); ++i ) {

//        op_->apply( **(i), *temp_ );
        op_->apply( x.getConstField(i), *temp_ );

//        if( j>=y.sInd_ )
        if( i>=op_->getSpace()->sInd(S,3) )
//          (*j)->add( 1., **j, 0.5, *temp_ );
          y.getFieldPtr(i)->add( 1., y.getConstField(i), 0.5, *temp_ );
//        ++j;
//        if( j<y.eInd_ )
        if( i+1<op_->getSpace()->eInd(S,3) )
          y.getFieldPtr(i+1)->add( 0., y.getConstField(i+1), 0.5, *temp_ );
      }

    }
    else{
//      typename RangeFieldT::Iter j = y.sInd_;
//      for( typename DomainFieldT::Iter i=x.sInd_; i<x.eInd_; ++i )
//        op_->apply( **i , **(j++) );
      typename RangeFieldT::Iter j = y.sInd_;
      for( Ordinal i=op_->getSpace()->sInd(S,3); i<op_->getSpace()->eInd(S,3); ++i )
        op_->apply( x.getConstField(i) , y.getField(i) );
    }
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

  Teuchos::RCP<Operator> getOperatorPtr() { return( op_ ); }


}; // end of class TimeOpWrap



/// \relates TimeOpWrap
template< class Operator, bool CNY=false >
Teuchos::RCP< TimeOpWrap<Operator,CNY> > createTimeOpWrap(
    const Teuchos::RCP<Operator>& op,
    const Teuchos::RCP<typename Operator::RangeFieldT>& temp ) {

    return( Teuchos::rcp( new TimeOpWrap<Operator,CNY>( op, temp ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMEOPWRAP_HPP
