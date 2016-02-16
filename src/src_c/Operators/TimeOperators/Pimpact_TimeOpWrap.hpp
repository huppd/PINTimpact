#pragma once
#ifndef PIMPACT_TIMEOPWRAP_HPP
#define PIMPACT_TIMEOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_TimeField.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact {


/// \brief Operator wrapper for \c TimeField's.
/// \ingroup TimeOperator
/// wraps and \c OperatorT and adds the functionality of handling \c TimeField.
template<class OperatorT, bool CNyes=false>
class TimeOpWrap  {

public:

  using DomainFieldT = TimeField<typename OperatorT::DomainFieldT>;
  using RangeFieldT = TimeField<typename OperatorT::RangeFieldT>;

  using SpaceT = typename DomainFieldT::SpaceT;

  using Ordinal = typename SpaceT::Ordinal;

protected:

  Teuchos::RCP<OperatorT> op_;
  Teuchos::RCP<typename OperatorT::RangeFieldT> temp_;

public:

  TimeOpWrap(
			const Teuchos::RCP<const SpaceT>& space ):
		op_( create<OperatorT>(space) ),
		temp_( create<typename OperatorT::RangeFieldT>(space) ) {};

  TimeOpWrap(
			const Teuchos::RCP<OperatorT>& op ):
		op_(op),
		temp_( create<typename OperatorT::RangeFieldT>(op->space()) ) {};

  /// \brief default apply
  void apply( const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {

    if( true==CNyes ) {

      x.exchange();
      for( int i=0; i < space()->eInd(S,3); ++i ) {

        op_->apply( x.getConstField(i), *temp_ );

        if( i>=space()->sInd(S,3) )
          y.getFieldPtr(i)->add( 1., y.getConstField(i), 0.5, *temp_ );

        if( i+1<space()->eInd(S,3) )
          y.getFieldPtr(i+1)->add( 0., y.getConstField(i+1), 0.5, *temp_ );
      }

    }
    else{

//      typename RangeFieldT::Iter j = y.sInd_;
      for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
        op_->apply( x.getConstField(i) , y.getField(i) );

    }
    y.changed();

  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

  Teuchos::RCP<OperatorT> getOperatorPtr() { return( op_ ); }
	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}

  void print( std::ostream& out=std::cout ) const {
		op_->print(out);
	}

	const std::string getLabel() const { return( "TimeOpWrap (" ); };

}; // end of class TimeOpWrap



/// \relates TimeOpWrap
template< class OpT, bool CNY=false >
Teuchos::RCP< TimeOpWrap<OpT,CNY> > createTimeOpWrap(
		const Teuchos::RCP<OpT>& op ) {

	return( Teuchos::rcp( new TimeOpWrap<OpT,CNY>( op ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMEOPWRAP_HPP
