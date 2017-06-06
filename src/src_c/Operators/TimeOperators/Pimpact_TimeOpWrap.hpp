#pragma once
#ifndef PIMPACT_TIMEOPWRAP_HPP
#define PIMPACT_TIMEOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_TimeField.hpp"
#include "Pimpact_Utils.hpp"




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

protected:

  using Ordinal = typename SpaceT::Ordinal;

  Teuchos::RCP<OperatorT> op_;

public:

  TimeOpWrap( const Teuchos::RCP<const SpaceT>& space ):
    op_( create<OperatorT>(space) ) {};

  TimeOpWrap( const Teuchos::RCP<OperatorT>& op ): op_(op) {};

  /// \brief default apply
  void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N ) const {

    if( true==CNyes ) {

      typename OperatorT::RangeFieldT temp( space() );

      x.exchange();
      for( int i=0; i <=space()->ei(F::S,3); ++i ) {

        op_->apply( x(i), temp );

        if( i>=space()->si(F::S,3) )
          y(i).add( 1., y(i), 0.5, temp );

        if( i+1<=space()->ei(F::S,3) )
          y(i+1).add( 0., y(i+1), 0.5, temp );
      }
    } else {
      for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
        op_->apply( x(i) , y(i), add );
    }
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const {
    return( op_->hasApplyTranspose() );
  }

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return(op_->space());
  };

  Teuchos::RCP<OperatorT> getOperatorPtr() {
    return( op_ );
  }
  void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
    op_->setParameter( para );
  }

  void print( std::ostream& out=std::cout ) const {
    op_->print(out);
  }

  const std::string getLabel() const {
    return( "TimeOpWrap (" );
  };

}; // end of class TimeOpWrap



/// \relates TimeOpWrap
template< class OpT, bool CNY=false >
Teuchos::RCP< TimeOpWrap<OpT,CNY> > createTimeOpWrap(
  const Teuchos::RCP<OpT>& op ) {

  return( Teuchos::rcp( new TimeOpWrap<OpT,CNY>( op ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMEOPWRAP_HPP
