#pragma once
#ifndef PIMPACT_MULTIOPWRAP_HPP
#define PIMPACT_MULTIOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_MultiField.hpp"




namespace Pimpact {


/// \brief Operator wrapper for \c MultiField's.
/// \ingroup MultiOperator
/// wraps and \c Operator and adds the functionality of handling \c MultiField.
template<class OT>
class MultiOpWrap  {

public:

  typedef OT OperatorT;

  typedef typename OperatorT::SpaceT SpaceT;

  typedef MultiField<typename OperatorT::DomainFieldT> DomainFieldT;
  typedef MultiField<typename OperatorT::RangeFieldT> RangeFieldT;

protected:

  Teuchos::RCP<OperatorT> op_;

public:

  MultiOpWrap( const Teuchos::RCP<OperatorT>& op ):op_(op) {}

  template<class IOperatorT>
  MultiOpWrap( const Teuchos::RCP<IOperatorT>& op ):
    op_( create<IOperatorT>(op) ) {}


  /// \brief default apply
  void apply( const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS) const {

    TEUCHOS_TEST_FOR_EXCEPT( x.getNumberVecs()!=y.getNumberVecs() );

    int m = x.getNumberVecs();

    for( int i=0; i<m; ++i )
      op_->apply( x.getConstField(i), y.getField(i) );
  }

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( mv.getConstField(0) );
  };

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

  Teuchos::RCP<OperatorT> getOperatorPtr() { return( op_ ); }

  Teuchos::RCP<const SpaceT> getSpace() const { return( op_->getSpace() ); };
  Teuchos::RCP<const SpaceT> space() const { return( op_->space() ); };

}; // end of class MultiOpWrap



/// \relates MultiOpWrap
/// \todo offer version where create<SpaceObject> is used
template<class OperatorT>
Teuchos::RCP< MultiOpWrap<OperatorT> > createMultiOpWrap( const Teuchos::RCP<OperatorT>& op ) {

    return( Teuchos::rcp( new MultiOpWrap<OperatorT>( op ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIOPWRAP_HPP
