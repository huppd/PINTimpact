#pragma once
#ifndef PIMPACT_MULTIOPUNWRAP_HPP
#define PIMPACT_MULTIOPUNWRAP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_MultiField.hpp"




namespace Pimpact {


/// \brief Operator wrapper for \c MultiField's.
/// \ingroup MultiOperator
/// wraps and \c Operator and adds the functionality of handling \c MultiField.
template<class MOperatorT>
class MultiOpUnWrap  {

public:

  typedef MOperatorT OperatorT;

  typedef typename MOperatorT::DomainFieldT::FieldT DomainFieldT;
  typedef typename MOperatorT::RangeFieldT::FieldT RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP<MOperatorT> op_;

public:

  MultiOpUnWrap( const Teuchos::RCP<MOperatorT>& op ):
    op_(op) {}

  template<class IOperatorT>
  MultiOpUnWrap( const Teuchos::RCP<IOperatorT>& op ):
    op_( create<MOperatorT>(op) ) {}

  template<class IOperatorT>
  MultiOpUnWrap( const Teuchos::RCP<IOperatorT>& op, Teuchos::RCP<Teuchos::ParameterList> pl ):
    op_( Teuchos::rcp( new MOperatorT( op, pl ) ) ) {}


  /// \brief default apply
  void apply( const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS) const {

    op_->apply(
        *createMultiField( Teuchos::rcp_const_cast<DomainFieldT>(Teuchos::rcpFromRef(x)) ),
        *createMultiField( Teuchos::rcpFromRef(y) ) );

  }

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( *createMultiField( Teuchos::rcp_const_cast<DomainFieldT>(Teuchos::rcpFromRef(mv))) );
  };

  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

  Teuchos::RCP<MOperatorT> getOperatorPtr() { return( op_ ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); }

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}


}; // end of class MultiOpUnWrap



/// \relates MultiOpUnWrap
/// \todo offer version where create<SpaceObject> is used
template<class MOperatorT>
Teuchos::RCP< MultiOpUnWrap<MOperatorT> > createMultiOpUnWrap( const Teuchos::RCP<MOperatorT>& op ) {

//    return( Teuchos::rcp( new MultiOpUnWrap<MOperatorT>( op ) ) );
    return( create<MultiOpUnWrap,MOperatorT>( op ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIOPUNWRAP_HPP
