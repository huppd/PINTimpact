#pragma once
#ifndef PIMPACT_MULTIOPUNWRAP_HPP
#define PIMPACT_MULTIOPUNWRAP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_MultiField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief Operator wrapper for \c MultiField's.
/// \ingroup MultiOperator
/// wraps and \c Operator and adds the functionality of handling \c MultiField.
template<class MOperatorT>
class MultiOpUnWrap  {

public:

  using OperatorT = MOperatorT;

  using DomainFieldT = typename MOperatorT::DomainFieldT::InnerFieldT;
  using RangeFieldT = typename MOperatorT::RangeFieldT::InnerFieldT;

  using SpaceT = typename DomainFieldT::SpaceT;

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
  void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N ) const {

    op_->apply(
      *wrapMultiField( Teuchos::rcp_const_cast<DomainFieldT>(Teuchos::rcpFromRef(x)) ),
      *wrapMultiField( Teuchos::rcpFromRef(y) ) );
  }

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( MultiField<DomainFieldT>( Teuchos::rcp_const_cast<DomainFieldT>(Teuchos::rcpFromRef(mv))) );
  };

  bool hasApplyTranspose() const {
    return( op_->hasApplyTranspose() );
  }

  Teuchos::RCP<MOperatorT> getOperatorPtr() {
    return( op_ );
  }

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return(op_->space());
  }

  void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
    op_->setParameter( para );
  }

  const std::string getLabel() const {
    return( op_->getLabel() );
  };

  void print( std::ostream& out=std::cout ) const {
    op_->print();
  }

}; // end of class MultiOpUnWrap



/// \relates MultiOpUnWrap
template<class MOperatorT>
Teuchos::RCP< MultiOpUnWrap<MOperatorT> > createMultiOpUnWrap( const Teuchos::RCP<MOperatorT>& op ) {

  return( create<MultiOpUnWrap,MOperatorT>( op ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIOPUNWRAP_HPP
