#pragma once
#ifndef PIMPACT_OPERATORBASE_HPP
#define PIMPACT_OPERATORBASE_HPP


#include "Pimpact_OperatorMV.hpp"



namespace Pimpact {


/// \brief Operator class for type erasure
///
/// for Belos preconditioner and operator has to be from the same type
template<class MultiVector>
class OperatorBase {
public:

  typedef MultiVector MV;
  virtual ~OperatorBase() {};

//  virtual void apply( const MV& x, MV& y, Belos::ETrans trans=Belos::NOTRANS ) const =0 ;
//  virtual bool hasApplyTranspose() const =0;
  virtual void apply( const MV& x, MV& y, Belos::ETrans trans=Belos::NOTRANS ) const {} ;
  virtual bool hasApplyTranspose() const {return( false );};
}; // end of class OperatorBase


template<class MV,class Op>
class OperatorPimpl : public virtual OperatorBase<MV> {
  Teuchos::RCP<OperatorMV<Op> > opm_;
public:

  OperatorPimpl( const Teuchos::RCP<OperatorMV<Op> >& opm ):opm_(opm) {};

  virtual ~OperatorPimpl() {opm_=Teuchos::null;};

  virtual void apply( const MV& x, MV& y, Belos::ETrans trans=Belos::NOTRANS ) const {
    opm_->apply( x, y, trans );
  }

  virtual bool hasApplyTranspose() const {
    return( opm_->hasApplyTranspose() );
  };

  Teuchos::RCP<Op> getOperatorPtr() { return( opm_->getOperatorPtr() ); }

}; // end of OperatorPimpl


template<class MV, class Op>
Teuchos::RCP<OperatorBase<MV> > createOperatorBase( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null( op ) )
    return( Teuchos::rcp( new OperatorPimpl<MV,Op>( Pimpact::createOperatorMV<Op>() ) ) );
  else
    return( Teuchos::rcp( new OperatorPimpl<MV,Op>( Pimpact::createOperatorMV<Op>( op ) ) ) );
}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORBASE_HPP
