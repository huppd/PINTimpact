#pragma once
#ifndef PIMPACT_OPERATORBASE_HPP
#define PIMPACT_OPERATORBASE_HPP


#include "Pimpact_OperatorMV.hpp"



namespace Pimpact {


/// \brief Operator base class for type erasure.
/// \ingroup Operator
template<class MultiVector>
class OperatorBase {
public:

  typedef MultiVector MV;
  virtual ~OperatorBase() {};

//  virtual void apply( const MV& x, MV& y, Belos::ETrans trans=Belos::NOTRANS ) const =0 ;
//  virtual bool hasApplyTranspose() const =0;
  virtual void apply( const MV& x, MV& y, Belos::ETrans trans=Belos::NOTRANS ) const {} ;
//  virtual void setVector( const Teuchos::RCP<MV>& x ) {};
  virtual bool hasApplyTranspose() const {return( false );};

}; // end of class OperatorBase



///// \deprecated
//template<class MV,class Op>
//class OperatorPimpldep : public virtual OperatorBase<MV> {
//  Teuchos::RCP<OperatorMV<Op> > opm_;
//public:
//
//  OperatorPimpldep( const Teuchos::RCP<OperatorMV<Op> >& opm ):opm_(opm) {};
//
//  virtual ~OperatorPimpldep() {opm_=Teuchos::null;};
//
//  virtual void apply( const MV& x, MV& y, Belos::ETrans trans=Belos::NOTRANS ) const {
//    opm_->apply( x, y, trans );
//  }
//
////  virtual void setVector( const Teuchos::RCP<MV>& x ) { opm_->set( x ) };
//
//  virtual bool hasApplyTranspose() const {
//    return( opm_->hasApplyTranspose() );
//  };
//
//  Teuchos::RCP<Op> getOperatorPtr() { return( opm_->getOperatorPtr() ); }
//
//}; // end of OperatorPimpl
//
//
///// \deprecated
//template<class MV, class Op>
//Teuchos::RCP<OperatorBase<MV> > createOperatorBasedep( const Teuchos::RCP<Op>& op=Teuchos::null ) {
//  if( Teuchos::is_null( op ) )
//    return( Teuchos::rcp( new OperatorPimpldep<MV,Op>( Pimpact::createOperatorMV<Op>() ) ) );
//  else
//    return( Teuchos::rcp( new OperatorPimpldep<MV,Op>( Pimpact::createOperatorMV<Op>( op ) ) ) );
//}



template<class MV,class Op>
class OperatorPimpl : public virtual OperatorBase<MV> {
  Teuchos::RCP<Op> opm_;
public:

  OperatorPimpl( const Teuchos::RCP<Op>& opm ):opm_(opm) {};

  virtual ~OperatorPimpl() {opm_=Teuchos::null;};

  virtual void apply( const MV& x, MV& y, Belos::ETrans trans=Belos::NOTRANS ) const {
    opm_->apply( x, y, trans );
  }

//  virtual void setVector( const Teuchos::RCP<MV>& x ) { opm_->set( x ) };

  virtual bool hasApplyTranspose() const {
    return( opm_->hasApplyTranspose() );
  };

  Teuchos::RCP<Op> getOperatorPtr() { return( opm_ ); }

}; // end of OperatorPimpl


template<class MV, class Op>
Teuchos::RCP<OperatorBase<MV> > createOperatorBase( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null( op ) )
    return( Teuchos::rcp( new OperatorPimpl<MV,Op>( Teuchos::rcp(new Op()) ) ) );
  else
    return( Teuchos::rcp( new OperatorPimpl<MV,Op>( op ) ) );
}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORBASE_HPP
