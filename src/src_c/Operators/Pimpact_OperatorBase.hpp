#pragma once
#ifndef PIMPACT_OPERATORBASE_HPP
#define PIMPACT_OPERATORBASE_HPP


//#include "Pimpact_OperatorMV.hpp"



namespace Pimpact {


/// \brief Operator base class for type erasure.
/// \ingroup Operator
template<class MultiVector>
class OperatorBase {

public:
  virtual ~OperatorBase() {};

  typedef MultiVector MF;

  virtual void apply( const MF& x, MF& y, Belos::ETrans trans=Belos::NOTRANS ) const {} ;

  virtual void assignField( const MF& mv ) {};

  virtual bool hasApplyTranspose() const {return( false );};


}; // end of class OperatorBase



template<class MF,class Op>
class OperatorPimpl : public virtual OperatorBase<MF> {
  Teuchos::RCP<Op> opm_;
public:

  OperatorPimpl( const Teuchos::RCP<Op>& opm ):opm_(opm) {};

  virtual ~OperatorPimpl() {opm_=Teuchos::null;};

  virtual void apply( const MF& x, MF& y, Belos::ETrans trans=Belos::NOTRANS ) const {
    opm_->apply( x, y, trans );
  }

  virtual void assignField( const MF& mv ) {
    opm_->assignField( mv );
  };

  virtual bool hasApplyTranspose() const {
    return( opm_->hasApplyTranspose() );
  };




  Teuchos::RCP<Op> getOperatorPtr() { return( opm_ ); }

}; // end of OperatorPimpl



template<class MF, class Op>
Teuchos::RCP<OperatorBase<MF> > createOperatorBase( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null( op ) )
    return( Teuchos::rcp( new OperatorPimpl<MF,Op>( Teuchos::rcp(new Op()) ) ) );
  else
    return( Teuchos::rcp( new OperatorPimpl<MF,Op>( op ) ) );
}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORBASE_HPP
