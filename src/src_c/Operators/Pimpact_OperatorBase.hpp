#pragma once
#ifndef PIMPACT_OPERATORBASE_HPP
#define PIMPACT_OPERATORBASE_HPP



namespace Pimpact {


/// \brief Operator base class for type erasure.
/// \ingroup Operator
template<class DomainField,class RangeField=DomainField>
class OperatorBase {

public:
  virtual ~OperatorBase() {};

  typedef DomainField MF;
  typedef DomainField DomainFieldT;
  typedef RangeField RangeFieldT;

  virtual void apply( const DomainField& x, RangeField& y, Belos::ETrans trans=Belos::NOTRANS ) const {} ;

  virtual void assignField( const DomainField& mv ) {};

  virtual bool hasApplyTranspose() const {return( false );};


}; // end of class OperatorBase



template< class DomainField, class Op, class RangeField=DomainField >
class OperatorPimpl : public virtual OperatorBase<DomainField,RangeField> {
  Teuchos::RCP<Op> opm_;
public:
  typedef DomainField DomainFieldT;
  typedef RangeField RangeFieldT;

  OperatorPimpl( const Teuchos::RCP<Op>& opm ):opm_(opm) {};

  virtual ~OperatorPimpl() {opm_=Teuchos::null;};

  virtual void apply( const DomainField& x, RangeField& y, Belos::ETrans trans=Belos::NOTRANS ) const {
    opm_->apply( x, y, trans );
  }

  virtual void assignField( const DomainField& field ) {
    opm_->assignField( field );
  };

  virtual bool hasApplyTranspose() const {
    return( opm_->hasApplyTranspose() );
  };




  Teuchos::RCP<Op> getOperatorPtr() { return( opm_ ); }

}; // end of OperatorPimpl



/// \relates OperatorBase
/// \relates OperatorPimpl
template<class DF, class Op, class RF=DF>
Teuchos::RCP<OperatorBase<DF,RF> > createOperatorBase( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null( op ) )
    return( Teuchos::rcp( new OperatorPimpl<DF,Op,RF>( Teuchos::rcp(new Op()) ) ) );
  else
    return( Teuchos::rcp( new OperatorPimpl<DF,Op,RF>( op ) ) );
}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORBASE_HPP
