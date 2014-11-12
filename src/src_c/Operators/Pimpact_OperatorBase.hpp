#pragma once
#ifndef PIMPACT_OPERATORBASE_HPP
#define PIMPACT_OPERATORBASE_HPP



namespace Pimpact {


/// \brief Operator base class for type erasure.
/// \ingroup Operator
template<class DomainField, class RangeField=DomainField>
class OperatorBase {

public:
  virtual ~OperatorBase() {};

  typedef DomainField MF;
  typedef DomainField DomainFieldT;
  typedef RangeField RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

  virtual void apply( const DomainField& x, RangeField& y, Belos::ETrans trans=Belos::NOTRANS ) const {} ;

  virtual void assignField( const DomainField& mv ) {};

  virtual bool hasApplyTranspose() const {return( false );};


}; // end of class OperatorBase



template<class Op>
class OperatorPimpl : public virtual OperatorBase<typename Op::DomainFieldT, typename Op::RangeFieldT> {

public:

  typedef typename Op::DomainFieldT DomainFieldT;
  typedef typename Op::RangeFieldT RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP<Op> opm_;

public:

  OperatorPimpl( const Teuchos::RCP<Op>& opm ):opm_(opm) {};

  virtual ~OperatorPimpl() {opm_=Teuchos::null;};

  virtual void apply( const DomainFieldT& x, RangeFieldT& y, Belos::ETrans trans=Belos::NOTRANS ) const {
    opm_->apply( x, y, trans );
  }

  virtual void assignField( const DomainFieldT& field ) {
    opm_->assignField( field );
  };

  virtual bool hasApplyTranspose() const {
    return( opm_->hasApplyTranspose() );
  };

  Teuchos::RCP<Op> getOperatorPtr() { return( opm_ ); }

}; // end of OperatorPimpl



/// \relates OperatorBase
/// \relates OperatorPimpl
template<class Op>
Teuchos::RCP< OperatorBase<typename Op::DomainFieldT, typename Op::RangeFieldT> >
createOperatorBase( const Teuchos::RCP<Op>& op ) {

    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<typename Op::DomainFieldT, typename Op::RangeFieldT>  >(
            Teuchos::rcp( new OperatorPimpl<Op>( op )
            )
        )
    );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORBASE_HPP
