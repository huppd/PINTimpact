#pragma once
#ifndef PIMPACT_OPERATORBASE_HPP
#define PIMPACT_OPERATORBASE_HPP


#include "Pimpact_OperatorMV.hpp"



namespace Pimpact {


/// \brief Operator class for type erasure
///
/// for Belos preconditioner and operator has to be from the same type
template<class MV>
class OperatorBase {
public:
  virtual ~OperatorBase() {};

  virtual void apply( const MV& x, MV& y, Belos::ETrans trans=Belos::NOTRANS ) const =0 ;
  virtual bool hasApplyTranspose() const =0;
};


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
};


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORBASE_HPP
