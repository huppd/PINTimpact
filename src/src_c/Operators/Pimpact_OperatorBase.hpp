/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_OPERATORBASE_HPP
#define PIMPACT_OPERATORBASE_HPP

#include "BelosTypes.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"




namespace Pimpact {


/// \brief Operator base class for type erasure.
/// \ingroup Operator
template<class DomainField, class RangeField=DomainField>
class OperatorBase {

public:
  virtual ~OperatorBase() {};

  using MF = DomainField;
  using DomainFieldT = DomainField;
  using RangeFieldT = RangeField;

  using GridT = typename DomainFieldT::GridT;

  virtual void apply(const DomainField& x, RangeField& y, const Belos::ETrans
      trans=Belos::NOTRANS) const {} ;

  virtual void assignField(const DomainField& mv) {};

  virtual const Teuchos::RCP<const GridT>& grid() const =0;

  virtual void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {}

  virtual bool hasApplyTranspose() const {
    return false;
  };

  virtual const std::string getLabel() const {
    return std::string("PImpact: ");
  };

  virtual void print(std::ostream& out=std::cout) const {  };

}; // end of class OperatorBase



template<class Op>
class OperatorPimpl : public virtual OperatorBase<typename Op::DomainFieldT, typename Op::RangeFieldT> {

public:

  using DomainFieldT = typename Op::DomainFieldT;
  using  RangeFieldT = typename Op::RangeFieldT;

  using GridT = typename DomainFieldT::GridT;

protected:

  Teuchos::RCP<Op> opm_;

public:

  OperatorPimpl(const Teuchos::RCP<Op>& opm):opm_(opm) {};

  virtual ~OperatorPimpl() {
    opm_=Teuchos::null;
  };

  virtual void apply(const DomainFieldT& x, RangeFieldT& y, const Belos::ETrans
      trans=Belos::NOTRANS) const {

    opm_->apply(x, y, trans);
  }

  virtual void assignField(const DomainFieldT& field) {
    opm_->assignField(field);
  };

  virtual bool hasApplyTranspose() const {
    return opm_->hasApplyTranspose();
  };

  virtual const Teuchos::RCP<const GridT>& grid() const {
    return opm_->grid();
  };

  virtual void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    opm_->setParameter(para);
  }

  virtual Teuchos::RCP<Op> getOperatorPtr() {
    return opm_;
  }

  virtual const std::string getLabel() const {
    return opm_->getLabel();
  };

  virtual void print(std::ostream& out=std::cout) const {
    opm_->print(out);
  };

}; // end of OperatorPimpl



/// \relates OperatorBase
/// \relates OperatorPimpl
template<class Op>
Teuchos::RCP<OperatorBase<typename Op::DomainFieldT, typename Op::RangeFieldT> >
createOperatorBase(const Teuchos::RCP<Op>& op) {

  return Teuchos::rcp_dynamic_cast<OperatorBase<typename Op::DomainFieldT, typename
    Op::RangeFieldT> >(Teuchos::rcp(new OperatorPimpl<Op>(op)));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_OPERATORBASE_HPP
