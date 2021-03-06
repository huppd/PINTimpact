/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_SIMPLEVECTORITERATION_HPP
#define PIMPACT_SIMPLEVECTORITERATION_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "BelosTypes.hpp"

//#include "Pimpact_Grid.hpp" // just for createOstream<>




namespace Pimpact {



template<class OperatorT>
class SimpleVectorIteration {

protected:

  using DomainFieldT = typename OperatorT::DomainFieldT;
  using RangeFieldT = typename OperatorT::RangeFieldT;

  using GridT = typename DomainFieldT::GridT;

  using ScalarT = typename GridT::Scalar;


  ScalarT lamMax_;

public:

  SimpleVectorIteration(
    const Teuchos::RCP<const OperatorT>& op,
    const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList()) {

    Teuchos::RCP<DomainFieldT> x = create<DomainFieldT>(op->grid());
    Teuchos::RCP<RangeFieldT>  r = create<RangeFieldT>(op->grid());

    x->random();

    int numIters = pl->get<int>("numIters", 200);
    ScalarT tol = pl->get<ScalarT>("tol", 1.e-3);

    ScalarT lamp = 0.;
    for(int i=0; i<numIters; ++i) {
      op->apply(*x, *r);
      ScalarT lam = r->dot(*x)/x->dot(*x);
      r->scale(1./r->norm());
      x.swap(r);
      if(std::fabs(lamp-lam)/std::fabs(lam) <tol)
        break;
      else
        lamp=lam;
    }
    //x->write(999);

    op->apply(*x, *r);
    lamMax_ = r->dot(*x)/x->dot(*x);
  }

  ScalarT getMaxEV() const {
    return lamMax_;
  };
  ScalarT getMinEV() const {
    return Teuchos::ScalarTraits<ScalarT>::zero();
  };

  bool computeMinEV() const {
    return false;
  };


}; // end of class SimpleVectorIteration


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SIMPLEVECTORITERATION_HPP
