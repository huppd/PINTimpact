#pragma once
#ifndef PIMPACT_FORCINGOP_HPP
#define PIMPACT_FORCINGOP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief forcing operator.
/// \ingroup Operator
template<class Field>
class ForcingOp {

public:

  using DomainFieldT = Field;
  using RangeFieldT = Field;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  using Scalar = typename SpaceT::Scalar;

  Teuchos::RCP<DomainFieldT> forcing_;
  Scalar mul_;

public:

  ForcingOp(const Teuchos::RCP<DomainFieldT>& forcing=Teuchos::null, Scalar mul=1):
    forcing_(forcing), mul_(mul) {};

  void setForcing(const Teuchos::RCP<DomainFieldT>& forcing) {
    forcing_ = forcing;
  }

  void setMultiplicator(const Scalar mul) {
    mul_ = mul;
  }

  /// \brief \f[ y = force*x \f]
  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    if(std::abs(mul_-1.) <Teuchos::ScalarTraits<Scalar>::eps()) {
      y = x;
    } else
      y.add(mul_, x, 0., y);
    y.scale(*forcing_);
  }

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return forcing_->space();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  void assignField(const DomainFieldT& field) {};

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "Forcing";
  };

  void print(std::ostream& out=std::cout) const {
    out <<getLabel() <<":\n";
  }

}; // end of ForcingOp



/// \relates ForcingOp
template<class F>
Teuchos::RCP<ForcingOp<F> > createForcingOp(
    const Teuchos::RCP<F>& forcing, typename F::SpaceT::Scalar mul=1.) {

  return Teuchos::rcp(new ForcingOp<F>(forcing, mul));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_FORCINGOP_HPP
