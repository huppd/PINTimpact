#pragma once
#ifndef PIMPACT_INVERSETRIANGULAROP_HPP
#define PIMPACT_INVERSETRIANGULAROP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_MultiField.hpp"




namespace Pimpact {



/// \ingroup CompoundOperator
///
/// \f[ \begin{bmatrix} opV2V^{-1} & opS2V \\ 0 & -opS2S^{-1} \end{bmatrix}^{-1} \mathbf{x} = \mathbf{y} \f]
template<class OpV2V, class OpS2V, class OpS2S>
class InverseTriangularOp {

  using VF = typename OpV2V::RangeFieldT;
  using SF = typename OpS2S::RangeFieldT;

public:

  using DomainFieldT = CompoundField<VF, SF>;
  using RangeFieldT = CompoundField<VF, SF>;

  using GridT = typename DomainFieldT::GridT;

protected:

  const Teuchos::RCP<OpV2V> opV2V_;
  const Teuchos::RCP<OpS2V> opS2V_;
  const Teuchos::RCP<OpS2S> opS2S_;

public:

  InverseTriangularOp(
    const Teuchos::RCP<OpV2V>& opV2V,
    const Teuchos::RCP<OpS2V>& opS2V,
    const Teuchos::RCP<OpS2S>& opS2S):
    opV2V_(opV2V),
    opS2V_(opS2V),
    opS2S_(opS2S) {};


  /// \brief apply
  ///
  /// \f[ \begin{bmatrix} opV2V^{-1} & opS2V \\ 0 & -opS2S^{-1} \end{bmatrix}^{-1} \mathbf{x} = \mathbf{y} \f]
  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    // upper triangular
    opS2S_->apply(x.getSField(),  y.getSField());
    y.getSField().scale(-1.);

    VF tempv(grid());

    opS2V_->apply(y.getSField(), tempv, Add::Y);

    tempv.add(-1., tempv, 1., x.getVField());

    opV2V_->apply(tempv, y.getVField());
  }


  void assignField(const DomainFieldT& mv) {
    opV2V_->assignField(mv.getVField());
  };

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return opV2V_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {
    opV2V_->setParameter(para);
    opS2S_->setParameter(para);
  }

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "InverseTriangularOp ";
  };

  void print(std::ostream& out=std::cout) const {
    out << getLabel() << ":\n";
    opV2V_->print(out);
    opS2V_->print(out);
    opS2S_->print(out);
  }

}; // end of class InverseTriangularOp



/// \relates InverseTriangularOp
template<class OpV2V, class OpS2V, class OpS2S >
Teuchos::RCP<InverseTriangularOp<OpV2V, OpS2V, OpS2S> >
createInverseTriangularOp(
    const Teuchos::RCP<OpV2V>& opV2V,
    const Teuchos::RCP<OpS2V>& opS2V,
    const Teuchos::RCP<OpS2S>& opS2S) {

  //  return Teuchos::null;
  return Teuchos::rcp(new InverseTriangularOp<OpV2V, OpS2V, OpS2S>(opV2V, opS2V, opS2S));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSETRIANGULAROP_HPP
