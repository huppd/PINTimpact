/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_COMPOUNDSMOOTHER_HPP
#define PIMPACT_COMPOUNDSMOOTHER_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_TripleCompositionOp.hpp"




namespace Pimpact {



/// \ingroup CompoundOperator
///
/// \tparam CompoundOpT should be of compatible to CompoundOpWrap
/// \f[ \begin{bmatrix} opV2V^{-1} & opS2V \\ 0 & -opS2S^{-1} \end{bmatrix}^{-1} \mathbf{x} = \mathbf{y} \f]
template<class CompoundOpT, template<class> class vSmoother, class sSmoother>
class CompoundSmoother {

public:

  using OpV2VT = typename CompoundOpT::OpV2VT;
  using OpS2VT = typename CompoundOpT::OpS2VT;
  using OpV2ST = typename CompoundOpT::OpV2ST;

  using DomainFieldT = typename CompoundOpT::DomainFieldT;
  using RangeFieldT = typename CompoundOpT::RangeFieldT;

  using GridT = typename CompoundOpT::GridT;

protected:

  using VF = typename OpS2VT::RangeFieldT;
  using SF = typename OpS2VT::DomainFieldT;

  using OpVSmoother = vSmoother<OpV2VT>;
  using OpSSmootherT = TripleCompositionOp<sSmoother, TripleCompositionOp<OpV2ST, OpV2VT, OpS2VT>, sSmoother>;

  const Teuchos::RCP<OpS2VT> opS2V_;

  const Teuchos::RCP<OpVSmoother> opVSmoother_;
  Teuchos::RCP<OpSSmootherT> opSSmoother_;

public:

  CompoundSmoother(
    const Teuchos::RCP<CompoundOpT >& op,
    Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::null):
    opS2V_(op->getOpS2V()),
    opVSmoother_(Teuchos::rcp(new OpVSmoother(op->getOpV2V(),
                                Teuchos::sublist(pl, "VSmoother")))) {

    Teuchos::RCP<TripleCompositionOp<OpV2ST, OpV2VT, OpS2VT> > opDHG =
      createTripleCompositionOp(op->getOpV2S(), op->getOpV2V(),
                                 op->getOpS2V());

    Teuchos::RCP<sSmoother> opDGi = Teuchos::rcp(new sSmoother(grid()));

    opSSmoother_ = createTripleCompositionOp(opDGi, opDHG, opDGi);
  };

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    Teuchos::RCP<VF> tempv =  create<VF>(grid());

    opSSmoother_->apply(x.getSField() ,  y.getSField());

    opS2V_->apply(y.getSField(), *tempv);

    tempv->add(-1., *tempv, 1., x.getVField());

    opVSmoother_->apply(*tempv, y.getVField());
  }


  void assignField(const DomainFieldT& mv) { };

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return opS2V_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const std::string getLabel() const {
    return "CompoundSmoother";
  };

  void print(std::ostream& out=std::cout) const {
    out << getLabel() << ":\n";
//    opV2V_->print(out);
//		opS2V_->print(out);
//		opS2S_->print(out);
  }

}; // end of class CompoundSmoother



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDSMOOTHER_HPP
