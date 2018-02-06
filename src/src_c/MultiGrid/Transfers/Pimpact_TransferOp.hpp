#pragma once
#ifndef PIMPACT_TRANSFEROP_HPP
#define PIMPACT_TRANSFEROP_HPP


#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief Transfers fields from "coarse" to "fine" grids, necessary when \c Grid::dimNC is  different.
///
/// Goes in both direction. If this is used a lot, it could be beneficial, to
/// seperate StencilWidths with data layout, so having same datalayout for all
/// stencile. Coping could be beneficial because Cash effects are bether
///
/// \tparam FGridT fine grid type in the sense of stencil order.
/// \tparam CGridT coase grid type
/// \ingroup BaseOperator
template<class FGT, class CGT>
class TransferOp {

public:

  using GridT = FGT;

  using FGridT = FGT;
  using CGridT = CGT;

  using Scalar = typename FGridT::Scalar;
  using Ordinal = typename FGridT::Ordinal;

  using DomainFieldT = ScalarField<FGridT>;
  using RangeFieldT = ScalarField<CGridT>;

protected:

  using TO = const Teuchos::Tuple<Scalar*, 3>;

  Teuchos::RCP<const FGridT> fGrid_;
  Teuchos::RCP<const CGridT> cGrid_;

public:

  TransferOp(
    const Teuchos::RCP<const FGridT>& fGrid,
    const Teuchos::RCP<const CGridT>& cGrid):
    fGrid_(fGrid), cGrid_(cGrid) {}


  template<class SP1T, class SP2T>
  void apply(const ScalarField<SP1T>& x, ScalarField<SP2T>& y) const {

    const F fType = x.getType();

    assert(fType == y.getType());

    for(int i=0; i<GridT::sdim; ++i)
      assert(x.grid()->nLoc(i) == y.grid()->nLoc(i));

    for(Ordinal k=x.grid()->si(fType, Z, B::Y); k<=x.grid()->ei(fType, Z, B::Y); ++k)
      for(Ordinal j=x.grid()->si(fType, Y, B::Y); j<=x.grid()->ei(fType, Y, B::Y); ++j)
        for(Ordinal i=x.grid()->si(fType, X, B::Y); i<=x.grid()->ei(fType, X, B::Y); ++i)
          y(i, j, k) = x(i, j, k);

    y.changed();
  }


  void assignField(const RangeFieldT& mv) {};

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  bool hasApplyTranspose() const {
    return false;
  }

  void print(std::ostream& out=std::cout) const {
    out << "--- " << getLabel() << " ---\n";
  }

  const std::string getLabel() const {
    return "TransferOp";
  };

}; // end of class TransferOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRANSFEROP_HPP
