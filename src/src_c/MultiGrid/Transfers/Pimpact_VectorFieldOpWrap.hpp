#pragma once
#ifndef PIMPACT_VECTORFIELDOPWRAP_HPP
#define PIMPACT_VECTORFIELDOPWRAP_HPP


#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




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
template<class SOpT>
class VectorFieldOpWrap {

public:


  using FGridT = typename SOpT::FGridT;
  using CGridT = typename SOpT::CGridT;

  using GridT = FGridT;

  using DomainFieldT = VectorField<FGridT>;
  using  RangeFieldT = VectorField<CGridT>;

protected:

  Teuchos::RCP<const SOpT> sop_;

public:

  template<class SP1T, class SP2T>
  VectorFieldOpWrap(
    const Teuchos::RCP<const SP1T>& fGrid,
    const Teuchos::RCP<const SP2T>& cGrid):
    sop_(create<SOpT>(fGrid, cGrid)) {}

  template<class SP1T, class SP2T>
  VectorFieldOpWrap(
    const Teuchos::RCP<const SP1T>& fGrid,
    const Teuchos::RCP<const SP2T>& cGrid,
    Teuchos::Tuple<int, SP1T::dimension> nb):
    sop_(Teuchos::rcp(new SOpT(fGrid, cGrid, nb))) {}


  template<class SP1T, class SP2T>
  void apply(const VectorField<SP1T>& x, VectorField<SP2T>& y) const {

    for(F i=F::U; i<GridT::sdim; ++i) {
      sop_->apply(x(i), y(i));
    }

  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return sop_->grid();
  };
  /// \note dirty
  Teuchos::RCP<const GridT> gridC() const {
    return sop_->gridC();
  };
  /// \note dirty
  Teuchos::RCP<const GridT> gridF() const {
    return sop_->gridF();
  };

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    sop_->setParameter(para);
  }

  void print(std::ostream& out=std::cout) const {
    sop_->print();
  }

  const std::string getLabel() const {
    return "VectorFieldOpWrap("+sop_->getLabel()+") ";
  };

}; // end of class VectorFieldOpWrap



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_VECTORFIELDOPWRAP_HPP
