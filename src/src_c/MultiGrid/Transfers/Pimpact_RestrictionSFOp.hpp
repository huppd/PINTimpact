#pragma once
#ifndef PIMPACT_RESTRICTIONSFOP_HPP
#define PIMPACT_RESTRICTIONSFOP_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Pimpact_RestrictionBaseOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Grid.hpp"
#include "Pimpact_Stencil.hpp"




namespace Pimpact {


/// \brief Opetartor that restricts from a fine grid to a coarse grid
/// \todo c++fy
/// \tparam GT type of the \c Grid
template<class GT>
class RestrictionSFOp : private RestrictionBaseOp<GT> {

  static const int sdim = GT::sdim;
  static const int dimension = GT::dimension;

  using Scalar = typename GT::Scalar;
  using Ordinal = typename GT::Ordinal;

public:

  using GridT = GT;

  using FGridT = GridT;
  using CGridT = GridT;

  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

  using Stenc = Stencil<Scalar, Ordinal, 1, -1, 1 >;

protected:

  Teuchos::Tuple<Stenc, 3> cRS_;

  void initSF() {

    for(int i=0; i<3; ++i) {

      cRS_[i] = Stenc(this->iimax_[i]);

      MG_getCRS(
        this->iimax_[i],
        (this->nGather_[i]>1)?
        gridF()->getBCLocal()->getBCL(i):
        gridC()->getBCLocal()->getBCL(i),
        (this->nGather_[i]>1)?
        gridF()->getBCLocal()->getBCU(i):
        gridC()->getBCLocal()->getBCU(i),
        this->dd_[i],
        gridF()->getGridSizeLocal()->get(i),
        gridF()->bl(i),
        gridF()->bu(i),
        gridF()->getCoordinatesLocal()->getX(F::S, i),
        cRS_[i].get());
    }
  }

public:

  RestrictionSFOp(
    const Teuchos::RCP<const GridT>& gridF,
    const Teuchos::RCP<const GridT>& gridC):
    RestrictionBaseOp<GT>(gridF, gridC) {

    initSF();
  }


  RestrictionSFOp(
    const Teuchos::RCP<const GridT>& gridF,
    const Teuchos::RCP<const GridT>& gridC,
    const Teuchos::Tuple<int, dimension>& np):
    RestrictionBaseOp<GT>(gridF, gridC, np) {

    initSF();
  }


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    assert(x.getType()==y.getType());
    assert(x.getType()==F::S);

    x.exchange();

    const int sdimens = sdim; // so ugly
    MG_restrictFW(
      sdimens,
      gridF()->nLoc(),
      gridF()->bl(),
      gridF()->bu(),
      gridC()->nLoc(),
      gridC()->bl(),
      gridC()->bu(),
      this->iimax_.getRawPtr(),
      this->dd_.getRawPtr(),
      cRS_[0].get(),
      cRS_[1].get(),
      cRS_[2].get(),
      x.getConstRawPtr(),
      y.getRawPtr());

    this->gather(y.getRawPtr());

    y.changed();
  }


  void print(std::ostream& out=std::cout) const {

    out << "=== Restriction OP ===\n";
    out << "nGather:\t" << this->nGather_ << "\n";
    out << "rankc2:\t" << this->rankc2_ << "\n";
    out << "comm2:\t" << this->comm2_ << "\n";

    out << " --- scalar stencil: ---";
    for(int j=0; j<3; ++j) {
      out << "\ndir: " << j << "\n";
      cRS_[j].print(out);
    }
  }


  Teuchos::Tuple<Ordinal, dimension> getDD() const {
    return this->dd_;
  };

  Teuchos::RCP<const GridT> gridC() const {
    return this->gridC_;
  };
  Teuchos::RCP<const GridT> gridF() const {
    return this->gridF_;
  };

  const std::string getLabel() const {
    return "Restriction SF";
  };

}; // end of class RestrictionSFOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_RESTRICTIONSFOP_HPP
