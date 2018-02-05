#pragma once
#ifndef PIMPACT_TIMESTOKESOP_HPP
#define PIMPACT_TIMESTOKESOP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_TimeField.hpp"




namespace Pimpact {

extern "C" {

  void OP_TimeStokes(
    const int& dimens,
    const int* const N,
    const int* const bl,
    const int* const bu,
    const int* const dl,
    const int* const du,
    const int* const gl,
    const int* const gu,
    const int* const ss,
    const int* const nn,
    const int* const su,
    const int* const nu,
    const int* const sv,
    const int* const nv,
    const int* const sw,
    const int* const nw,
    const double* const c11p,
    const double* const c22p,
    const double* const c33p,
    const double* const c11u,
    const double* const c22v,
    const double* const c33w,
    const double* const cD1,
    const double* const cD2,
    const double* const cD3,
    const double* const cG1,
    const double* const cG2,
    const double* const cG3,
    const double& mulI,
    const double& mulL,
    //      const double* const velp,
    const double* const veln,
    const double* const pn,
    double* const r_vel,
    double* const r_p);
}



/// \ingroup CompoundOperator
///
// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class ST>
class TimeStokesOp {

public:

  using GridT = ST;

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  using DomainFieldT = CompoundField<TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >;
  using RangeFieldT = CompoundField<TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > > ;


protected:

  Teuchos::RCP<const HelmholtzOp<ST> > helm_;
  Teuchos::RCP<const GradOp<ST> > grad_;
  Teuchos::RCP<const DivOp<ST> > div_;

public:

  TimeStokesOp(const Teuchos::RCP<const GridT>& grid):
    helm_(create<HelmholtzOp<ST> >(grid)),
    grad_(create<GradOp<ST> >(grid)),
    div_(create<DivOp<ST> >(grid)) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    Scalar pi = 4.*std::atan(1.);
    Scalar idt = ((Scalar)grid()->nGlo(3))/2./pi;
    Scalar re = grid()->getDomainSize()->getRe();
    Scalar mulI = grid()->getDomainSize()->getAlpha2()*idt/re;

    auto& xu = x.getVField();
    auto& xp = x.getSField();
    auto& yu = y.getVField();
    auto& yp = y.getSField();

    xu.exchange();

    for(Ordinal i=grid()->si(F::S, 3)-1; i<grid()->ei(F::S, 3); ++i) {
      //xu(i-1).exchange();
      xu(i).exchange();
      xp(i).exchange();
    }

    int dimens = GridT::sdim;

    OP_TimeStokes(
      dimens,
      grid()->nLoc(),
      grid()->bl(),
      grid()->bu(),
      grid()->dl(),
      grid()->du(),
      grid()->gl(),
      grid()->gu(),
      grid()->sInd(F::S),
      grid()->eInd(F::S),
      grid()->sInd(F::U),
      grid()->eInd(F::U),
      grid()->sInd(F::V),
      grid()->eInd(F::V),
      grid()->sInd(F::W),
      grid()->eInd(F::W),
      helm_->getC(X, F::S),
      helm_->getC(Y, F::S),
      helm_->getC(Z, F::S),
      helm_->getC(X, F::U),
      helm_->getC(Y, F::V),
      helm_->getC(Z, F::W),
      div_->getC(X),
      div_->getC(Y),
      div_->getC(Z),
      grad_->getC(X),
      grad_->getC(Y),
      grad_->getC(Z),
      mulI,
      1./re,
      xu.getConstRawPtr(),
      xp.getConstRawPtr(),
      yu.getRawPtr(),
      yp.getRawPtr());

    for(Ordinal i=grid()->si(F::S, 3)-1; i<grid()->ei(F::S, 3); ++i) {
      yu(i).changed();
      yp(i).changed();
    }

    yu.changed();
    yp.changed();
  }

  void assignField(const DomainFieldT& mv) { };

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return helm_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  bool hasApplyTranspose() const {
    return false;
  }



  Teuchos::RCP<const HelmholtzOp<ST> > getHelmholtzOp() const {
    return helm_;
  }
  Teuchos::RCP<const GradOp<ST> > getGradOp() const {
    return grad_;
  }
  Teuchos::RCP<const DivOp<ST> > getDivOp() const {
    return div_;
  }

  const std::string getLabel() const {
    return "TimeStokesOp ";
  };

  void print(std::ostream& out=std::cout) const {
    out <<getLabel() <<":\n";
  }

}; // end of class TimeStokesOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMESTOKESOP_HPP
