/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_TIMENSOp_HPP
#define PIMPACT_TIMENSOp_HPP


#include "Teuchos_RCP.hpp"
#include "Pimpact_ScalarField.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_DiffusionOp.hpp"
#include "Pimpact_InterpolateS2VOp.hpp"
#include "Pimpact_TimeField.hpp"




namespace Pimpact {

extern "C" {

  void OP_TimeNS(
    const int& dimens,
    const int* const N,
    const int* const bl, const int* const bu,
    const int* const cL, const int* const cU,
    const int* const dl, const int* const du,
    const int* const gl, const int* const gu,
    const int* const ss, const int* const nn,
    const int* const su, const int* const nu,
    const int* const sv, const int* const nv,
    const int* const sw, const int* const nw,
    const double* const c1uD, const double* const c2vD, const double* const c3wD,
    const double* const c1uU, const double* const c2vU,	const double* const c3wU,
    const double* const c1pD,	const double* const c2pD,	const double* const c3pD,
    const double* const c1pU,	const double* const c2pU,	const double* const c3pU,
    const double* const c11p,	const double* const c22p,	const double* const c33p,
    const double* const c11u,	const double* const c22v,	const double* const c33w,
    const double* const cD1,
    const double* const cD2,
    const double* const cD3,
    const double* const cG1,
    const double* const cG2,
    const double* const cG3,
    const double& mulI,
    const double& mulL,
    const double* const windU,
    const double* const windV,
    const double* const windW,
    const double* const veln,
    const double* const pn,
    double* const r_vel,
    double* const r_p);
}


/// \ingroup CompoundOperator
///
// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class ST>
class TimeNSOp {

public:

  using GridT = ST;

  using DomainFieldT = CompoundField<TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >;
  using RangeFieldT = CompoundField<TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >;

protected:

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  static const int sdim = GridT::sdim;
  static const int dimension = GridT::dimension;

  static const int dimNC = GridT::dimNC;

  Teuchos::RCP<TimeField<VectorField<ST>>> windU_;
  Teuchos::RCP<TimeField<VectorField<ST>>> windV_;
  Teuchos::RCP<TimeField<VectorField<ST>>> windW_;

  Teuchos::RCP<const InterpolateS2V<ST>> interpolateS2V_;
  Teuchos::RCP<const InterpolateV2S<Scalar, Ordinal, sdim, dimension, dimNC>> interpolateV2S_;
  Teuchos::RCP<const ConvectionSOp<ST>> conv_;
  Teuchos::RCP<const DiffusionOp<ST>> helm_;
  Teuchos::RCP<const GradOp<ST>> grad_;
  Teuchos::RCP<const DivOp<ST>>  div_;

public:

  TimeNSOp(const Teuchos::RCP<const GridT>& grid):
    windU_(create<TimeField<VectorField<ST> > >(grid)),
    windV_(create<TimeField<VectorField<ST> > >(grid)),
    windW_(create<TimeField<VectorField<ST> > >(grid)),
    interpolateS2V_(create<InterpolateS2V>(grid)),
    interpolateV2S_(createInterpolateV2S(grid)),
    conv_(create<ConvectionSOp<ST> >(grid)),
    helm_(create<DiffusionOp<ST>   >(grid)),
    grad_(create<GradOp<ST>        >(grid)),
    div_ (create<DivOp<ST>         >(grid)) {};


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    Scalar pi = 4.*std::atan(1.);
    Scalar idt = ((Scalar)grid()->nGlo()[3])/2./pi;
    Scalar re = grid()->getDomainSize()->getRe();
    Scalar mulI = grid()->getDomainSize()->getAlpha2()*idt/re;

    auto& xu = x.getVField();
    auto& xp = x.getSField();
    auto& yu = y.getVField();
    auto& yp = y.getSField();

    xu.exchange();

    for(Ordinal i=grid()->si(F::S, 3)-1; i<grid()->ei(F::S, 3); ++i) {
      xu(i).exchange();
      xp(i).exchange();
    }

    int dimens = GridT::sdim;

    OP_TimeNS(
      dimens,
      grid()->nLoc(),
      grid()->bl(),
      grid()->bu(),
      grid()->nl(),
      grid()->nu(),
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
      conv_->getCD(X, F::U),
      conv_->getCD(Y, F::V),
      conv_->getCD(Z, F::W),
      conv_->getCU(X, F::U),
      conv_->getCU(Y, F::V),
      conv_->getCU(Z, F::W),
      conv_->getCD(X, F::S),
      conv_->getCD(Y, F::S),
      conv_->getCD(Z, F::S),
      conv_->getCU(X, F::S),
      conv_->getCU(Y, F::S),
      conv_->getCU(Z, F::S),
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
      windU_->getConstRawPtr(),
      windV_->getConstRawPtr(),
      windW_->getConstRawPtr(),
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


  void computeResidual(const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res) const {
    apply(x, res);
    res.add(1., b, -1., res);
  }

  void assignField(const DomainFieldT& cmv) {

    Ordinal nt = grid()->nLoc(3) + grid()->bu(3) - grid()->bl(3);

    auto& mv = cmv.getVField();

    ScalarField<ST> temp(grid());

    mv.exchange();

    for(Ordinal it=0; it<nt; ++it) {

      interpolateV2S_->apply(mv(it)(F::U), temp);
      for(F j=F::U; j<GridT::sdim; ++j) {
        interpolateS2V_->apply(temp, windU_->operator()(it)(j));
      }
      interpolateV2S_->apply(mv(it)(F::V), temp);
      for(F j=F::U; j<GridT::sdim; ++j) {
        interpolateS2V_->apply(temp, windV_->operator()(it)(j));
      }
      if(3==GridT::sdim) {
        interpolateV2S_->apply(mv(it)(F::W), temp);
        for(F j=F::U; j<GridT::sdim; ++j) {
          interpolateS2V_->apply(temp, windW_->operator()(it)(j));
        }
      }
    }
    windU_->changed();
    windV_->changed();
    if(3==GridT::sdim)
      windW_->changed();

  };

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return conv_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  bool hasApplyTranspose() const {
    return false;
  }

  Teuchos::RCP<const DiffusionOp<ST> > getDiffusionOp() const {
    return helm_;
  }
  Teuchos::RCP<const GradOp<ST> > getGradOp() const {
    return grad_;
  }
  Teuchos::RCP<const DivOp<ST> > getDivOp() const {
    return div_;
  }
  Teuchos::RCP<const ConvectionSOp<ST> > getConvOp() const {
    return conv_;
  }

  Teuchos::RCP<TimeField<VectorField<ST> > > getWindU_() const {
    return windU_;
  }
  Teuchos::RCP<TimeField<VectorField<ST> > > getWindV_() const {
    return windV_;
  }
  Teuchos::RCP<TimeField<VectorField<ST> > > getWindW_() const {
    return windW_;
  }



  void print(std::ostream& out=std::cout) const {
    out << getLabel() << ":\n";
  }

  const std::string getLabel() const {
    return "TimeNSOp";
  };

}; // end of class TimeNSOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMENSOp_HPP
