/// Pimpact 
/// \author Pietro Benedusi
/// \date 2016


#pragma once
#ifndef PIMPACT_TIMENS4DBSMOOTHER_HPP
#define PIMPACT_TIMENS4DBSMOOTHER_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_TimeField.hpp"



namespace Pimpact {

extern "C" void OP_TimeNS4DBSmoother(
  const int* const N,
  const int* const bl,
  const int* const bu,
  const int* const cl,
  const int* const cu,
  const int* const BCL,
  const int* const BCU,
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
  const double* const c1uD,
  const double* const c2vD,
  const double* const c3wD,
  const double* const c1uU,
  const double* const c2vU,
  const double* const c3wU,
  const double* const c1pD,
  const double* const c2pD,
  const double* const c3pD,
  const double* const c1pU,
  const double* const c2pU,
  const double* const c3pU,
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
  const double* const windU,
  const double* const windV,
  const double* const windW,
  const double* const rhs_vel,
  const double* const rhs_p,
  double* const vel,
  double* const p,
  const int&    direction_flag);


/// \ingroup CompoundOperator
///
// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class OperatorT>
class TimeNS4DBSmoother {

public:

  using GridT = typename OperatorT::GridT;

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  using DomainFieldT = CompoundField<TimeField<VectorField<GridT>>,
        TimeField<ScalarField<GridT>>>;

  using RangeFieldT = CompoundField<TimeField<VectorField<GridT>>,
        TimeField<ScalarField<GridT>>>;


protected:

  const Teuchos::RCP<const OperatorT> op_;
  int numIters_;

public:

  /// \note todo constructor from grid
  TimeNS4DBSmoother(const Teuchos::RCP<const OperatorT>& op, Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList()):
    op_(op),
    numIters_(pl->get<int>("numIters", 4))	{};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    Scalar pi = 4.*std::atan(1.);
    Scalar idt = ((Scalar)grid()->nGlo()[3])/2./pi;
    Scalar re = grid()->getDomainSize()->getRe();
    Scalar mulI = grid()->getDomainSize()->getAlpha2()*idt/re;

    int direction_flag = 0;

    for(int iters=0; iters<numIters_; ++iters) {

      // this is for alternating directions
      direction_flag++;

      auto& xu = x.getVField();
      auto& xp = x.getSField();
      auto& yu = y.getVField();
      auto& yp = y.getSField();

      xu.exchange();
      //		xp->exchange();

      for(Ordinal i=grid()->si(F::S, 3); i<grid()->ei(F::S, 3); ++i) {
        xu(i-1).exchange();
        xu(i).exchange();
        xp(i).exchange();
      }


      OP_TimeNS4DBSmoother(
        grid()->nLoc(),
        grid()->bl(),
        grid()->bu(),
        grid()->getBCLocal()->getBCL(),
        grid()->getBCLocal()->getBCU(),
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
        op_->getConvOp()->getCD(X, F::U),
        op_->getConvOp()->getCD(Y, F::V),
        op_->getConvOp()->getCD(Z, F::W),
        op_->getConvOp()->getCU(X, F::U),
        op_->getConvOp()->getCU(Y, F::V),
        op_->getConvOp()->getCU(Z, F::W),
        op_->getConvOp()->getCD(X, F::S),
        op_->getConvOp()->getCD(Y, F::S),
        op_->getConvOp()->getCD(Z, F::S),
        op_->getConvOp()->getCU(X, F::S),
        op_->getConvOp()->getCU(Y, F::S),
        op_->getConvOp()->getCU(Z, F::S),
        op_->getDiffusionOp()->getC(X, F::S),
        op_->getDiffusionOp()->getC(Y, F::S),
        op_->getDiffusionOp()->getC(Z, F::S),
        op_->getDiffusionOp()->getC(X, F::U),
        op_->getDiffusionOp()->getC(Y, F::V),
        op_->getDiffusionOp()->getC(Z, F::W),
        op_->getDivOp()->getC(X),
        op_->getDivOp()->getC(Y),
        op_->getDivOp()->getC(Z),
        op_->getGradOp()->getC(X),
        op_->getGradOp()->getC(Y),
        op_->getGradOp()->getC(Z),
        mulI,
        1./re,
        op_->getWindU_()->getConstRawPtr(),
        op_->getWindV_()->getConstRawPtr(),
        op_->getWindW_()->getConstRawPtr(),
        xu.getConstRawPtr(),
        xp.getConstRawPtr(),
        yu.getRawPtr(),
        yp.getRawPtr(),
        direction_flag);

      for(Ordinal i=grid()->si(F::S, 3); i<grid()->ei(F::S, 3); ++i) {
        yu(i).changed();
        yp(i).changed();
      }

      yu.changed();
      yp.changed();
    }
  }

  void assignField(const DomainFieldT& mv) { };

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "TimeNS4DBSmoother ";
  };

}; // end of class TimeNS4DBSmoother


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMENS4DBSMOOTHER_HPP
