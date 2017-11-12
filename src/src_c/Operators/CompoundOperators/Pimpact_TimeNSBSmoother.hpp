#pragma once
#ifndef PIMPACT_TIMENSBSMOOTHER_HPP
#define PIMPACT_TIMENSBSMOOTHER_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_TimeField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"



namespace Pimpact {

extern "C" {

  void OP_TimeNSBSmoother(
    const int& dimens,
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
    const int&    direction_flag );

}

/// \ingroup CompoundOperator
///
// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class OperatorT>
class TimeNSBSmoother {

public:

  using SpaceT = typename OperatorT::SpaceT;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using DomainFieldT = CompoundField< TimeField<VectorField<SpaceT> >, TimeField<ScalarField<SpaceT> > >;
  using RangeFieldT = CompoundField< TimeField<VectorField<SpaceT> >, TimeField<ScalarField<SpaceT> > >;

protected:

  const Teuchos::RCP<const OperatorT> op_;
  int numIters_;

public:

  /// \note todo constructor from space
  TimeNSBSmoother( const Teuchos::RCP<const OperatorT>& op , Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList()):
    op_( op ),
    numIters_( pl->get<int>("numIters",4) )	{};

  void apply(const DomainFieldT& x, RangeFieldT& y ) const {

    Scalar pi = 4.*std::atan(1.);
    Scalar idt = ((Scalar)space()->nGlo()[3])/2./pi;
    Scalar re = space()->getDomainSize()->getRe();
    Scalar mulI = space()->getDomainSize()->getAlpha2()*idt/re;

    int direction_flag = 0;

    auto& xu = x.getVField();
    auto& xp = x.getSField();
    auto& yu = y.getVField();
    auto& yp = y.getSField();

    for( int iters=0; iters<numIters_; ++iters ) {

      direction_flag++;

      xu.exchange();
      //xp.exchange();

      for( Ordinal i=space()->si(F::S,3); i<space()->ei(F::S,3); ++i ) {
        xu(i-1).exchange();
        xu(i).exchange();
        xp(i).exchange();
      }

      int dimens = SpaceT::sdim;
 
      OP_TimeNSBSmoother(
          dimens,
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->getBCLocal()->getBCL(),
          space()->getBCLocal()->getBCU(),
          space()->nl(),
          space()->nu(),
          space()->dl(),
          space()->du(),
          space()->gl(),
          space()->gu(),
          space()->sInd(F::S),
          space()->eInd(F::S),
          space()->sInd(F::U),
          space()->eInd(F::U),
          space()->sInd(F::V),
          space()->eInd(F::V),
          space()->sInd(F::W),
          space()->eInd(F::W),
          op_->getConvOp()->getCD(X,F::U),
          op_->getConvOp()->getCD(Y,F::V),
          op_->getConvOp()->getCD(Z,F::W),
          op_->getConvOp()->getCU(X,F::U),
          op_->getConvOp()->getCU(Y,F::V),
          op_->getConvOp()->getCU(Z,F::W),
          op_->getConvOp()->getCD(X,F::S),
          op_->getConvOp()->getCD(Y,F::S),
          op_->getConvOp()->getCD(Z,F::S),
          op_->getConvOp()->getCU(X,F::S),
          op_->getConvOp()->getCU(Y,F::S),
          op_->getConvOp()->getCU(Z,F::S),
          op_->getHelmholtzOp()->getC(X,F::S),
          op_->getHelmholtzOp()->getC(Y,F::S),
          op_->getHelmholtzOp()->getC(Z,F::S),
          op_->getHelmholtzOp()->getC(X,F::U),
          op_->getHelmholtzOp()->getC(Y,F::V),
          op_->getHelmholtzOp()->getC(Z,F::W),
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
          direction_flag );

      for( Ordinal i=space()->si(F::S,3); i<space()->ei(F::S,3); ++i ) {
        yu(i).changed();
        yp(i).changed();
      }

      yu.changed();
      yp.changed();
    }
  }

  void assignField( const DomainFieldT& mv ) { };

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op_->space();
  };

  void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "TimeNSBSmoother ";
  };

}; // end of class TimeNSBSmoother


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMENSBSMOOTHER_HPP
