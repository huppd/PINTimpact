#pragma once
#ifndef PIMPACT_TIMESTOKESLSMOOTHER_HPP
#define PIMPACT_TIMESTOKESLSMOOTHER_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_TimeField.hpp"


namespace Pimpact {

extern "C" {

  void OP_TimeStokesLSmoother(
    const int& dimens,
    const int* const N,
    const int* const bl,
    const int* const bu,
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
    const double* const vel,
    const double* const p,
    double* const r_vel,
    double* const r_p,
    const int& L );

}

/// \ingroup CompoundOperator
///
// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class OperatorT>
class TimeStokesLSmoother {

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
  TimeStokesLSmoother( const Teuchos::RCP<const OperatorT>& op, Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList() ):
    op_( op ),
    numIters_( pl->get<int>("numIters",4) ) {};

  void apply(const DomainFieldT& x, RangeFieldT& y, const Ordinal L=1, const Belos::ETrans
      trans=Belos::NOTRANS  ) const {

    Scalar pi = 4.*std::atan(1.);
    Scalar idt = ((Scalar)space()->nGlo()[3])/2./pi;
    Scalar re = space()->getDomainSize()->getRe();
    Scalar mulI = space()->getDomainSize()->getAlpha2()*idt/re;

    auto& xu = x.getVField();
    auto& xp = x.getSField();
    auto& yu = y.getVField();
    auto& yp = y.getSField();

    for( int iters=0; iters<numIters_; ++iters ) {


      xu.exchange();
      //xp->exchange();

      for( Ordinal i=space()->si(F::S,3); i<space()->ei(F::S,3); ++i ) {
        xu(i-1).exchange();
        xu(i).exchange();
        xp(i).exchange();
      }

      int dimens = SpaceT::sdim;

      OP_TimeStokesLSmoother(
        dimens,
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->getBCLocal()->getBCL(),
        space()->getBCLocal()->getBCU(),
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
        xu.getConstRawPtr(),
        xp.getConstRawPtr(),
        yu.getRawPtr(),
        yp.getRawPtr(),
        L );

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


}; // end of class TimeStokesLSmoother



} // end of namespace Pimpact



#ifdef COMPILE_ETI
#include "Pimpact_TimeStokesOp.hpp"
extern template class Pimpact::TimeStokesLSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::TimeStokesLSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,4> > >;
#endif

#endif // end of #ifndef PIMPACT_TIMESTOKESBSMOOTHER_HPP
