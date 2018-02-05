#pragma once
#ifndef PIMPACT_DIVGRADOP_HPP
#define PIMPACT_DIVGRADOP_HPP


#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \brief "laplace" for pressure.
/// \ingroup BaseOperator
/// \note todo not workin properly?
/// \warning does not hold test.
template<class SpT>
class DivGradOp {

public:

  using GridT = SpT;

  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  Teuchos::RCP<DivOp<GridT> > div_;
  Teuchos::RCP<GradOp<GridT> > grad_;

public:

  DivGradOp(const Teuchos::RCP<const GridT>& grid):
    div_ (create<DivOp>(grid)),
    grad_(create<GradOp>(grid)) {};

  DivGradOp(
    const Teuchos::RCP<DivOp<GridT> >& div,
    const Teuchos::RCP<GradOp<GridT> >& grad):
    div_ (div),
    grad_(grad) {};

  void apply(const DomainFieldT& x, RangeFieldT& y,
      const Belos::ETrans trans=Belos::NOTRANS, const Add add=Add::N) const {

    VectorField<GridT> temp(grid());

    switch(trans) {
    case Belos::NOTRANS : {
      grad_->applyG(x, temp, add);
      grad_->applyJ(temp);
      div_->apply(temp, y, add);
      break;
    }
    case Belos::TRANS : {
      div_->apply(x, temp, add);
      grad_->apply(temp, y, add);
      break;
    }
    case Belos::CONJTRANS : {
      div_->apply(x, temp, add);
      grad_->apply(temp, y, add);
      break;
    }
    }
  }


  void computeResidual(const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res) const {
    apply(x, res);
    res.add(1., b, -1., res);
  }


  /// \todo fix for periodic BC
  void applyInvDiag(const DomainFieldT& x, RangeFieldT& y) const {


    for(OT k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
      for(OT j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
        for(OT i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {

          ST diag = ((3==GridT::sdim)?innerDiag3D(i, j, k):innerDiag2D(i, j, k));

          assert(0!=diag);
          y(i, j, k) = x(i, j, k)/std::fabs(diag);
        }

    y.changed();
  }

  Teuchos::RCP<const DivOp<GridT> > getDivOp() const {
    return div_;
  }
  Teuchos::RCP<const GradOp<GridT> > getGradOp() const {
    return grad_;
  }

  void assignField(const DomainFieldT& mv) {};

  bool hasApplyTranspose() const {
    return true;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return div_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  constexpr const std::string getLabel() const {
    return "DivGrad";
  };

  void print(std::ostream& out=std::cout) const {
    out <<"---" <<getLabel() <<"---\n";
    div_->print(out);
    grad_->print(out);
    grid()->getInterpolateV2S()->print();
  }


  ST innerDiag3D(const OT i, const OT j, const OT k) const {

    ST diag = 0.;
    //const ST eps = 0.1;
    const ST eps = 1./static_cast<ST>(GradOp<GridT>::epsI);

    const bool bcX = (grid()->getBCLocal()->getBCL(X) > 0 && i==grid()->si(F::S, X)) ||
                     (              grid()->getBCLocal()->getBCU(X) > 0 && i==grid()->ei(F::S, X)) ;
    const bool bcY = (grid()->getBCLocal()->getBCL(Y) > 0 && j==grid()->si(F::S, Y)) ||
                     (              grid()->getBCLocal()->getBCU(Y) > 0 && j==grid()->ei(F::S, Y)) ;
    const bool bcZ = (grid()->getBCLocal()->getBCL(Z) > 0 && k==grid()->si(F::S, Z)) ||
                     (              grid()->getBCLocal()->getBCU(Z) > 0 && k==grid()->ei(F::S, Z)) ;

    const ST epsX = ((bcY||bcZ)?eps:1.);
    const ST epsY = ((bcX||bcZ)?eps:1.);
    const ST epsZ = ((bcX||bcY)?eps:1.);

    // X direction
    for(OT ii=grid()->dl(X); ii<=grid()->du(X); ++ii) {
      if(0<grid()->getBCLocal()->getBCL(X) && i+ii==grid()->si(F::U, X, B::Y)) {
        for(OT iii=0; iii<=grid()->du(X); ++iii)
          diag -= div_->getC(X, i, ii) * epsX * grad_->getC(X, 1+iii, -iii-ii-1)
                  * grid()->getInterpolateV2S()->getC(X, 1, iii) /
                  grid()->getInterpolateV2S()->getC(X, 1, -1);
      } else if(0<grid()->getBCLocal()->getBCU(X) && i+ii==grid()->ei(F::U, X, B::Y)) {
        for(OT iii=grid()->dl(X); iii<=-1; ++iii)
          diag -= div_->getC(X, i, ii) * epsX * grad_->getC(X, grid()->ei(F::U, X, B::Y)+iii, -iii-ii)
                  * grid()->getInterpolateV2S()->getC(X, grid()->ei(F::U, X, B::Y), iii) /
                  grid()->getInterpolateV2S()->getC(X, grid()->ei(F::U, X, B::Y), 0);
      } else if(i+ii>=0 && i+ii<=grid()->nLoc(X))
        diag += div_->getC(X, i, ii) * epsX * grad_->getC(X, i+ii, -ii);
    }

    // Y direction
    for(OT jj=grid()->dl(Y); jj<=grid()->du(Y); ++jj) {
      if(0<grid()->getBCLocal()->getBCL(Y) && j+jj==grid()->si(F::V, Y, B::Y)) {
        for(OT jjj=0; jjj<=grid()->du(Y); ++jjj)
          diag -= div_->getC(Y, j, jj) * epsY * grad_->getC(Y, 1+jjj, -jjj+j-1)
                  * grid()->getInterpolateV2S()->getC(Y, 1, jjj) /
                  grid()->getInterpolateV2S()->getC(Y, 1, -1);
      } else if(0<grid()->getBCLocal()->getBCU(Y) && j+jj==grid()->ei(F::V, Y, B::Y)) {
        for(OT jjj=grid()->dl(Y); jjj<=-1; ++jjj)
          diag -= div_->getC(Y, j, jj) * epsY * grad_->getC(Y, grid()->ei(F::V, Y, B::Y)+jjj, -jjj-jj)
                  * grid()->getInterpolateV2S()->getC(Y, grid()->ei(F::V, Y, B::Y), jjj) /
                  grid()->getInterpolateV2S()->getC(Y, grid()->ei(F::V, Y, B::Y), 0);
      } else if(j+jj>=0 && j+jj<=grid()->nLoc(Y))
        diag += div_->getC(Y, j, jj)*epsY*grad_->getC(Y, j+jj, -jj);
    }

    // Z direction
    for(OT kk=grid()->dl(Z); kk<=grid()->du(Z); ++kk) {
      if(0<grid()->getBCLocal()->getBCL(Z) && k+kk==grid()->si(F::W, Z, B::Y)) {
        for(OT kkk=0; kkk<=grid()->du(Z); ++kkk)
          diag -= div_->getC(Z, k, kk) * epsZ * grad_->getC(Z, 1+kkk, -kkk+k-1)
                  * grid()->getInterpolateV2S()->getC(Z, 1, kkk) /
                  grid()->getInterpolateV2S()->getC(Z, 1, -1);
      } else if(0<grid()->getBCLocal()->getBCU(Z) && k+kk==grid()->ei(F::W, Z, B::Y)) {
        for(OT kkk=grid()->dl(Z); kkk<=-1; ++kkk)
          diag -= div_->getC(Z, k, kk) * epsZ * grad_->getC(Z, grid()->ei(F::W, Z, B::Y)+kkk, -kkk-kk)
                  * grid()->getInterpolateV2S()->getC(Z, grid()->ei(F::W, Z, B::Y), kkk) /
                  grid()->getInterpolateV2S()->getC(Z, grid()->ei(F::W, Z, B::Y), 0);
      } else if(k+kk>=0 && k+kk<=grid()->nLoc(Z))
        diag += div_->getC(Z, k, kk) * epsZ * grad_->getC(Z, k+kk, -kk);
    }

    return diag;
  }

  ST innerDiag2D(const OT i, const OT j, const OT k) const {

    ST diag = 0.;
    //const ST eps = 0.1;
    const ST eps = 1./static_cast<ST>(GradOp<GridT>::epsI);

    const bool bcX = (grid()->getBCLocal()->getBCL(X) > 0 && i==grid()->si(F::S, X)) ||
                     (              grid()->getBCLocal()->getBCU(X) > 0 && i==grid()->ei(F::S, X)) ;
    const bool bcY = (grid()->getBCLocal()->getBCL(Y) > 0 && j==grid()->si(F::S, Y)) ||
                     (              grid()->getBCLocal()->getBCU(Y) > 0 && j==grid()->ei(F::S, Y)) ;

    const ST epsX = (bcY?eps:1.);
    const ST epsY = (bcX?eps:1.);

    // X direction
    for(OT ii=grid()->dl(X); ii<=grid()->du(X); ++ii) {
      if(0<grid()->getBCLocal()->getBCL(X) && i+ii==grid()->si(F::U, X, B::Y)) {
        for(OT iii=0; iii<=grid()->du(X); ++iii)
          diag -= div_->getC(X, i, ii) * epsX * grad_->getC(X, 1+iii, -iii-ii-1)
                  * grid()->getInterpolateV2S()->getC(X, 1, iii) /
                  grid()->getInterpolateV2S()->getC(X, 1, -1);
      } else if(0<grid()->getBCLocal()->getBCU(X) && i+ii==grid()->ei(F::U, X, B::Y)) {
        for(OT iii=grid()->dl(X); iii<=-1; ++iii)
          diag -= div_->getC(X, i, ii) * epsX * grad_->getC(X, grid()->ei(F::U, X, B::Y)+iii, -iii-ii)
                  * grid()->getInterpolateV2S()->getC(X, grid()->ei(F::U, X, B::Y), iii) /
                  grid()->getInterpolateV2S()->getC(X, grid()->ei(F::U, X, B::Y), 0);
      } else if(i+ii>=0 && i+ii<=grid()->nLoc(X))
        diag += div_->getC(X, i, ii) * epsX * grad_->getC(X, i+ii, -ii);
    }

    // Y direction
    for(OT jj=grid()->dl(Y); jj<=grid()->du(Y); ++jj) {
      if(0<grid()->getBCLocal()->getBCL(Y) && j+jj==grid()->si(F::V, Y, B::Y)) {
        for(OT jjj=0; jjj<=grid()->du(Y); ++jjj)
          diag -= div_->getC(Y, j, jj) * epsY * grad_->getC(Y, 1+jjj, -jjj+j-1)
                  * grid()->getInterpolateV2S()->getC(Y, 1, jjj) /
                  grid()->getInterpolateV2S()->getC(Y, 1, -1);
      } else if(0<grid()->getBCLocal()->getBCU(Y) && j+jj==grid()->ei(F::V, Y, B::Y)) {
        for(OT jjj=grid()->dl(Y); jjj<=-1; ++jjj)
          diag -= div_->getC(Y, j, jj) * epsY * grad_->getC(Y, grid()->ei(F::V, Y, B::Y)+jjj, -jjj-jj)
                  * grid()->getInterpolateV2S()->getC(Y, grid()->ei(F::V, Y, B::Y), jjj) /
                  grid()->getInterpolateV2S()->getC(Y, grid()->ei(F::V, Y, B::Y), 0);
      } else if(j+jj>=0 && j+jj<=grid()->nLoc(Y))
        diag += div_->getC(Y, j, jj)*epsY*grad_->getC(Y, j+jj, -jj);
    }

    return diag;
  }

}; // end of DivGradOp



/// \relates DivGradOp
template<class GridT>
Teuchos::RCP<DivGradOp<GridT> > createDivGradOp(
    const Teuchos::RCP<DivOp<GridT> >& div,
    const Teuchos::RCP<GradOp<GridT> >& grad) {

  return Teuchos::rcp(new DivGradOp<GridT>(div, grad));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DIVGRADOP_HPP
