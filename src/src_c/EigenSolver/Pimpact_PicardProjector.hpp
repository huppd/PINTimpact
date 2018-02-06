#pragma once
#ifndef PIMPACT_PICARDPROJECTOR_HPP
#define PIMPACT_PICARDPROJECTOR_HPP


#include "Pimpact_CompoundField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_DivGradNullSpace.hpp"


namespace Pimpact {



template<class OperatorT>
class PicardProjector {

  using GridT = typename OperatorT::GridT;

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using RangeFieldT = typename OperatorT::RangeFieldT;

  CompoundField<VectorField<GridT>, ScalarField<GridT>> nullspace_;
  CompoundField<VectorField<GridT>, ScalarField<GridT>> projection_;

  ST dotNP_;

  void setCornersZero(ScalarField<GridT>& rhs) const {

    auto grid = nullspace_.grid();
    F f = rhs.getType();

    // BC XY
    if(grid->getBCLocal()->getBCL(X)>0 && grid->getBCLocal()->getBCL(Y)>0) {
      OT i = grid->si(f, X, B::Y);
      OT j = grid->si(f, Y, B::Y);
      for(OT k=grid->si(f, Z, B::Y); k<=grid->ei(f, Z, B::Y); ++k)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCL(X)>0 && grid->getBCLocal()->getBCU(Y)>0) {
      OT i = grid->si(f, X, B::Y);
      OT j = grid->ei(f, Y, B::Y);
      for(OT k=grid->si(f, Z, B::Y); k<=grid->ei(f, Z, B::Y); ++k)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(X)>0 && grid->getBCLocal()->getBCL(Y)>0) {
      OT i = grid->ei(f, X, B::Y);
      OT j = grid->si(f, Y, B::Y);
      for(OT k=grid->si(f, Z, B::Y); k<=grid->ei(f, Z, B::Y); ++k)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(X)>0 && grid->getBCLocal()->getBCU(Y)>0) {
      OT i = grid->ei(f, X, B::Y);
      OT j = grid->ei(f, Y, B::Y);
      for(OT k=grid->si(f, Z, B::Y); k<=grid->ei(f, Z, B::Y); ++k)
        rhs(i, j, k) = 0.;
    }

    // BC XZ
    if(grid->getBCLocal()->getBCL(X)>0 && grid->getBCLocal()->getBCL(Z)>0) {
      OT i = grid->si(f, X, B::Y);
      OT k = grid->si(f, Z, B::Y);
      for(OT j=grid->si(f, Y, B::Y); j<=grid->ei(f, Y, B::Y); ++j)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCL(X)>0 && grid->getBCLocal()->getBCU(Z)>0) {
      OT i = grid->si(f, X, B::Y);
      OT k = grid->ei(f, Z, B::Y);
      for(OT j=grid->si(f, Y, B::Y); j<=grid->ei(f, Y, B::Y); ++j)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(X)>0 && grid->getBCLocal()->getBCL(Z)>0) {
      OT i = grid->ei(f, X, B::Y);
      OT k = grid->si(f, Z, B::Y);
      for(OT j=grid->si(f, Y, B::Y); j<=grid->ei(f, Y, B::Y); ++j)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(X)>0 && grid->getBCLocal()->getBCU(Z)>0) {
      OT i = grid->ei(f, X, B::Y);
      OT k = grid->ei(f, Z, B::Y);
      for(OT j=grid->si(f, Y, B::Y); j<=grid->ei(f, Y, B::Y); ++j)
        rhs(i, j, k) = 0.;
    }

    // BC YZ
    if(grid->getBCLocal()->getBCL(Y)>0 && grid->getBCLocal()->getBCL(Z)>0) {
      OT j = grid->si(f, Y, B::Y);
      OT k = grid->si(f, Z, B::Y);
      for(OT i=grid->si(f, X, B::Y); i<=grid->ei(f, X, B::Y); ++i)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCL(Y)>0 && grid->getBCLocal()->getBCU(Z)>0) {
      OT j = grid->si(f, Y, B::Y);
      OT k = grid->ei(f, Z, B::Y);
      for(OT i=grid->si(f, X, B::Y); i<=grid->ei(f, X, B::Y); ++i)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(Y)>0 && grid->getBCLocal()->getBCL(Z)>0) {
      OT j = grid->ei(f, Y, B::Y);
      OT k = grid->si(f, Z, B::Y);
      for(OT i=grid->si(f, X, B::Y); i<=grid->ei(f, X, B::Y); ++i)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(Y)>0 && grid->getBCLocal()->getBCU(Z)>0) {
      OT j = grid->ei(f, Y, B::Y);
      OT k = grid->ei(f, Z, B::Y);
      for(OT i=grid->si(f, X, B::Y); i<=grid->ei(f, X, B::Y); ++i)
        rhs(i, j, k) = 0.;
    }
    rhs.changed();
  }


  void project(VectorField<GridT>& rhs_v, ScalarField<GridT>& rhs_s) const {

    //setCornersZero(rhs_s);

    auto grid = nullspace_.grid();

    ST bla = -(nullspace_.getVField().dot(rhs_v) + nullspace_.getSField().dot(rhs_s)
       )/dotNP_;

    if(0==grid->rankST())
      std::cout << "Picard^-1"<< ": nullspace contributtion: " << std::abs(bla)  << "\n";

    if(std::abs(bla) >= Teuchos::ScalarTraits<ST>::eps()) {
      rhs_v.add(1., rhs_v, bla, projection_.getVField());
      rhs_s.add(1., rhs_s, bla, projection_.getSField());
    }

    setCornersZero(rhs_s);
  }

public:

  PicardProjector() {}

  PicardProjector(const Teuchos::RCP<const OperatorT>& op):
    nullspace_(op->grid()), projection_(op->grid()) {

    auto grid = nullspace_.grid();

    DivGradNullSpace<DivOp<GridT> > compNullspace;

    compNullspace.computeNullSpace(op->getOpV2S()->getOperatorPtr(),
        nullspace_.getSField(), true);

    // U inflow
    if(Pimpact::BC::Dirichlet==grid->bcl(Pimpact::X)) {
      OT si = grid()->si(F::U, X, B::Y);
      for(OT k=grid()->si(F::U, Z, B::Y); k<=grid()->ei(F::U, Z, B::Y); ++k)
        for(OT j=grid()->si(F::U, Y, B::Y); j<=grid()->ei(F::U, Y, B::Y); ++j)
          nullspace_.getVField()(F::U)(si, j, k) = -1.;
    }
    // U outflow
    if(Pimpact::BC::Dirichlet==grid->bcu(Pimpact::X)) {
      OT ei = grid()->ei(F::U, X, B::Y);
      for(OT k=grid()->si(F::U, Z, B::Y); k<=grid()->ei(F::U, Z, B::Y); ++k)
        for(OT j=grid()->si(F::U, Y, B::Y); j<=grid()->ei(F::U, Y, B::Y); ++j)
          nullspace_.getVField()(F::U)(ei, j, k) =  1.;
    }

    // V inflow
    if(Pimpact::BC::Dirichlet==grid->bcl(Pimpact::Y)) {
      OT sj = grid()->si(F::V, Y, B::Y);
      for(OT k=grid()->si(F::V, Z, B::Y); k<=grid()->ei(F::V, Z, B::Y); ++k)
        for(OT i=grid()->si(F::V, X, B::Y); i<=grid()->ei(F::V, X, B::Y); ++i)
          nullspace_.getVField()(F::V)(i, sj, k) = -1.;
    }
    // V outflow
    if(Pimpact::BC::Dirichlet==grid->bcu(Pimpact::Y)) {
      OT ej = grid()->ei(F::V, Y, B::Y);
      for(OT k=grid()->si(F::V, Z, B::Y); k<=grid()->ei(F::V, Z, B::Y); ++k)
        for(OT i=grid()->si(F::V, X, B::Y); i<=grid()->ei(F::V, X, B::Y); ++i)
          nullspace_.getVField()(F::V)(i, ej, k) = 1.;
    }

    // W in/outflow
    if(Pimpact::BC::Dirichlet==grid->bcl(Pimpact::Z)) {
      OT sk = grid()->si(F::W, Z, B::Y);
      for(OT j=grid()->si(F::W, Y, B::Y); j<=grid()->ei(F::W, Y, B::Y); ++j)
        for(OT i=grid()->si(F::W, X, B::Y); i<=grid()->ei(F::W, X, B::Y); ++i)
          nullspace_.getVField()(F::W)(i, j, sk) = -1.;
    }
    // W in/outflow
    if(Pimpact::BC::Dirichlet==grid->bcu(Pimpact::Z)) {
      OT ek = grid()->ei(F::W, Z, B::Y);
      for(OT j=grid()->si(F::W, Y, B::Y); j<=grid()->ei(F::W, Y, B::Y); ++j)
        for(OT i=grid()->si(F::W, X, B::Y); i<=grid()->ei(F::W, X, B::Y); ++i)
          nullspace_.getVField()(F::W)(i, j, ek) = 1.;
    }
    nullspace_.getVField().changed();

    setCornersZero(projection_.getSField());

    //ST blup =  1./nullspace_.norm();
    //nullspace_.scale(blup);
    //nullspace_.write(777);

    // outflow projection
    const ST pi = 4.*std::atan(1.);
    const ST width = 0.95;
    const ST eps = 0.;

    auto scalefunc = [=](ST x, ST y, ST z) ->ST {
      return (y<=width)?0.:(std::cos(pi*(y-width)/(1.-width) + pi)/2. + 0.5); };

    if(Pimpact::BC::Dirichlet==grid->bcu(Pimpact::Y)) {
      // outflow velocity
      OT ej = grid()->ei(F::V, Y, B::Y);
      for(OT k=grid()->si(F::V, Z, B::Y); k<=grid()->ei(F::V, Z, B::Y); ++k)
        for(OT i=grid()->si(F::V, X, B::Y); i<=grid()->ei(F::V, X, B::Y); ++i)
          projection_.getVField()(F::V)(i, ej, k) = 1.;

      //// pressure "outflow" 
      //OT ej = grid()->ei(F::S, Y, B::Y);
      //for(OT k=grid()->si(F::S, Z, B::Y); k<=grid()->ei(F::S, Z, B::Y); ++k)
        //for(OT i=grid()->si(F::S, X, B::Y); i<=grid()->ei(F::S, X, B::Y); ++i)
          //projection_.getSField()(i, ej, k) = 1.;
    }

    projection_.getSField().initFromFunction(scalefunc);
    //setCornersZero(projection_.getVField());
    setCornersZero(projection_.getSField());

    dotNP_ = nullspace_.dot(projection_);

    if(grid->rankST()==0) std::cout << "dotNP: " << dotNP_ << "\n";
    assert(dotNP_!=0.);
  }


  void operator()(RangeFieldT& rhs) const {

    auto grid = nullspace_.grid();

    if(0==grid->si(F::U, 3))
      project(rhs.getVField().get0Field(), rhs.getSField().get0Field());

    for(OT i=std::max(grid->si(F::U, 3), 1); i<=grid->ei(F::U, 3); ++i) {
      project(rhs.getVField().getCField(i), rhs.getSField().getCField(i));
      project(rhs.getVField().getSField(i), rhs.getSField().getSField(i));
    }
  }

}; // end of class EmptyProjector

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_PICARDPROJECTOR_HPP
