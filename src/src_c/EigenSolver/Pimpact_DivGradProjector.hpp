#pragma once
#ifndef PIMPACT_DGPROJECTOR_HPP
#define PIMPACT_DGPROJECTOR_HPP


#include "Pimpact_DivOp.hpp"
#include "Pimpact_DivGradNullSpace.hpp"


namespace Pimpact {



template<class OperatorT>
class DGProjector {

  using RangeFieldT = typename OperatorT::RangeFieldT;
  using GridT = typename RangeFieldT::GridT;
  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  //RangeFieldT nullspace_;

  void setCornersZero(RangeFieldT& rhs) const {

    auto grid = rhs.grid();

    // BC XY
    if(grid->getBCLocal()->getBCL(X)>0 && grid->getBCLocal()->getBCL(Y)>0) {
      OT i = grid->si(F::S, X, B::Y);
      OT j = grid->si(F::S, Y, B::Y);
      for(OT k=grid->si(F::S, Z, B::Y); k<=grid->ei(F::S, Z, B::Y); ++k)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCL(X)>0 && grid->getBCLocal()->getBCU(Y)>0) {
      OT i = grid->si(F::S, X, B::Y);
      OT j = grid->ei(F::S, Y, B::Y);
      for(OT k=grid->si(F::S, Z, B::Y); k<=grid->ei(F::S, Z, B::Y); ++k)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(X)>0 && grid->getBCLocal()->getBCL(Y)>0) {
      OT i = grid->ei(F::S, X, B::Y);
      OT j = grid->si(F::S, Y, B::Y);
      for(OT k=grid->si(F::S, Z, B::Y); k<=grid->ei(F::S, Z, B::Y); ++k)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(X)>0 && grid->getBCLocal()->getBCU(Y)>0) {
      OT i = grid->ei(F::S, X, B::Y);
      OT j = grid->ei(F::S, Y, B::Y);
      for(OT k=grid->si(F::S, Z, B::Y); k<=grid->ei(F::S, Z, B::Y); ++k)
        rhs(i, j, k) = 0.;
    }

    // BC XZ
    if(grid->getBCLocal()->getBCL(X)>0 && grid->getBCLocal()->getBCL(Z)>0) {
      OT i = grid->si(F::S, X, B::Y);
      OT k = grid->si(F::S, Z, B::Y);
      for(OT j=grid->si(F::S, Y, B::Y); j<=grid->ei(F::S, Y, B::Y); ++j)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCL(X)>0 && grid->getBCLocal()->getBCU(Z)>0) {
      OT i = grid->si(F::S, X, B::Y);
      OT k = grid->ei(F::S, Z, B::Y);
      for(OT j=grid->si(F::S, Y, B::Y); j<=grid->ei(F::S, Y, B::Y); ++j)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(X)>0 && grid->getBCLocal()->getBCL(Z)>0) {
      OT i = grid->ei(F::S, X, B::Y);
      OT k = grid->si(F::S, Z, B::Y);
      for(OT j=grid->si(F::S, Y, B::Y); j<=grid->ei(F::S, Y, B::Y); ++j)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(X)>0 && grid->getBCLocal()->getBCU(Z)>0) {
      OT i = grid->ei(F::S, X, B::Y);
      OT k = grid->ei(F::S, Z, B::Y);
      for(OT j=grid->si(F::S, Y, B::Y); j<=grid->ei(F::S, Y, B::Y); ++j)
        rhs(i, j, k) = 0.;
    }

    // BC YZ
    if(grid->getBCLocal()->getBCL(Y)>0 && grid->getBCLocal()->getBCL(Z)>0) {
      OT j = grid->si(F::S, Y, B::Y);
      OT k = grid->si(F::S, Z, B::Y);
      for(OT i=grid->si(F::S, X, B::Y); i<=grid->ei(F::S, X, B::Y); ++i)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCL(Y)>0 && grid->getBCLocal()->getBCU(Z)>0) {
      OT j = grid->si(F::S, Y, B::Y);
      OT k = grid->ei(F::S, Z, B::Y);
      for(OT i=grid->si(F::S, X, B::Y); i<=grid->ei(F::S, X, B::Y); ++i)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(Y)>0 && grid->getBCLocal()->getBCL(Z)>0) {
      OT j = grid->ei(F::S, Y, B::Y);
      OT k = grid->si(F::S, Z, B::Y);
      for(OT i=grid->si(F::S, X, B::Y); i<=grid->ei(F::S, X, B::Y); ++i)
        rhs(i, j, k) = 0.;
    }
    if(grid->getBCLocal()->getBCU(Y)>0 && grid->getBCLocal()->getBCU(Z)>0) {
      OT j = grid->ei(F::S, Y, B::Y);
      OT k = grid->ei(F::S, Z, B::Y);
      for(OT i=grid->si(F::S, X, B::Y); i<=grid->ei(F::S, X, B::Y); ++i)
        rhs(i, j, k) = 0.;
    }
    rhs.changed();
  }

public:

  DGProjector() {}
  DGProjector(const Teuchos::RCP<const OperatorT>& op)//:
    //nullspace_(op->grid())
  {

    //DivGradNullSpace<DivOp<SpaceT> > compNullspace;

    //compNullspace.computeNullSpace(op->getDivOp(), nullspace_);
    
    //setCornersZero(nullspace_);
    
    //Scalar blup = std::sqrt(1./nullspace_.dot(nullspace_));
    //nullspace_.scale(blup);
  }

  void operator()(RangeFieldT& rhs) const {

    //auto grid = rhs.grid();
    
    //ST bla = -nullspace_.dot(rhs);

    //if(0==grid->rankST())
      //std::cout << "DivGrad^-1"<< ": nullspace contributtion: " << std::abs(bla)  << "\n";

    //if(std::abs(bla)>0.)
      //rhs.add(1., rhs, bla, nullspace_);

    setCornersZero(rhs);
  }

}; // end of class EmptyProjector

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DGPROJECTOR_HPP
