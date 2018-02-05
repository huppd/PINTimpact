#pragma once
#ifndef PIMPACT_GRADOP_HPP
#define PIMPACT_GRADOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \ingroup BaseOperator
template<class SpT>
class GradOp {

public:

  using GridT = SpT;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using SW = typename GridT::SW;

  static const int dimNC = SpT::dimNC;
  static const int dim = SpT::dimension;

  using StencD = Stencil<ST, OT, 0, SW::DL(0), SW::DU(0) >;
  using StencG = Stencil<ST, OT, 0, SW::GL(0), SW::GU(0) >;

  using TD = const Teuchos::Tuple<StencD, SpT::sdim>;
  using TG = const Teuchos::Tuple<StencG, SpT::sdim>;

  Teuchos::RCP<const GridT> grid_;

  TG c_;
  TD cT_;

public:

  static constexpr int epsI = 1e3;

  using DomainFieldT =  ScalarField<GridT>;
  using RangeFieldT =  VectorField<GridT>;

  /// \todo for the transposed stencil it would be better to compute locally all global stencils instead of MPI_allreduce them
  /// \todo make it MG readey (participant, all reduce)
  GradOp(const Teuchos::RCP<const GridT>& grid):
    grid_(grid) {

    //const bool mapping=true; // order ~2
    const bool mapping=false; // order ~6

    for(int dir=0; dir<SpT::sdim; ++dir) {
      // Gradient stencil

      F fdir = static_cast<F>(dir);

      c_[dir] = StencG(grid_->nLoc(dir));

      FD_getDiffCoeff(
        0,
        grid_->nLoc(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->gl(dir),
        grid_->gu(dir),
        grid_->getBCLocal()->getBCL(dir),
        grid_->getBCLocal()->getBCU(dir),
        grid_->getShift(dir),
        2,
        dir+1,
        1,
        0,
        mapping,
        grid_->getStencilWidths()->getDimNcbG(dir),
        grid_->getStencilWidths()->getNcbG(dir),
        grid_->getCoordinatesLocal()->getX(F::S, dir),
        grid_->getCoordinatesLocal()->getX(fdir, dir),
        c_[dir].get());

      // transposed Gradient stencil
      cT_[dir] = StencD(grid_->nLoc(dir));

      OT nTempG = (grid_->nGlo(dir) + grid_->bu(dir) - grid_->bl(dir) + 1)
                  *(grid_->bu(dir) - grid_->bl(dir) + 1);


      Stencil<ST, OT, SW::BL(0), SW::BL(0), SW::BU(0) >
      cG1(grid_->nGlo(dir) + grid_->bu(dir));
      Stencil<ST, OT, SW::BL(0), SW::BL(0), SW::BU(0) >
      cG2(grid_->nGlo(dir) + grid_->bu(dir));

      for(OT i = grid_->si(fdir, dir, B::Y); i<=grid_->ei(fdir, dir, B::Y); ++i)
        for(OT ii = grid_->gl(dir); ii<=grid_->gu(dir); ++ii)
          cG1(i+grid_->getShift(dir), ii) = getC(static_cast<ECoord>(dir), i, ii);

      MPI_Allreduce(
        cG1.get(),	                              // const void *sendbuf,
        cG2.get(),                                // void *recvbuf,
        nTempG,			                              // int count,
        MPI_REAL8,	                              // MPI_Datatype datatype,
        MPI_SUM,		                              // MPI_Op op,
        grid_->getProcGrid()->getCommBar(dir)); // MPI_Comm comm)

      if(-1==grid_->getBCGlobal()->getBCL(dir)) {

        OT ls1 = grid_->getStencilWidths()->getLS(dir);
        OT M1 = grid_->nGlo(dir);

        for(OT i=grid->bl(dir); i<=-1; ++i)
          for(OT ii=grid->bl(dir); ii<=grid->bu(dir); ++ii)
            cG2(2+ls1+i, ii) = cG2(M1+1+ls1+i, ii);

        for(OT i=1; i<=grid->bu(dir); ++i)
          for(OT ii=grid->bl(dir); ii<=grid->bu(dir); ++ii)
            cG2(M1+ls1+i, ii) = cG2(1+ls1+i, ii);
      }

      for(OT
           i =grid_->si(F::S, static_cast<ECoord>(dir), B::Y);
           i<=grid_->ei(F::S, static_cast<ECoord>(dir), B::Y);
           ++i)
        for(OT ii=grid->dl(dir); ii<=grid->du(dir); ++ii) {
          cT_[dir](i, ii) = cG2(i+ii+grid_->getShift(dir), -ii);
        }
    }
  };


  void apply(const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N) const {

    applyG(x, y, add);

    if(Add::N==add)
      applyJ(y);
  }


  void applyG(const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N) const {

    const B& b = ((Add::N==add) ? B::Y : B::N);
    //const B& b = B::N;

    x.exchange(X);
    for(OT k=grid()->si(F::U, Z, b); k<=grid()->ei(F::U, Z, b); ++k)
      for(OT j=grid()->si(F::U, Y, b); j<=grid()->ei(F::U, Y, b); ++j)
        for(OT i=grid()->si(F::U, X, b); i<=grid()->ei(F::U, X, b); ++i) {
          if(Add::N==add) y(F::U)(i, j, k) = 0.;
          y(F::U)(i, j, k) += innerStencU(x, i, j, k);
        }

    x.exchange(Y);
    for(OT k=grid()->si(F::V, Z, b); k<=grid()->ei(F::V, Z, b); ++k)
      for(OT j=grid()->si(F::V, Y, b); j<=grid()->ei(F::V, Y, b); ++j)
        for(OT i=grid()->si(F::V, X, b); i<=grid()->ei(F::V, X, b); ++i) {
          if(Add::N==add) y(F::V)(i, j, k) = 0.;
          y(F::V)(i, j, k) += innerStencV(x, i, j, k);
        }

    if(3==GridT::sdim)  {
      x.exchange(Z);
      for(OT k=grid()->si(F::W, Z, b); k<=grid()->ei(F::W, Z, b); ++k)
        for(OT j=grid()->si(F::W, Y, b); j<=grid()->ei(F::W, Y, b); ++j)
          for(OT i=grid()->si(F::W, X, b); i<=grid()->ei(F::W, X, b); ++i) {
            if(Add::N==add) y(F::W)(i, j, k) = 0.;
            y(F::W)(i, j, k) += innerStencW(x, i, j, k);
          }
    }
    y.changed();
  }


  void applyJ(RangeFieldT& y) const {

    // BC scaling
    const ST eps = 1./static_cast<ST>(epsI);

    for(F dir=F::U; dir<GridT::sdim; ++dir) {
      B bc2 = B::Y;
      if(F::U!=dir) {
        if(grid()->getBCLocal()->getBCL(X) > 0) {
          OT i = grid()->si(dir, X, B::Y);
          for(OT k=grid()->si(dir, Z, bc2); k<=grid()->ei(dir, Z, bc2); ++k)
            for(OT j=grid()->si(dir, Y, bc2); j<=grid()->ei(dir, Y, bc2); ++j)
              y(dir)(i, j, k) *= eps;
        }
        if(grid()->getBCLocal()->getBCU(X) > 0) {
          OT i = grid()->ei(dir, X, B::Y);
          for(OT k=grid()->si(dir, Z, bc2); k<=grid()->ei(dir, Z, bc2); ++k)
            for(OT j=grid()->si(dir, Y, bc2); j<=grid()->ei(dir, Y, bc2); ++j)
              y(dir)(i, j, k) *= eps;
        }
        bc2 = B::N;
      }

      if(F::V!=dir) {
        if(grid()->getBCLocal()->getBCL(Y) > 0) {
          OT j = grid()->si(dir, Y, B::Y);
          for(OT k=grid()->si(dir, Z, bc2); k<=grid()->ei(dir, Z, bc2); ++k)
            for(OT i=grid()->si(dir, X, bc2); i<=grid()->ei(dir, X, bc2); ++i)
              y(dir)(i, j, k) *= eps;
        }
        if(grid()->getBCLocal()->getBCU(Y) > 0) {
          OT j = grid()->ei(dir, Y, B::Y);
          for(OT k=grid()->si(dir, Z, bc2); k<=grid()->ei(dir, Z, bc2); ++k)
            for(OT i=grid()->si(dir, X, bc2); i<=grid()->ei(dir, X, bc2); ++i)
              y(dir)(i, j, k) *= eps;
        }
        bc2 = B::N;
      }

      if(F::W!=dir) {
        if(grid()->getBCLocal()->getBCL(Z) > 0) {
          OT k = grid()->si(dir, Z, B::Y);
          for(OT j=grid()->si(dir, Y, bc2); j<=grid()->ei(dir, Y, bc2); ++j)
            for(OT i=grid()->si(dir, X, bc2); i<=grid()->ei(dir, X, bc2); ++i)
              y(dir)(i, j, k) *= eps;
        }
        if(grid()->getBCLocal()->getBCU(Z) > 0) {
          OT k = grid()->ei(dir, Z, B::Y);
          for(OT j=grid()->si(dir, Y, bc2); j<=grid()->ei(dir, Y, bc2); ++j)
            for(OT i=grid()->si(dir, X, bc2); i<=grid()->ei(dir, X, bc2); ++i)
              y(dir)(i, j, k) *= eps;
        }
        bc2 = B::N;
      }
    }
    y.extrapolateBC();
    y.changed();
  }


  void apply(const RangeFieldT& x, DomainFieldT& y, const Add add=Add::N) const {

    for(int dir=0; dir<GridT::sdim; ++dir)
      x.exchange(dir, dir);

    if(3==GridT::sdim)  {

      for(OT k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
        for(OT j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
          for(OT i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
            if(Add::N==add) y(i, j, k) = 0.;
            y(i, j, k) += innerStenc3D(x, i, j, k);
          }
    } else {

      for(OT k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
        for(OT j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
          for(OT i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
            if(Add::N==add) y(i, j, k) = 0.;
            y(i, j, k) += innerStenc2D(x, i, j, k);
          }
    }

    y.changed();
  }



  void assignField(const RangeFieldT& mv) {};
  void assignField(const DomainFieldT& mv) {};

  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return grid_;
  };

  constexpr const ST* getC(const ECoord dir) const {
    return c_[dir].get();
  }

  constexpr ST getC(const ECoord dir, OT i, OT off) const {
    return c_[dir](i, off);
  }

  constexpr ST getCTrans(const ECoord dir, OT i, OT off) const {
    return cT_[dir](i, off);
  }

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {}

  void print(std::ostream& out=std::cout) const {
    out <<"\n--- " <<getLabel() <<" ---\n";
    //out <<" --- stencil: ---";
    for(int dir=0; dir<SpT::sdim; ++dir) {
      out <<"\ndir: " <<static_cast<ECoord>(dir) <<"\n";

      c_[dir].print(out);
    }

    //out <<"--- " <<getLabel() <<"^T ---\n";
    //out <<" --- stencil: ---";
    //for(int dir=0; dir<SpT::sdim; ++dir) {
    //out <<"\ndir: " <<static_cast<ECoord>(dir) <<"\n";

    //cT_[dir].print(out);
    //}
  }

  const std::string getLabel() const {
    return "Grad";
  };

protected:

  constexpr ST innerStencU(const DomainFieldT& x, const OT i, const OT j, const OT k) {

    ST grad = 0.;

    for(int ii=SW::GL(X); ii<=SW::GU(X); ++ii)
      grad += getC(X, i, ii)*x(i+ii, j, k);

    return grad;
  }

  constexpr ST innerStencV(const DomainFieldT& x, const OT i, const OT j, const OT k) {

    ST grad = 0.;

    for(int jj=SW::GL(Y); jj<=SW::GU(Y); ++jj)
      grad += getC(Y, j, jj)*x(i, j+jj, k);

    return grad;
  }

  constexpr ST innerStencW(const DomainFieldT& x, const OT i, const OT j, const OT k) {

    ST grad = 0.;

    for(int kk=SW::GL(Z); kk<=SW::GU(Z); ++kk)
      grad += getC(Z, k, kk)*x(i, j, k+kk);

    return grad;
  }

  constexpr ST innerStenc3D(const RangeFieldT& x, const OT i, const OT j, const OT k) {

    ST gradT = 0.;

    for(int ii=SW::DL(X); ii<=SW::DU(X); ++ii)
      gradT += getCTrans(X, i, ii)*x(F::U)(i+ii, j, k);

    for(int jj=SW::DL(Y); jj<=SW::DU(Y); ++jj)
      gradT += getCTrans(Y, j, jj)*x(F::V)(i, j+jj, k);

    for(int kk=SW::DL(Z); kk<=SW::DU(Z); ++kk)
      gradT += getCTrans(Z, k, kk)*x(F::W)(i, j, k+kk);

    return gradT;
  }

  constexpr ST innerStenc2D(const RangeFieldT& x, const OT i, const OT j, const OT k) {

    ST gradT = 0.;

    for(int ii=SW::DL(X); ii<=SW::DU(X); ++ii)
      gradT += getCTrans(X, i, ii)*x(F::U)(i+ii, j, k);

    for(int jj=SW::DL(Y); jj<=SW::DU(Y); ++jj)
      gradT += getCTrans(Y, j, jj)*x(F::V)(i, j+jj, k);

    return gradT;
  }


}; // end of class GradOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRADOP_HPP
