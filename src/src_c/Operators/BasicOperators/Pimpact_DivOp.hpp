#pragma once
#ifndef PIMPACT_DIVOP_HPP
#define PIMPACT_DIVOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \brief Divergence operator.
/// \ingroup BaseOperator
template<class ST>
class DivOp {

public:

  using GridT = ST;

  using DomainFieldT = VectorField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

protected:

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  static const int dimNC = ST::dimNC;
  static const int dim = ST::dimension;

  using StencD = Stencil<Scalar, Ordinal, 0, ST::SW::DL(0), ST::SW::DU(0) >;
  using StencG = Stencil<Scalar, Ordinal, 0, ST::SW::GL(0), ST::SW::GU(0) >;

  using TD = const Teuchos::Tuple<StencD, ST::sdim >;
  using TG = const Teuchos::Tuple<StencG, ST::sdim >;

  Teuchos::RCP<const GridT> grid_;

  TD c_;
  TG cT_;

public:

  /// \todo make it MG readey (participant, all reduce)
  DivOp(const Teuchos::RCP<const GridT>& grid):
    grid_(grid) {

    //const bool mapping = true;  // order ~2
    const bool mapping = false;   // order ~6

    for(int dir=0; dir<ST::sdim; ++dir) {
      F fdir = static_cast<F>(dir);

      // Divergence stencil
      c_[dir] = StencD(grid_->nLoc(dir));

      FD_getDiffCoeff(
        1,
        grid_->nLoc(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->dl(dir),
        grid_->du(dir),
        grid_->getBCLocal()->getBCL(dir),
        grid_->getBCLocal()->getBCU(dir),
        grid_->getShift(dir),
        3,
        dir+1,
        1,
        0,
        mapping, // mapping
        grid_->getStencilWidths()->getDimNcbD(dir),
        grid_->getStencilWidths()->getNcbD(dir),
        grid_->getCoordinatesLocal()->getX(fdir, dir),
        grid_->getCoordinatesLocal()->getX(F::S, dir),
        c_[dir].get());

      // Divergence stencil transposed
      cT_[dir] = StencG(grid_->nLoc(dir));

      Ordinal nTempG = (grid_->nGlo(dir) + grid_->bu(dir) - grid_->bl(dir) + 1)
                       *(grid_->bu(dir) - grid_->bl(dir) + 1);

      Stencil<Scalar, Ordinal, ST::SW::BL(0), ST::SW::BL(0), ST::SW::BU(0) >
      cG1(grid_->nGlo(dir) + grid_->bu(dir));
      Stencil<Scalar, Ordinal, ST::SW::BL(0), ST::SW::BL(0), ST::SW::BU(0) >
      cG2(grid_->nGlo(dir) + grid_->bu(dir));

      for(Ordinal i = grid_->si(F::S, dir); i<=grid_->ei(F::S, dir); ++i)
        for(Ordinal ii = grid_->dl(dir); ii<=grid_->du(dir); ++ii)
          cG1(i+grid_->getShift(dir), ii)= getC(static_cast<ECoord>(dir), i, ii);

      MPI_Allreduce(
        cG1.get(),    		                        // const void *sendbuf,
        cG2.get(),    		                        // void *recvbuf,
        nTempG,			                              // int count,
        MPI_REAL8,	                              // MPI_Datatype datatype,
        MPI_SUM,		                              // MPI_Op op,
        grid_->getProcGrid()->getCommBar(dir)); // MPI_Comm comm)


      if(-1==grid_->getBCGlobal()->getBCL(dir)) {

        Ordinal ls1 = grid_->getStencilWidths()->getLS(dir);
        Ordinal M1 = grid_->nGlo(dir);

        for(Ordinal i=grid->bl(dir); i<=-1; ++i)
          for(Ordinal ii=grid->bl(dir); ii<=grid->bu(dir); ++ii)
            cG2(2+ls1+i, ii) = cG2(M1+1+ls1+i, ii);

        for(Ordinal i=1; i<=grid->bu(dir); ++i)
          for(Ordinal ii=grid->bl(dir); ii<=grid->bu(dir); ++ii)
            cG2(M1+ls1+i, ii) = cG2(1+ls1+i, ii);
      }

      for(Ordinal i = grid_->si(fdir, dir, B::Y);
           i<=grid_->ei(fdir, dir, B::Y); ++i)
        for(Ordinal ii=grid->gl(dir); ii<=grid->gu(dir); ++ii)
          cT_[dir](i, ii) = cG2(i+ii+grid_->getShift(dir), -ii);
    }
  };



  void apply(const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N) const {

    for(int dir=0; dir<ST::sdim; ++dir)
      x.exchange(dir, dir);

    if(3==ST::sdim) {

      for(Ordinal k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
        for(Ordinal j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
          for(Ordinal i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
            if(Add::N==add) y(i, j, k) = 0.;
            y(i, j, k) += innerStenc3D(x, i, j, k);
          }
    } else {

      for(Ordinal k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
        for(Ordinal j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
          for(Ordinal i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
            if(Add::N==add) y(i, j, k) = 0.;
            y(i, j, k) += innerStenc2D(x, i, j, k);
          }
    }

    y.changed();
  }


  void apply(const RangeFieldT& x, DomainFieldT& y, const Add add=Add::N) const {

    x.exchange(X);
    for(Ordinal k=grid()->si(F::U, Z, B::Y); k<=grid()->ei(F::U, Z, B::Y); ++k)
      for(Ordinal j=grid()->si(F::U, Y, B::Y); j<=grid()->ei(F::U, Y, B::Y); ++j)
        for(Ordinal i=grid()->si(F::U, X, B::Y); i<=grid()->ei(F::U, X, B::Y); ++i) {
          if(Add::N==add) y(F::U)(i, j, k) = 0.;
          y(F::U)(i, j, k) += innerStencU(x, i, j, k);
        }

    x.exchange(Y);
    for(Ordinal k=grid()->si(F::V, Z, B::Y); k<=grid()->ei(F::V, Z, B::Y); ++k)
      for(Ordinal j=grid()->si(F::V, Y, B::Y); j<=grid()->ei(F::V, Y, B::Y); ++j)
        for(Ordinal i=grid()->si(F::V, X, B::Y); i<=grid()->ei(F::V, X, B::Y); ++i) {
          if(Add::N==add) y(F::V)(i, j, k) = 0.;
          y(F::V)(i, j, k) += innerStencV(x, i, j, k);
        }

    if(3==GridT::sdim)  {

      x.exchange(Z);
      for(Ordinal k=grid()->si(F::W, Z, B::Y); k<=grid()->ei(F::W, Z, B::Y); ++k)
        for(Ordinal j=grid()->si(F::W, Y, B::Y); j<=grid()->ei(F::W, Y, B::Y); ++j)
          for(Ordinal i=grid()->si(F::W, X, B::Y); i<=grid()->ei(F::W, X, B::Y); ++i) {
            if(Add::N==add) y(F::W)(i, j, k) = 0.;
            y(F::W)(i, j, k) += innerStencW(x, i, j, k);
          }
    }

    y.extrapolateBC(Belos::TRANS);

    // BC scaling
    const Scalar eps = 0.1;

    for(F dir=F::U; dir<GridT::sdim; ++dir) {
      B bc2 = B::Y;
      if(F::U!=dir) {
        if(grid()->getBCLocal()->getBCL(X) > 0) {
          Ordinal i = grid()->si(dir, X, B::Y);
          for(Ordinal k=grid()->si(dir, Z, bc2); k<=grid()->ei(dir, Z, bc2); ++k)
            for(Ordinal j=grid()->si(dir, Y, bc2); j<=grid()->ei(dir, Y, bc2); ++j)
              y(dir)(i, j, k) *= eps;
        }
        if(grid()->getBCLocal()->getBCU(X) > 0) {
          Ordinal i = grid()->ei(dir, X, B::Y);
          for(Ordinal k=grid()->si(dir, Z, bc2); k<=grid()->ei(dir, Z, bc2); ++k)
            for(Ordinal j=grid()->si(dir, Y, bc2); j<=grid()->ei(dir, Y, bc2); ++j)
              y(dir)(i, j, k) *= eps;
        }
        bc2 = B::N;
      }

      if(F::V!=dir) {
        if(grid()->getBCLocal()->getBCL(Y) > 0) {
          Ordinal j = grid()->si(dir, Y, B::Y);
          for(Ordinal k=grid()->si(dir, Z, bc2); k<=grid()->ei(dir, Z, bc2); ++k)
            for(Ordinal i=grid()->si(dir, X, bc2); i<=grid()->ei(dir, X, bc2); ++i)
              y(dir)(i, j, k) *= eps;
        }
        if(grid()->getBCLocal()->getBCU(Y) > 0) {
          Ordinal j = grid()->ei(dir, Y, B::Y);
          for(Ordinal k=grid()->si(dir, Z, bc2); k<=grid()->ei(dir, Z, bc2); ++k)
            for(Ordinal i=grid()->si(dir, X, bc2); i<=grid()->ei(dir, X, bc2); ++i)
              y(dir)(i, j, k) *= eps;
        }
        bc2 = B::N;
      }

      if(F::W!=dir) {
        if(grid()->getBCLocal()->getBCL(Z) > 0) {
          Ordinal k = grid()->si(dir, Z, B::Y);
          for(Ordinal j=grid()->si(dir, Y, bc2); j<=grid()->ei(dir, Y, bc2); ++j)
            for(Ordinal i=grid()->si(dir, X, bc2); i<=grid()->ei(dir, X, bc2); ++i)
              y(dir)(i, j, k) *= eps;
        }
        if(grid()->getBCLocal()->getBCU(Z) > 0) {
          Ordinal k = grid()->ei(dir, Z, B::Y);
          for(Ordinal j=grid()->si(dir, Y, bc2); j<=grid()->ei(dir, Y, bc2); ++j)
            for(Ordinal i=grid()->si(dir, X, bc2); i<=grid()->ei(dir, X, bc2); ++i)
              y(dir)(i, j, k) *= eps;
        }
        bc2 = B::N;
      }
    }

    y.changed();
  }


  void assignField(const RangeFieldT& mv) const {};
  void assignField(const DomainFieldT& mv) const {};

  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return grid_;
  };

  constexpr const Scalar* getC(const ECoord dir) const {
    return c_[dir].get();
  }

  constexpr const Scalar getC(const ECoord dir, Ordinal i, Ordinal off) const {
    return c_[dir](i, off);
  }

  constexpr const Scalar getCTrans(const ECoord dir, Ordinal i, Ordinal off) const {
    return cT_[dir](i, off);
  }

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  void print(std::ostream& out=std::cout) const {
    out <<"\n--- " <<getLabel() <<" ---\n";
    out <<" --- stencil: ---";
    for(int dir=0; dir<ST::sdim; ++dir) {
      out <<"\ndir: " <<toString(static_cast<ECoord>(dir)) <<"\n";
      c_[dir].print(out);
    }

    //out <<"--- " <<getLabel() <<"^T ---\n";
    //out <<" --- stencil: ---";
    //for(int dir=0; dir<ST::sdim; ++dir) {
    //out <<"\ndir: " <<toString(static_cast<ECoord>(dir)) <<"\n\n";
    //cT_[dir].print(out);
    //}
  }

  const std::string getLabel() const {
    return "Div";
  };

protected:

  constexpr Scalar innerStenc3D(const DomainFieldT& x, const Ordinal i, const Ordinal j,
      const Ordinal k) const {

    Scalar div = 0.;

    for(int ii=grid_->dl(X); ii<=grid_->du(X); ++ii)
      div += getC(X, i, ii)*x(F::U)(i+ii, j, k);

    for(int jj=grid_->dl(Y); jj<=grid_->du(Y); ++jj)
      div += getC(Y, j, jj)*x(F::V)(i, j+jj, k);

    for(int kk=grid_->dl(Z); kk<=grid_->du(Z); ++kk)
      div += getC(Z, k, kk)*x(F::W)(i, j, k+kk);

    return div;
  }

  constexpr Scalar innerStenc2D(const DomainFieldT& x,
                                 const Ordinal i, const Ordinal j, const Ordinal k) const {

    Scalar div = 0.;

    for(int ii=grid_->dl(X); ii<=grid_->du(X); ++ii)
      div += getC(X, i, ii)*x(F::U)(i+ii, j, k);

    for(int jj=grid_->dl(Y); jj<=grid_->du(Y); ++jj)
      div += getC(Y, j, jj)*x(F::V)(i, j+jj, k);

    return div;
  }

  constexpr Scalar innerStencU(const RangeFieldT& x,
                                const Ordinal i, const Ordinal j, const Ordinal k) const {

    Scalar divT = 0.;

    for(int ii=grid_->gl(X); ii<=grid_->gu(X); ++ii)
      divT += getCTrans(X, i, ii)*x(i+ii, j, k);

    return divT;
  }

  constexpr Scalar innerStencV(const RangeFieldT& x,
                                const Ordinal i, const Ordinal j, const Ordinal k) const {

    Scalar divT = 0.;

    for(int jj=grid_->gl(Y); jj<=grid_->gu(Y); ++jj)
      divT += getCTrans(Y, j, jj)*x(i, j+jj, k);

    return divT;
  }

  constexpr Scalar innerStencW(const RangeFieldT& x,
                                const Ordinal i, const Ordinal j, const Ordinal k) const {

    Scalar divT = 0.;

    for(int kk=grid_->gl(Z); kk<=grid_->gu(Z); ++kk)
      divT += getCTrans(Z, k, kk)*x(i, j, k+kk);

    return divT;
  }

}; // end of class DivOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DIVOP_HPP
