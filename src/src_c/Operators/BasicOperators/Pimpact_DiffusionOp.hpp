/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_HELMHOLTZOP_HPP
#define PIMPACT_HELMHOLTDOP_HPP


#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {


/// \brief Helmholtz operator
/// \ingroup BaseOperator
///
/// computes \f$ y = (mulI_ I - mulL_ \Delta) x \f$
template<class SpT>
class DiffusionOp {

public:

  using GridT = SpT;

  using DomainFieldT = VectorField<GridT>;
  using RangeFieldT = VectorField<GridT>;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  static const int dimNC = SpT::dimNC;
  static const int dim = SpT::dimension;

  using SW = typename SpT::SW;

  using Stenc = Stencil<ST, OT, 0, SW::BL(0), SW::BU(0) >;

  using TO = const Teuchos::Tuple<Stenc, SpT::sdim>;

  const Teuchos::RCP<const GridT> grid_;

  ST mulI_;
  ST mulL_;

  TO cS_;
  TO cV_;

public:

  DiffusionOp(
    const Teuchos::RCP<const GridT>& grid):
    grid_(grid),
    mulI_(static_cast<ST>(0.)),
    mulL_(1./grid_->getDomainSize()->getRe()) {

    //const bool mapping = true; // order ~2
    const bool mapping = false; // order ~6

    for(int dir=0; dir<GridT::sdim; ++dir) {

      F fdir = static_cast<F>(dir);

      ST mulBC_ = 1.e+2; /// good idea, but has to be considered also in Force(what happens in Schurcomplement? Hopefully nothing

      // scalar stencil
      cS_[dir] = Stenc(grid_->nLoc(dir));
      FD_getDiffCoeff(
        1,
        grid_->nLoc(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->getBCLocal()->getBCL(dir),
        grid_->getBCLocal()->getBCU(dir),
        grid_->getShift(dir),
        5,
        dir+1,
        2,
        0,
        mapping, // mapping
        grid_->getStencilWidths()->getDimNcbC(dir),
        grid_->getStencilWidths()->getNcbC(dir),
        grid_->getCoordinatesLocal()->getX(F::S, dir),
        grid_->getCoordinatesLocal()->getX(F::S, dir),
        cS_[dir].get());

      // lower BC
      if(BC::Dirichlet==grid_->bcl(dir)) {
        for(int ii=SW::BL(X); ii<=SW::BU(X); ++ii)
          cS_[dir](1, ii) = 0.;

        cS_[dir](1, 0) = mulBC_;

      } else if(BC::Neumann==grid_->bcl(dir)) {
        FD_getDiffCoeff(
          1,
          1,
          grid_->bl(dir),
          grid_->bu(dir),
          grid_->bl(dir),
          grid_->bu(dir),
          grid_->getBCLocal()->getBCL(dir),
          0,
          grid_->getShift(dir),
          5,
          dir+1,
          1,
          0,
          mapping, // mapping
          grid_->getStencilWidths()->getDimNcbC(dir),
          grid_->getStencilWidths()->getNcbC(dir),
          grid_->getCoordinatesLocal()->getX(F::S, dir),
          grid_->getCoordinatesLocal()->getX(F::S, dir),
          cS_[dir].get());
      }
      // upper BC
      if(BC::Dirichlet==grid_->bcu(dir)) {
        for(int ii=SW::BL(X); ii<=SW::BU(X); ++ii)
          cS_[dir](grid_->nLoc(dir), ii) = 0.;
          cS_[dir](grid_->nLoc(dir), 0) = mulBC_;
      }
      if(BC::Neumann==grid_->bcu(dir)) {
        FD_getDiffCoeff(
          grid_->nLoc(dir),
          grid_->nLoc(dir),
          grid_->bl(dir),
          grid_->bu(dir),
          grid_->bl(dir),
          grid_->bu(dir),
          0,
          grid_->getBCLocal()->getBCU(dir),
          grid_->getShift(dir),
          5,
          dir+1,
          1,
          0,
          mapping, // mapping
          grid_->getStencilWidths()->getDimNcbC(dir),
          grid_->getStencilWidths()->getNcbC(dir),
          grid_->getCoordinatesLocal()->getX(F::S, dir),
          grid_->getCoordinatesLocal()->getX(F::S, dir),
          cS_[dir].get());
      }


      // velocity stencil
      cV_[dir] = Stenc(grid_->nLoc(dir));
      FD_getDiffCoeff(
        0,
        grid_->nLoc(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->getBCLocal()->getBCL(dir),
        grid_->getBCLocal()->getBCU(dir),
        grid_->getShift(dir),
        1,
        dir+1,
        2,
        0,
        mapping,
        grid_->getStencilWidths()->getDimNcbC(dir),
        grid_->getStencilWidths()->getNcbC(dir),
        grid_->getCoordinatesLocal()->getX(fdir, dir),
        grid_->getCoordinatesLocal()->getX(fdir, dir),
        cV_[dir].get());

      if(BC::Dirichlet==grid_->bcl(dir)) {
        for(OT ii=SW::BL(dir); ii<=SW::BU(dir); ++ii)
          cV_[dir](0, ii) = 0.;

        for(OT ii=SW::DL(dir); ii<=SW::DU(dir); ++ii)
          cV_[dir](0, ii+1) = mulBC_*grid()->getInterpolateV2S()->getC(dir, 1, ii);
          //cV_[dir](0, ii+1) = grid()->getInterpolateV2S()->getC(dir, 1, ii);
      } else if(BC::Neumann==grid_->bcl(dir)) {

        using StencD = Stencil<ST, OT, 0, SW::DL(0), SW::DU(0) >;

        StencD c_(grid_->nLoc(dir));

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
          c_.get());

        for(OT ii=SW::BL(dir); ii<=SW::BU(dir); ++ii)
          cV_[dir](0, ii) = 0.;

        for(OT ii=SW::DL(dir); ii<=SW::DU(dir); ++ii)
          cV_[dir](0, ii+1) = c_(1, ii);
      }


      if(BC::Dirichlet==grid_->bcu(dir)) {
        for(OT ii=SW::BL(dir); ii<=SW::BU(dir); ++ii)
          cV_[dir](grid()->ei(fdir, dir, B::Y), ii) = 0.;

        for(OT ii=SW::DL(dir); ii<=SW::DU(dir); ++ii)
          cV_[dir](grid()->ei(fdir, dir, B::Y), ii) =
            mulBC_*grid()->getInterpolateV2S()->getC(dir, grid()->ei(F::S, dir, B::Y), ii);
            //grid()->getInterpolateV2S()->getC(dir, grid()->ei(F::S, dir, B::Y), ii);
      } else if(BC::Neumann==grid_->bcu(dir)) {
        using StencD = Stencil<ST, OT, 0, SW::DL(0), SW::DU(0) >;

        StencD c_(grid_->nLoc(dir));

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
          c_.get());

        for(OT ii=SW::BL(dir); ii<=SW::BU(dir); ++ii)
          cV_[dir](grid()->ei(fdir, dir, B::Y), ii) = 0.;

        for(OT ii=SW::DL(dir); ii<=SW::DU(dir); ++ii)
          cV_[dir](grid()->ei(fdir, dir, B::Y), ii) = c_(grid()->ei(F::S, dir, B::Y), ii);
      }
    }
  };



  void apply(const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N ) const {

    const B nb = B::N;

    for(int dir=0; dir<GridT::sdim; ++dir) {

      F f = static_cast<F>(dir);

      x(f).exchange();

      if(3==GridT::sdim) {
        for(OT k=grid()->si(f, Z, nb); k<=grid()->ei(f, Z, nb); ++k)
          for(OT j=grid()->si(f, Y, nb); j<=grid()->ei(f, Y, nb); ++j)
            for(OT i=grid()->si(f, X, nb); i<=grid()->ei(f, X, nb); ++i) {
              if(Add::N==add) y(f)(i, j, k) = 0.;
              y(f)(i, j, k) +=
                mulI_*x(f)(i, j, k) - mulL_*innerStenc3D(x(f), f, i, j, k);
            }
      } else {
        for(OT k=grid()->si(f, Z, nb); k<=grid()->ei(f, Z, nb); ++k)
          for(OT j=grid()->si(f, Y, nb); j<=grid()->ei(f, Y, nb); ++j)
            for(OT i=grid()->si(f, X, nb); i<=grid()->ei(f, X, nb); ++i) {
              if(Add::N==add) y(f)(i, j, k) = 0.;
              y(f)(i, j, k) +=
                mulI_*x(f)(i, j, k) - mulL_*innerStenc2D(x(f), f, i, j, k);
            }
      }
    }
    if(Add::N==add) applyBC(x, y);
    y.changed();
  }


  void applyBC(const VectorField<GridT>& x, VectorField<GridT>& y) const {
    for(F field=F::U; field<GridT::sdim; ++field)
      applyBC(x(field), y(field));
  }


  /// \brief implements Dirichlet boundary conditions as identity in tangential
  /// velocity direction or interpolation in wand normal direction
  void applyBC(const ScalarField<GridT>& x, ScalarField<GridT>& y) const {

    assert(x.getType()==y.getType());

    const F f = x.getType();

    // U-field
    if(F::U==f) {

      // tangential direction: Y
      if(0<grid()->bcl(Y)) {
        OT j = grid()->si(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            y(i, j, k) = 0.;
            for(OT jj=0; jj<=SW::BU(Y); ++jj)
              y(i, j, k) += getC(Y, f, j, jj)*x(i, j+jj, k);
          }
      }
      if(0<grid()->bcu(Y)) {
        OT j = grid()->ei(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            y(i, j, k) = 0.;
            for(OT jj=SW::BL(Y); jj<=0; ++jj)
              y(i, j, k) += getC(Y, f, j, jj)*x(i, j+jj, k);
          }
      }

      // tangential direction: Z
      if(0<grid()->bcl(Z)) {
        OT k = grid()->si(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            y(i, j, k) = 0.;
            for(OT kk=0; kk<=SW::BU(Z); ++kk)
              y(i, j, k) += getC(Z, f, k, kk)*x(i, j, k+kk);
          }
      }
      if(0<grid()->bcu(Z)) {
        OT k = grid()->ei(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            y(i, j, k) = 0.;
            for(OT kk=SW::BL(Z); kk<=0; ++kk)
              y(i, j, k) += getC(Z, f, k, kk)*x(i, j, k+kk);
          }
      }

      // normal direction: X
      if(0<grid()->bcl(X)) {
        OT i = grid()->si(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            y(i, j, k) = 0.;
            for(OT ii=0; ii<=SW::BU(X); ++ii)
              y(i, j, k) += getC(X, f, i, ii)*x(i+ii, j, k);
          }
      }
      if(0<grid()->bcu(X)) {
        OT i = grid()->ei(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            y(i, j, k) = 0.;
            for(OT ii=SW::BL(X); ii<=0; ++ii)
              y(i, j, k) += getC(X, f, i, ii)*x(i+ii, j, k);
          }
      }
    }

    // V-field
    if(F::V==f) {

      // tangential direction: X
      if(0<grid()->bcl(X)) {
        OT i = grid()->si(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j) {
            y(i, j, k) = 0.;
            for(OT ii=0; ii<=SW::BU(X); ++ii)
              y(i, j, k) += getC(X, f, i, ii)*x(i+ii, j, k);
          }
      }
      if(0<grid()->bcu(X)) {
        OT i = grid()->ei(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j) {
            y(i, j, k) = 0.;
            for(OT ii=SW::BL(X); ii<=0; ++ii)
              y(i, j, k) += getC(X, f, i, ii)*x(i+ii, j, k);
          }
      }

      // tangential direction: Z
      if(0<grid()->bcl(Z)) {
        OT k = grid()->si(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            y(i, j, k) = 0.;
            for(OT kk=0; kk<=SW::BU(Z); ++kk)
              y(i, j, k) += getC(Z, f, k, kk)*x(i, j, k+kk);
          }
      }
      if(0<grid()->bcu(Z)) {
        OT k = grid()->ei(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            y(i, j, k) = 0.;
            for(OT kk=SW::BL(Z); kk<=0; ++kk)
              y(i, j, k) += getC(Z, f, k, kk)*x(i, j, k+kk);
          }
      }

      // normal direction: Y
      if(0<grid()->bcl(Y)) {
        OT j = grid()->si(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            y(i, j, k) = 0.;
            for(OT jj=0; jj<=SW::BU(Y); ++jj)
              y(i, j, k) += getC(Y, f, j, jj)*x(i, j+jj, k);
          }
      }
      if(0<grid()->bcu(Y)) {
        OT j = grid()->ei(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            y(i, j, k) = 0.;
            for(OT jj=SW::BL(Y); jj<=0; ++jj)
              y(i, j, k) += getC(Y, f, j, jj)*x(i, j+jj, k);
          }
      }

    }

    // W-field
    if(F::W==f) {

      // tangential direction: X
      if(0<grid()->bcl(X)) {
        OT i = grid()->si(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            y(i, j, k) = 0.;
            for(OT ii=0; ii<=SW::BU(X); ++ii)
              y(i, j, k) += getC(X, f, i, ii)*x(i+ii, j, k);
          }
      }
      if(0<grid()->bcu(X)) {
        OT i = grid()->ei(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            y(i, j, k) = 0.;
            for(OT ii=SW::BL(X); ii<=0; ++ii)
              y(i, j, k) += getC(X, f, i, ii)*x(i+ii, j, k);
          }
      }

      // tangential direction: Y
      if(0<grid()->bcl(Y)) {
        OT j = grid()->si(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            y(i, j, k) = 0.;
            for(OT jj=0; jj<=SW::BU(Y); ++jj)
              y(i, j, k) += getC(Y, f, j, jj)*x(i, j+jj, k);
          }
      }
      if(0<grid()->bcu(Y)) {
        OT j = grid()->ei(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            y(i, j, k) = 0.;
            for(OT jj=SW::BL(Y); jj<=0; ++jj)
              y(i, j, k) += getC(Y, f, j, jj)*x(i, j+jj, k);
          }
      }

      // normal direction: Z
      if(0<grid()->bcl(Z)) {
        OT k = grid()->si(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            y(i, j, k) = 0.;
            for(OT kk=0; kk<=SW::BU(Z); ++kk)
              y(i, j, k) += getC(Z, f, k, kk)*x(i, j, k+kk);
          }
      }
      if(0<grid()->bcu(Z)) {
        OT k = grid()->ei(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            y(i, j, k) = 0.;
            for(OT kk=SW::BL(Z); kk<=0; ++kk)
              y(i, j, k) += getC(Z, f, k, kk)*x(i, j, k+kk);
          }
      }
    }
  }


  void assignField(const DomainFieldT& mv) {};

  bool hasApplyTranspose() const {
    return false;
  }

  void print(std::ostream& out=std::cout) const {

    out << "--- " << getLabel() << " ---\n";
    out << " --- scalar stencil: ---";

    for(int dir=0; dir<GridT::sdim; ++dir) {

      out << "\ncoord: " << static_cast<ECoord>(dir) << "\n";

      cS_[dir].print(out);
    }
    out << " --- velocity stencil: ---";

    for(int dir=0; dir<GridT::sdim; ++dir) {

      out << "\ncoord: " << static_cast<ECoord>(dir) << "\n";

      cV_[dir].print(out);
    }
  }


  constexpr const Teuchos::RCP<const GridT>& grid() {
    return grid_;
  };

  constexpr const ST* getC(const int dir, const F ftype) {
    return (dir==ftype)?cV_[dir].get():cS_[dir].get();
  }

  constexpr ST getC(const ECoord dir, const F ftype, OT i, OT ii) {
    return (static_cast<int>(dir)==static_cast<int>(ftype))?
      cV_[dir](i, ii) : cS_[dir](i, ii);
  }

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    if(para->name()!="Linear Solver") {
      mulI_ = para->get<ST>("mulI", 0.);
      mulL_ = para->get<ST>("mulL", 1./grid_->getDomainSize()->getRe());
    }
  }

  const std::string getLabel() const {
    return "Helmholtz";
  };

protected:

  constexpr ST innerStenc3DU(const ScalarField<GridT>& x,
      const OT i, const OT j, const OT k) {

    ST lap = 0.;

    for(int ii=SW::BL(X); ii<=SW::BU(X); ++ii)
      lap += cV_[X](i, ii)*x(i+ii, j, k);

    for(int jj=SW::BL(Y); jj<=SW::BU(Y); ++jj)
      lap += cS_[Y](j, jj)*x(i, j+jj, k);

    for(int kk=SW::BL(Z); kk<=SW::BU(Z); ++kk)
      lap += cS_[Z](k, kk)*x(i, j, k+kk);

    return lap;
  }

  constexpr ST innerStenc3DV(const ScalarField<GridT>& x,
      const OT i, const OT j, const OT k) {

    ST lap = 0.;

    for(int ii=SW::BL(X); ii<=SW::BU(X); ++ii)
      lap += cS_[X](i, ii)*x(i+ii, j, k);

    for(int jj=SW::BL(Y); jj<=SW::BU(Y); ++jj)
      lap += cV_[Y](j, jj)*x(i, j+jj, k);

    for(int kk=SW::BL(Z); kk<=SW::BU(Z); ++kk)
      lap += cS_[Z](k, kk)*x(i, j, k+kk);

    return lap;
  }

  constexpr ST innerStenc3DW(const ScalarField<GridT>& x,
      const OT i, const OT j, const OT k) {

    ST lap = 0.;

    for(int ii=SW::BL(X); ii<=SW::BU(X); ++ii)
      lap += cS_[X](i, ii)*x(i+ii, j, k);

    for(int jj=SW::BL(Y); jj<=SW::BU(Y); ++jj)
      lap += cS_[Y](j, jj)*x(i, j+jj, k);

    for(int kk=SW::BL(Z); kk<=SW::BU(Z); ++kk)
      lap += cV_[Z](k, kk)*x(i, j, k+kk);

    return lap;
  }

public:

  constexpr ST innerStenc3D(const ScalarField<GridT>& x, const F fType, const OT i,
      const OT j, const OT k) {

    return (F::U==fType)?innerStenc3DU(x, i, j, k):((F::V==fType)?innerStenc3DV(x, i, j, k):innerStenc3DW(x, i, j, k));
  }

  constexpr ST innerStenc2D(const ScalarField<GridT>& x, const F fType, const OT i,
      const OT j, const OT k) {

    ST lap = 0.;

    for(int ii=SW::BL(X); ii<=SW::BU(X); ++ii)
      lap += getC(X, fType, i, ii)*x(i+ii, j, k);

    for(int jj=SW::BL(Y); jj<=SW::BU(Y); ++jj)
      lap += getC(Y, fType, j, jj)*x(i, j+jj, k);

    return lap;
  }

  constexpr ST innerDiag3D(const F fType,
                            const OT i, const OT j, const OT k) const {

    return getC(X, fType, i, 0) + getC(Y, fType, j, 0) + getC(Z, fType, k, 0);
  }

  constexpr ST innerDiag2D(const F fType, const OT i, const OT j, const OT k) {

    return getC(X, fType, i, 0) + getC(Y, fType, j, 0);
  }

}; // end of class DiffusionOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_HELMHOLTZOP_HPP
