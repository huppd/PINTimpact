#pragma once
#ifndef PIMPACT_CONVECTIONSOP_HPP
#define PIMPACT_CONVECTIONSOP_HPP


#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
template<class ST>
class ConvectionSOp {

public:

  using GridT = ST;

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  using FluxFieldT = ScalarField<GridT>[3];
  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

protected:

  static const int dimNC = ST::dimNC;
  static const int dim = ST::dimension;
  static const int sdim = ST::sdim;

  using SW = StencilWidths<dim, dimNC>;

  using Stenc = Stencil<Scalar, Ordinal, 0, SW::NL(0), SW::NU(0) >;
  using TO = const Teuchos::Tuple<Stenc, sdim >;

  const Teuchos::RCP<const GridT> grid_;

  TO cSD_;
  TO cVD_;

  TO cSU_;
  TO cVU_;

public:

  ConvectionSOp(const Teuchos::RCP<const GridT>& grid ):
    grid_(grid) {

    //const bool mapping = true; // order: ~3
    const bool mapping = false; //order: ~5

    for(int dir=0; dir<sdim; ++dir) {

      cSD_[dir] = Stenc(grid_->nLoc(dir));

      FD_getDiffCoeff(
        1,
        grid_->nLoc(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->nl(dir),
        grid_->nu(dir),
        grid_->getBCLocal()->getBCL(dir),
        grid_->getBCLocal()->getBCU(dir),
        grid_->getShift(dir),
        5,
        dir+1,
        1,
        -1,
        mapping,
        grid_->getStencilWidths()->getDimNcbC(dir),
        grid_->getStencilWidths()->getNcbC(dir),
        grid_->getCoordinatesLocal()->getX(F::S, dir),
        grid_->getCoordinatesLocal()->getX(F::S, dir),
        cSD_[dir].get());

      if(BC::Dirichlet==grid_->bcl(dir)) {
        Ordinal i = grid_->si(F::S, dir, B::Y);
        for(int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii)
          cSD_[dir](i, ii) = 0.;
      }
      if(BC::Dirichlet==grid_->bcu(dir)) {
        Ordinal i = grid_->ei(F::S, dir, B::Y);
        for(int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii)
          cSD_[dir](i, ii) = 0.;
      }


      cSU_[dir] = Stenc(grid_->nLoc(dir));

      FD_getDiffCoeff(
        1,
        grid_->nLoc(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->nl(dir),
        grid_->nu(dir),
        grid_->getBCLocal()->getBCL(dir),
        grid_->getBCLocal()->getBCU(dir),
        grid_->getShift(dir),
        5,
        dir+1,
        1,
        +1,
        mapping,
        grid_->getStencilWidths()->getDimNcbC(dir),
        grid_->getStencilWidths()->getNcbC(dir),
        grid_->getCoordinatesLocal()->getX(F::S, dir),
        grid_->getCoordinatesLocal()->getX(F::S, dir),
        cSU_[dir].get());

      if(BC::Dirichlet==grid_->bcl(dir)) {
        Ordinal i = grid_->si(F::S, dir, B::Y);
        for(int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii)
          cSU_[dir](i, ii) = 0.;
      }
      if(BC::Dirichlet==grid_->bcu(dir)) {
        Ordinal i = grid_->ei(F::S, dir, B::Y);
        for(int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii)
          cSU_[dir](i, ii) = 0.;
      }


      F fdir = static_cast<F>(dir);

      cVD_[dir] = Stenc(grid_->nLoc(dir));

      FD_getDiffCoeff(
        0,
        grid_->nLoc(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->nl(dir),
        grid_->nu(dir),
        grid_->getBCLocal()->getBCL(dir),
        grid_->getBCLocal()->getBCU(dir),
        grid_->getShift(dir),
        1,
        dir+1,
        1,
        -1,
        mapping,
        grid_->getStencilWidths()->getDimNcbC(dir),
        grid_->getStencilWidths()->getNcbC(dir),
        grid_->getCoordinatesLocal()->getX(fdir, dir),
        grid_->getCoordinatesLocal()->getX(fdir, dir),
        cVD_[dir].get());

      if(BC::Dirichlet==grid_->bcl(dir)) {
        Ordinal i = grid_->si(fdir, dir, B::Y);
        for(int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii)
          cVD_[dir](i, ii) = 0.;
      }
      if(BC::Dirichlet==grid_->bcu(dir)) {
        Ordinal i = grid_->ei(fdir, dir, B::Y);
        for(int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii)
          cVD_[dir](i, ii) = 0.;
      }


      cVU_[dir] = Stenc(grid_->nLoc(dir));

      FD_getDiffCoeff(
        0,
        grid_->nLoc(dir),
        grid_->bl(dir),
        grid_->bu(dir),
        grid_->nl(dir),
        grid_->nu(dir),
        grid_->getBCLocal()->getBCL(dir),
        grid_->getBCLocal()->getBCU(dir),
        grid_->getShift(dir),
        1,
        dir+1,
        1,
        +1,
        mapping,
        grid_->getStencilWidths()->getDimNcbC(dir),
        grid_->getStencilWidths()->getNcbC(dir),
        grid_->getCoordinatesLocal()->getX(fdir, dir),
        grid_->getCoordinatesLocal()->getX(fdir, dir),
        cVU_[dir].get());

      if(BC::Dirichlet==grid_->bcl(dir)) {
        Ordinal i = grid_->si(fdir, dir, B::Y);
        for(int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii)
          cVU_[dir](i, ii) = 0.;
      }
      if(BC::Dirichlet==grid_->bcu(dir)) {
        Ordinal i = grid_->ei(fdir, dir, B::Y);
        for(int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii)
          cVU_[dir](i, ii) = 0.;
      }
    }
  };


  void assignField(const RangeFieldT& mv) {};


  void apply(const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, const Scalar mulI, const Scalar mulC, const Scalar mulL, const Add add=Add::N) const {
    apply(x, y, z, mulC, add);
  }


  void apply(const FluxFieldT& wind, const DomainFieldT& y, RangeFieldT& z, const Add add=Add::N) const {

    const B b = ((Add::N==add) ? B::Y : B::N);

    const Scalar mulC = 1.;
    F m = z.getType();

    assert(z.getType() == y.getType());

    for(int i=0; i<GridT::sdim; ++i)
      assert(wind[i].getType()==y.getType());


    for(int vel_dir=0; vel_dir<GridT::sdim; ++vel_dir)
      wind[vel_dir].exchange();

    y.exchange();

    if(3==GridT::sdim)
      for(Ordinal k=grid()->si(m, Z, b); k<=grid()->ei(m, Z, b); ++k)
        for(Ordinal j=grid()->si(m, Y, b); j<=grid()->ei(m, Y, b); ++j)
          for(Ordinal i=grid()->si(m, X, b); i<=grid()->ei(m, X, b); ++i) {
            if(Add::N==add) z(i, j, k) = 0.;
            z(i, j, k) += mulC*innerStenc3D(wind[0](i, j, k), wind[1](i, j, k), wind[2](i, j, k), y, i, j, k);
          }
    else
      for(Ordinal k=grid()->si(m, Z, b); k<=grid()->ei(m, Z, b); ++k)
        for(Ordinal j=grid()->si(m, Y, b); j<=grid()->ei(m, Y, b); ++j)
          for(Ordinal i=grid()->si(m, X, b); i<=grid()->ei(m, X, b); ++i) {
            if(Add::N==add) z(i, j, k) = 0.;
            z(i, j, k) += mulC*innerStenc2D(wind[0](i, j, k), wind[1](i, j, k), y, i, j, k);
          }

    z.changed();
  }

  void print(std::ostream& out=std::cout) const {
    out <<" --- ConvectioSOp ---\n";
    for(int i=0; i<GridT::sdim; ++i) {
      out <<"dir: " <<static_cast<ECoord>(i) <<"\n ";
      out <<"cSD:\n";
      cSD_[i].print(out);
      out <<"cSU:\n";
      cSU_[i].print(out);
      out <<"cVD:\n";
      cVD_[i].print(out);
      out <<"cVU:\n";
      cVU_[i].print(out);
    }
  }


  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const Teuchos::RCP<const GridT>&  grid() const {
    return grid_;
  }

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}
 

  constexpr const Scalar* getCU(const ECoord dir, const F ftype) {
    return (((int)dir)==((int)ftype))?cVU_[dir].get():cSU_[dir].get();
  }

  constexpr const Scalar* getCD(const ECoord dir, const F ftype) {
    return (((int)dir)==((int)ftype))?cVD_[dir].get():cSD_[dir].get();
  }

protected:

  constexpr Scalar getC(const Scalar wind, const ECoord dir, const F ftype, const int i, const int ii) {
    return (static_cast<int>(dir)==static_cast<int>(ftype))?
      (wind>=0? cVU_[dir](i, ii):cVD_[dir](i, ii)) :
      (wind>=0? cSU_[dir](i, ii):cSD_[dir](i, ii));
  }
 
  constexpr Scalar innerStenc3DU(const Scalar u, const Scalar v, const Scalar w, const
      RangeFieldT& x, const Ordinal i, const Ordinal j, const Ordinal k) const {

    Scalar dx = 0.;
    for(int ii=SW::NL(X); ii<=SW::NU(X); ++ii) {
      dx += (u>=0? cVU_[X](i, ii):cVD_[X](i, ii))*x(i+ii, j, k);
    }

    Scalar dy = 0.;
    for(int jj=SW::NL(Y); jj<=SW::NU(Y); ++jj)
      dy += (v>=0? cSU_[Y](j, jj):cSD_[Y](j, jj))*x(i, j+jj, k);

    Scalar dz = 0.;
    for(int kk=SW::NL(Z); kk<=SW::NU(Z); ++kk)
      dz += (w>=0? cSU_[Z](k, kk):cSD_[Z](k, kk))*x(i, j, k+kk);

    return u*dx+v*dy+w*dz;
  }

  constexpr Scalar innerStenc3DV(const Scalar u, const Scalar v, const Scalar w, const
      RangeFieldT& x, const Ordinal i, const Ordinal j, const Ordinal k) const {

    Scalar dx = 0.;
    for(int ii=SW::NL(X); ii<=SW::NU(X); ++ii) {
      dx += (u>=0? cSU_[X](i, ii):cSD_[X](i, ii))*x(i+ii, j, k);
    }

    Scalar dy = 0.;
    for(int jj=SW::NL(Y); jj<=SW::NU(Y); ++jj)
      dy += (v>=0? cVU_[Y](j, jj):cVD_[Y](j, jj))*x(i, j+jj, k);

    Scalar dz = 0.;
    for(int kk=SW::NL(Z); kk<=SW::NU(Z); ++kk)
      dz += (w>=0? cSU_[Z](k, kk):cSD_[Z](k, kk))*x(i, j, k+kk);

    return u*dx+v*dy+w*dz;
  }

  constexpr Scalar innerStenc3DW(const Scalar u, const Scalar v, const Scalar w, const
      RangeFieldT& x, const Ordinal i, const Ordinal j, const Ordinal k) const {

    Scalar dx = 0.;
    for(int ii=SW::NL(X); ii<=SW::NU(X); ++ii) {
      dx += (u>=0? cSU_[X](i, ii):cSD_[X](i, ii))*x(i+ii, j, k);
    }

    Scalar dy = 0.;
    for(int jj=SW::NL(Y); jj<=SW::NU(Y); ++jj)
      dy += (v>=0? cSU_[Y](j, jj):cSD_[Y](j, jj))*x(i, j+jj, k);

    Scalar dz = 0.;
    for(int kk=SW::NL(Z); kk<=SW::NU(Z); ++kk)
      dz += (w>=0? cVU_[Z](k, kk):cVD_[Z](k, kk))*x(i, j, k+kk);

    return u*dx+v*dy+w*dz;
  }

public:

  constexpr Scalar innerStenc2D(const Scalar u, const Scalar v, const RangeFieldT& x,
      const Ordinal i, const Ordinal j, const Ordinal k)
  const {

    Scalar dx = 0.;
    for(int ii=SW::NL(X); ii<=SW::NU(X); ++ii)
      dx += getC(u, X, x.getType(), i, ii)*x(i+ii, j, k);

    Scalar dy = 0.;
    for(int jj=SW::NL(Y); jj<=SW::NU(Y); ++jj)
      dy += getC(v, Y, x.getType(), j, jj)*x(i, j+jj, k);

    return u*dx+v*dy;
  }

  constexpr Scalar innerStenc3D(const Scalar u, const Scalar v, const Scalar w, const
      RangeFieldT& x, const Ordinal i, const Ordinal j, const Ordinal k) const {

    return (F::U==x.getType())?innerStenc3DU(u, v, w, x, i, j, k):((F::V==x.getType())?innerStenc3DV(u, v, w, x, i, j, k):innerStenc3DW(u, v, w, x, i, j, k));
  }

  constexpr Scalar innerDiag3D(const Scalar u, const Scalar v, const Scalar w, const F
      fType, const Ordinal i, const Ordinal j, const Ordinal k) const {

    return u*getC(u, X, fType, i, 0) + v*getC(v, Y, fType, j, 0) + w*getC(w, Z, fType, k, 0);
  }

  constexpr Scalar innerDiag2D(const Scalar u, const Scalar v, const F fType, const
      Ordinal i, const Ordinal j, const Ordinal k) const {

    return u*getC(u, X, fType, i, 0) + v*getC(v, Y, fType, j, 0);
  }

  constexpr const std::string getLabel() const {
    return "Convection";
  };


}; // end of class ConvectionSOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONSOP_HPP
