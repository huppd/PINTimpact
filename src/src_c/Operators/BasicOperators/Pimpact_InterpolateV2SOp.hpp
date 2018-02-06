#pragma once
#ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
#define PIMPACT_INTERPOLATEVTOSOP_HPP


#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_CoordinatesLocal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {




template<class S, class O, int sd, int d, int dimNC>
class Grid;



template<class GridT>
class ScalarField;



/// \brief Interpolation operator.
/// \ingroup BaseOperator
/// \ingroup Grid
///
/// is used in the \c ScalarField::write method to interpolate the velocity to the pressure points, also used in \c ConvectionVOp
template<class Scalar, class Ordinal, int sdim, int dimension, int dimNC>
class InterpolateV2S {

public:

  using GridT = Grid<Scalar, Ordinal, sdim, dimension, dimNC>;

  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

protected:

  using SW = StencilWidths<dimension, dimNC>;

  using Stenc = Stencil<Scalar, Ordinal, 0, SW::DL(0), SW::DU(0) >;
  using TS = const Teuchos::Tuple<Stenc, sdim >;

  TS c_;
  //TS cm_;

public:

  InterpolateV2S(
    const Teuchos::RCP<const IndexSpace<Ordinal, dimension>>& indexSpace,
    const Teuchos::RCP<const GridSizeLocal<Ordinal, sdim, dimension> >& gridSizeLocal,
    const Teuchos::RCP<const StencilWidths<dimension, dimNC> >& stencilWidths,
    const Teuchos::RCP<const DomainSize<Scalar, sdim> >& domainSize,
    const Teuchos::RCP<const BoundaryConditionsLocal<dimension> >& boundaryConditionsLocal,
    const Teuchos::RCP<const CoordinatesLocal<Scalar, Ordinal, dimension, dimNC> >& coordinatesLocal) {

    //const bool mapping = true; // order ~4
    const bool mapping = false; // order ~6

    for(int i=0; i<sdim; ++i) {

      c_[i]  = Stenc(gridSizeLocal->get(i));

      F fi = static_cast<F>(i);

      FD_getDiffCoeff(
        1,
        gridSizeLocal->get(i),
        stencilWidths->getBL(i),
        stencilWidths->getBU(i),
        stencilWidths->getDL(i),
        stencilWidths->getDU(i),
        boundaryConditionsLocal->getBCL(i),
        boundaryConditionsLocal->getBCU(i),
        indexSpace->getShift(i),
        3,
        i+1,  // direction
        0,    // 0-derivative
        0,    // central
        mapping, // mapping, works with interpolateV2S
        stencilWidths->getDimNcbD(i),
        stencilWidths->getNcbD(i),
        coordinatesLocal->getX(fi, i),
        coordinatesLocal->getX(F::S, i),
        c_[i].get());
    }
  };


  void apply(const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N) const {

    assert(x.getType() != F::S);
    assert(y.getType() == F::S);

    Teuchos::RCP<const GridT> grid = x.grid();

    ECoord m = static_cast<ECoord>(x.getType());

    x.exchange(m);

    if(X==m) {
      for(Ordinal k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
        for(Ordinal j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
          for(Ordinal i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
            if(Add::N==add) y(i, j, k) = 0.;
            for(Ordinal ii=c_[m].bl(); ii<=c_[m].bu(); ++ii)
              y(i, j, k) += getC(m, i, ii)*x(i+ii, j, k);
          }
    }

    if(Y==m) {
      for(Ordinal k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
        for(Ordinal j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
          for(Ordinal i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
            if(Add::N==add) y(i, j, k) = 0.;
            for(Ordinal jj=c_[m].bl(); jj<=c_[m].bu(); ++jj)
              y(i, j, k) += getC(m, j, jj)*x(i, j+jj, k);
          }
    }

    if(Z==m) {
      for(Ordinal k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
        for(Ordinal j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
          for(Ordinal i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
            if(Add::N==add) y(i, j, k) = 0.;
            for(Ordinal kk=c_[m].bl(); kk<=c_[m].bu(); ++kk)
              y(i, j, k) += getC(m, k, kk)*x(i, j, k+kk);
          }
    }

    y.changed();
  }


  void assignField(const RangeFieldT& mv) {};

  bool hasApplyTranspose() const {
    return false;
  }

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  void print(std::ostream& out=std::cout) const {
    out << "--- " << getLabel() << " ---\n";
    for(int dir=0; dir<sdim; ++dir) {
      out << "\ndir: " << toString(static_cast<ECoord>(dir)) << "\n";
      c_[dir].print(out);
    }
  }

  constexpr const Scalar* getC(const int dir) const  {
    return c_[dir].get();
  }
  //constexpr const Scalar* getCM(const int dir) const  {
  //return cm_[dir].get();
  //}
  constexpr const Scalar getC(const int dir, Ordinal i, Ordinal off) const {
    return c_[dir](i, off);
  }
  //constexpr const Scalar getCM(const int dir, Ordinal i, Ordinal off) const {
  //return cm_[dir](i, off);
  //}

  const std::string getLabel() const {
    return "InterpolateV2S";
  };

};



/// \relates InterpolateV2S
template<class S, class O, int sd, int d, int dimNC >
Teuchos::RCP<const InterpolateV2S<S, O, sd, d, dimNC> > createInterpolateV2S(
    const Teuchos::RCP<const IndexSpace<O, d> >&  iS,
    const Teuchos::RCP<const GridSizeLocal<O, sd, d> >& gridSizeLocal,
    const Teuchos::RCP<const StencilWidths<d, dimNC> >& stencilWidths,
    const Teuchos::RCP<const DomainSize<S, sd> >& domainSize,
    const Teuchos::RCP<const BoundaryConditionsLocal<d> >& boundaryConditionsLocal,
    const Teuchos::RCP<const CoordinatesLocal<S, O, d, dimNC> >& coordinatesLocal) {

  return Teuchos::rcp(
      new InterpolateV2S<S, O, sd, d, dimNC>(
        iS,
        gridSizeLocal,
        stencilWidths,
        domainSize,
        boundaryConditionsLocal,
        coordinatesLocal));
}



/// \relates InterpolateV2S
template<class S, class O, int sd, int d, int dimNC >
Teuchos::RCP<const InterpolateV2S<S, O, sd, d, dimNC> > createInterpolateV2S(
    const Teuchos::RCP<const Grid<S, O, sd, d, dimNC> >& grid) {

  return grid->getInterpolateV2S();
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
