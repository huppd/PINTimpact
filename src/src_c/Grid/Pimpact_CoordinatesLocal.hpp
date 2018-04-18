/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_COORDINATESLOCAL_HPP
#define PIMPACT_COORDINATESLOCAL_HPP


#include <cmath>
#include <ostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"
#include "Pimpact_CoordinatesGlobal.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



extern "C"
void PI_getLocalCoordinates(
  const int& time,
  const double& L,
  const int& M,
  const int& N,
  const int& bL,
  const int& bU,
  const int& BC_L_global,
  const int& BC_U_global,
  const int& BC_L,
  const int& BC_U,
  const int& iB,
  const double* const ys,
  const double* const yv,
  double* const xs,
  double* const xv,
  double* const dxs,
  double* const dxv);



/// \brief local grid coordinates
/// \tparam ScalarT
/// \tparam OrdinalT
/// \tparam dim as soon as Time is own class ->sdim
/// \tparam dimNC
///
/// Coordinate | index
/// -----------| --------------------------
/// xS         | bl..nLoc+bu+(ib-1)*(NB-1)
/// xV         | bl..nLoc+bu+(ib-1)*(NB-1)
///
/// \todo make nice interface for getter
/// \ingroup GridObject
template<class ScalarT, class OrdinalT, int dim, int dimNC>
class CoordinatesLocal {

  template<class ST, class OT, int sdT, int dT, int dNC>
  friend Teuchos::RCP<const CoordinatesLocal<ST, OT, dT, dNC> > createCoordinatesLocal(
    const Teuchos::RCP<const StencilWidths<dT, dNC> >& stencilWidth,
    const Teuchos::RCP<const DomainSize<ST, sdT> >& domainSize,
    const Teuchos::RCP<const GridSizeGlobal<OT, sdT> >& gridSizeGlobal,
    const Teuchos::RCP<const GridSizeLocal<OT, sdT, dT> >& gridSizeLocal,
    const Teuchos::RCP<const BoundaryConditionsGlobal<dT> >& bcGlobal,
    const Teuchos::RCP<const BoundaryConditionsLocal<dT> >& bcLocal,
    const Teuchos::RCP<const ProcGrid<OT, dT> >& procGrid,
    const Teuchos::RCP<const CoordinatesGlobal<ST, OT, dT> >& coordGlobal);

protected:

  using SW = StencilWidths<dim, dimNC>;

  using DAS = Array<ScalarT, OrdinalT, 1>;
  using DAV = Array<ScalarT, OrdinalT, 0>;

  using TDAS = const Teuchos::Tuple<DAS, dim >;
  using TDAV = const Teuchos::Tuple<DAV, dim >;

  using AS = Array<ScalarT, OrdinalT, SW::BL(0)>;
  using AV = Array<ScalarT, OrdinalT, SW::BL(0)>;

  using TAS = const Teuchos::Tuple<AS, dim >;
  using TAV = const Teuchos::Tuple<AV, dim >;

  using TO = const Teuchos::Tuple<Teuchos::ArrayRCP<ScalarT>, dim >;

  TAS xS_;
  TAV xV_;

  TDAS dxS_;
  TDAV dxV_;

  Teuchos::RCP<const StencilWidths<dim, dimNC> > stencilWidths_;

  //template<int dimNC>
  template<int sdim>
  CoordinatesLocal(
    const Teuchos::RCP<const StencilWidths<dim, dimNC> >& stencilWidths,
    const Teuchos::RCP<const DomainSize<ScalarT, sdim> >& domainSize,
    const Teuchos::RCP<const GridSizeGlobal<OrdinalT, sdim> >& gridSizeGlobal,
    const Teuchos::RCP<const GridSizeLocal<OrdinalT, sdim, dim> >& gridSizeLocal,
    const Teuchos::RCP<const BoundaryConditionsGlobal<dim> >& bcGlobal,
    const Teuchos::RCP<const BoundaryConditionsLocal<dim> >& bcLocal,
    const Teuchos::RCP<const ProcGrid<OrdinalT, dim> >& procGrid,
    const Teuchos::RCP<const CoordinatesGlobal<ScalarT, OrdinalT, dim> >& coordGlobal):
    stencilWidths_(stencilWidths) {

    for(int i=0; i<dim; ++i) {

      xS_[i]  = AS(gridSizeLocal->get(i) + SW::BU(i));
      xV_[i]  = AV(gridSizeLocal->get(i) + SW::BU(i));
      dxS_[i] = DAS(std::max(gridSizeLocal->get(i), 1));
      dxV_[i] = DAV(std::max(gridSizeLocal->get(i), 1));

      F fi = static_cast<F>(i);

      if(i<3)
        PI_getLocalCoordinates(
          0,
          domainSize->getSize(i),
          gridSizeGlobal->get(i),
          gridSizeLocal->get(i),
          stencilWidths->getBL(i),
          stencilWidths->getBU(i),
          bcGlobal->getBCL(i),
          bcGlobal->getBCU(i),
          bcLocal->getBCL(i),
          bcLocal->getBCU(i),
          procGrid->getIB(i),
          coordGlobal->getX(F::S, i),
          coordGlobal->getX(fi, i),
          xS_[i].get(),
          xV_[i].get(),
          dxS_[i].get(),
          dxV_[i].get());
      else if(3==i) {
        ScalarT nt = gridSizeGlobal->get(i);
        OrdinalT offset = (procGrid->getIB(i)-1)*gridSizeLocal->get(i);
        ScalarT pi2 = 8.*std::atan(1.);
        for(OrdinalT ii = SW::BL(i); ii<=gridSizeLocal->get(i) + SW::BU(i); ++ii) {
          xS_[i][ii] = pi2/nt*(ii+offset -1);
          xV_[i][ii] = pi2/nt*(ii+offset -1);
        }
      }
    }
  }

public:

  /// \name getter
  /// @{

  constexpr const ScalarT* operator()(const F ftype, const int dir) {
    return getX(ftype, dir);
  }

  constexpr const ScalarT& operator()(const F ftype, const int dir, const OrdinalT& i) {
    return getX(ftype, dir, i);
  }

  /// \deprecated
  constexpr const ScalarT* getX(const F ftype, const int dir) {
    return (F::S==ftype || dir!=ftype) ?
      xS_[dir].get() : xV_[dir].get();
  }

  /// \deprecated
  constexpr const ScalarT& getX(const F ftype, const int dir, const OrdinalT& i) {
    return (F::S==ftype || dir!=ftype)?
      xS_[dir][i] : xV_[dir][i];
  }

  constexpr const ScalarT& dx(const F ftype, const int dir, const OrdinalT& i) {
    return (F::S==ftype || dir!=ftype)?
      dxS_[dir][i] : dxV_[dir][i];
  }

  /// \deprecated
  constexpr const AS& getS(const int dir) {
    return xS_[dir];
  }

  /// \deprecated
  constexpr const AV& getV(const int dir) {
    return xV_[dir];
  }

  ///  @}

  void print(std::ostream& out=std::cout) const {

    for(int i=0; i<dim; ++i) {
      out << "Local coordinates of scalars in dir: " << static_cast<ECoord>(i) << "\n";
      xS_[i].print(out);
    }

    for(int i=0; i<dim; ++i) {
      out << "Local coordinates of velocity in dir: " << static_cast<ECoord>(i) << "\n";
      xV_[i].print(out);
    }
  };

}; // end of class CoordinatesLocal



/// \brief create Grid coordinates Global
/// \relates CoordinatesLocal
template<class ST, class OT, int sd, int d, int dNC>
Teuchos::RCP<const CoordinatesLocal<ST, OT, d, dNC> >
createCoordinatesLocal(
  const Teuchos::RCP<const StencilWidths<d, dNC> >& stencilWidths,
  const Teuchos::RCP<const DomainSize<ST, sd> >& domainSize,
  const Teuchos::RCP<const GridSizeGlobal<OT, sd> >& gridSizeGlobal,
  const Teuchos::RCP<const GridSizeLocal<OT, sd, d> >& gridSizeLocal,
  const Teuchos::RCP<const BoundaryConditionsGlobal<d> >& bcGlobal,
  const Teuchos::RCP<const BoundaryConditionsLocal<d> >& bcLocal,
  const Teuchos::RCP<const ProcGrid<OT, d> >& procGrid,
  const Teuchos::RCP<const CoordinatesGlobal<ST, OT, d> >& coordGlobal) {

  return Teuchos::rcp(
      new CoordinatesLocal<ST, OT, d, dNC>(
        stencilWidths,
        domainSize,
        gridSizeGlobal,
        gridSizeLocal,
        bcGlobal,
        bcLocal,
        procGrid,
        coordGlobal));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COORDINATESLOCAL_HPP
