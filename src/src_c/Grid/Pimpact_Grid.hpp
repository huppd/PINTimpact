/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_SPACE_HPP
#define PIMPACT_SPACE_HPP


#include <iostream>

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"


#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"
#include "Pimpact_CoordinatesGlobal.hpp"
#include "Pimpact_CoordinatesLocal.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_InterpolateV2SOp.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_GridFactory.hpp"




/// \defgroup GridObject Grid Objects
///
/// overloaded class managing indexing, grid ...



namespace Pimpact {



extern "C" {

  void openH5F();
  void closeH5F();

}



/// \brief Grid in the sense of a vector grid, it is the connection between Field and Operators
///
/// \tparam ST scalar type normaly double
/// \tparam OT ordinal tye
/// \tparam sd spatial dimension 2 or 3
/// \tparam d  computational grid dimension 3 or 4
/// \tparam dNC dimension stencil
/// \ingroup GridObject
template<class ST, class OT, int sd, int d, int dNC >
class Grid {

public:

  using Scalar  = ST;
  using Ordinal = OT;

  static const int sdim = sd;
  static const int dimension = d;
  static const int dimNC = dNC;

  using SW = StencilWidths<dimension, dimNC>;

  /// \brief most important object
  ///
  /// \param pl parameters
  /// \note todo document parameter choices
  Grid(Teuchos::RCP<Teuchos::ParameterList> pl) {

    static_assert(sd!=2 || sd!=3, "spatial dimension not valid");

    pl->validateParametersAndSetDefaults(*getValidParameters(), 0);

    stencilWidths_ = Pimpact::createStencilWidths<dimension, dimNC>(
                       pl->get<bool>("spectral in time"));

    domainSize_ = Pimpact::createDomainSize<ST, sdim>(
                    pl->get<ST>("Re"),
                    pl->get<ST>("alpha2"),
                    pl->get<ST>("lx"),
                    pl->get<ST>("ly"),
                    pl->get<ST>("lz"),
                    pl->get<ST>("origin x", 0.),
                    pl->get<ST>("origin y", 0.),
                    pl->get<ST>("origin z", 0.));

    // are all template paramter needed here?
    //int domain = pl->get<int>("domain");
    if(2==sd) {
      pl->sublist("boundary conditions").set<int>("lower Z", -1);
      pl->sublist("boundary conditions").set<int>("upper Z", -1);
    }

    boundaryConditionsGlobal_ =
      Pimpact::createBoudaryConditionsGlobal<d>(
        Teuchos::sublist(pl, "boundary conditions"));

    gridSizeGlobal_ =
      Pimpact::createGridSizeGlobal<OT, sdim>(
        pl->get<OT>("nx"),
        pl->get<OT>("ny"),
        (2==sd)?2:pl->get<OT>("nz"),
        pl->get<OT>("nf"));

    Teuchos::Tuple<OT, d> procGridSize;
    procGridSize[0] = pl->get<OT>("npx");
    procGridSize[1] = pl->get<OT>("npy");
    procGridSize[2] = pl->get<OT>("npz");
    if(d>3)
      procGridSize[3] = pl->get<OT>("npf");

    procGrid_ =
      Pimpact::createProcGrid<OT, d>(
        procGridSize,
        boundaryConditionsGlobal_);

    boundaryConditionsLocal_ =
      Pimpact::createBoudaryConditionsLocal(
        boundaryConditionsGlobal_,
        procGrid_);

    gridSizeLocal_ =
      Pimpact::createGridSizeLocal<OT, sdim, d, dNC>(
        gridSizeGlobal_,
        procGrid_,
        stencilWidths_);

    indexSpace_ =
      Pimpact::createIndexSpace<OT, sdim, d, dNC>(
        stencilWidths_,
        gridSizeLocal_,
        boundaryConditionsLocal_,
        procGrid_);

    coordGlobal_ =
      Pimpact::createCoordinatesGlobal<ST, OT, sdim, d>(
        gridSizeGlobal_,
        domainSize_,
        Teuchos::tuple<Teuchos::RCP<Teuchos::ParameterList> >(
          Teuchos::sublist(pl, "Stretching in X"),
          Teuchos::sublist(pl, "Stretching in Y"),
          Teuchos::sublist(pl, "Stretching in Z")));

    coordLocal_ =
      Pimpact::createCoordinatesLocal<ST, OT, sdim, d, dNC>(
        stencilWidths_,
        domainSize_,
        gridSizeGlobal_,
        gridSizeLocal_,
        boundaryConditionsGlobal_,
        boundaryConditionsLocal_,
        procGrid_,
        coordGlobal_);

    interV2S_ =
      Pimpact::createInterpolateV2S<ST, OT, sdim, dimension, dimNC>(
        indexSpace_,
        gridSizeLocal_,
        stencilWidths_,
        domainSize_,
        boundaryConditionsLocal_,
        coordLocal_);

    openH5F();

    Teuchos::writeParameterListToXmlFile(*pl, "grid.xml");
  }


  Grid(
    const Teuchos::RCP<const StencilWidths<dimension, dimNC> >& stencilWidths,
    const Teuchos::RCP<const IndexSpace<OT, dimension> >& indexSpace,
    const Teuchos::RCP<const GridSizeGlobal<OT, sd> >& gridSizeGlobal,
    const Teuchos::RCP<const GridSizeLocal<OT, sd, dimension> >& gridSizeLocal,
    const Teuchos::RCP<const ProcGrid<OT, dimension> >& procGrid,
    const Teuchos::RCP<const CoordinatesGlobal<ST, OT, dimension> >& coordGlobal,
    const Teuchos::RCP<const CoordinatesLocal<ST, OT, dimension, dimNC> >& coordLocal,
    const Teuchos::RCP<const DomainSize<ST, sd> > domainSize,
    const Teuchos::RCP<const BoundaryConditionsGlobal<dimension> > boundaryConditionsGlobal,
    const Teuchos::RCP<const BoundaryConditionsLocal<dimension> > boundaryConditionsLocal,
    const Teuchos::RCP<const InterpolateV2S<ST, OT, sd, dimension, dimNC> >& interV2S):
    stencilWidths_(stencilWidths),
    indexSpace_(indexSpace),
    gridSizeGlobal_(gridSizeGlobal),
    gridSizeLocal_(gridSizeLocal),
    procGrid_(procGrid),
    coordGlobal_(coordGlobal),
    coordLocal_(coordLocal),
    domainSize_(domainSize),
    boundaryConditionsGlobal_(boundaryConditionsGlobal),
    boundaryConditionsLocal_(boundaryConditionsLocal),
    interV2S_(interV2S) {}

  //  ~Grid(){ closeH5F(); }

protected:

  Teuchos::RCP<const StencilWidths<dimension, dimNC> >
  stencilWidths_;

  Teuchos::RCP<const IndexSpace<OT, dimension> >
  indexSpace_;

  Teuchos::RCP<const GridSizeGlobal<OT, sd> >
  gridSizeGlobal_;

  Teuchos::RCP<const GridSizeLocal<OT, sd, dimension> >
  gridSizeLocal_;

  Teuchos::RCP<const ProcGrid<OT, dimension> >
  procGrid_;

  Teuchos::RCP<const CoordinatesGlobal<ST, OT, dimension> >
  coordGlobal_;

  Teuchos::RCP<const CoordinatesLocal<ST, OT, dimension, dimNC> >
  coordLocal_;

  Teuchos::RCP<const DomainSize<ST, sd> > domainSize_;

  Teuchos::RCP<const BoundaryConditionsGlobal<dimension> >
  boundaryConditionsGlobal_;

  Teuchos::RCP<const BoundaryConditionsLocal<dimension> >
  boundaryConditionsLocal_;

  Teuchos::RCP<const InterpolateV2S<ST, OT, sd, dimension, dimNC> >
  interV2S_;

public:

  /// \name getter methods
  /// \{

  constexpr const Teuchos::RCP<const StencilWidths<dimension, dimNC> >&
  getStencilWidths() const {
    return stencilWidths_;
  }

  constexpr const Teuchos::RCP<const IndexSpace<OT, dimension> >&
  getIndexSpace() const {
    return indexSpace_;
  }

  constexpr const Teuchos::RCP<const GridSizeGlobal<OT, sd> >&
  getGridSizeGlobal() const {
    return gridSizeGlobal_;
  }

  constexpr const Teuchos::RCP<const GridSizeLocal<OT, sd, dimension> >&
  getGridSizeLocal() const {
    return gridSizeLocal_;
  }

  constexpr const Teuchos::RCP<const ProcGrid<OT, dimension> >&
  getProcGrid() const {
    return procGrid_;
  }

  constexpr const Teuchos::RCP<const CoordinatesGlobal<ST, OT, dimension> >&
  getCoordinatesGlobal() const {
    return coordGlobal_;
  }

  constexpr const Teuchos::RCP<const CoordinatesLocal<ST, OT, dimension, dimNC> >&
  getCoordinatesLocal() const {
    return coordLocal_;
  }

  constexpr const Teuchos::RCP<const DomainSize<ST, sd> >&
  getDomainSize() const {
    return domainSize_;
  }

  constexpr const Teuchos::RCP<const BoundaryConditionsGlobal<dimension> >&
  getBCGlobal()   const {
    return boundaryConditionsGlobal_;
  }

  constexpr const Teuchos::RCP<const BoundaryConditionsLocal<dimension> >&
  getBCLocal()    const {
    return boundaryConditionsLocal_;
  }

  constexpr const Teuchos::RCP<const InterpolateV2S<ST, OT, sd, dimension, dimNC> >&
  getInterpolateV2S() const {
    return interV2S_;
  }

  /// \}

  /// \name getter methods IMPACT style
  /// \{

  constexpr const MPI_Comm& comm()  const {
    return procGrid_->getCommS();
  }

  constexpr const MPI_Comm& commST()  const {
    return procGrid_->getCommWorld();
  }

  constexpr const int rankST() const {
    return procGrid_->getRank();
  }
  constexpr const int rankS () const {
    return procGrid_->getRankS();
  }

  constexpr const Ordinal* nGlo()        const {
    return gridSizeGlobal_->getRawPtr();
  }
  constexpr const Ordinal nGlo(const int i) const {
    return gridSizeGlobal_->get(i);
  }

  constexpr const Ordinal* nLoc()        const {
    return gridSizeLocal_->getRawPtr();
  }
  constexpr const Ordinal nLoc(const int i) const {
    return gridSizeLocal_->get(i);
  }

  constexpr const Ordinal* bl()         const {
    return stencilWidths_->getBL();
  }
  constexpr const Ordinal bl(const int i)  const {
    return stencilWidths_->getBL(i);
  }

  constexpr const Ordinal* bu()         const {
    return stencilWidths_->getBU();
  }
  constexpr const Ordinal bu(const int i)  const {
    return stencilWidths_->getBU(i);
  }

  constexpr const Ordinal* dl()         const {
    return stencilWidths_->getDL();
  }
  constexpr const Ordinal dl(const int i)  const {
    return stencilWidths_->getDL(i);
  }

  constexpr const Ordinal* du()         const {
    return stencilWidths_->getDU();
  }
  constexpr const Ordinal du(const int i)  const {
    return stencilWidths_->getDU(i);
  }

  constexpr const Ordinal* gl()         const {
    return stencilWidths_->getGL();
  }
  constexpr const Ordinal gl(const int i)  const {
    return stencilWidths_->getGL(i);
  }

  constexpr const Ordinal* gu()         const {
    return stencilWidths_->getGU();
  }
  constexpr const Ordinal gu(const int i)  const {
    return stencilWidths_->getGU(i);
  }

  constexpr const Ordinal* nl()         const {
    return stencilWidths_->getNL();
  }
  constexpr const Ordinal nl(const int i)  const {
    return stencilWidths_->getNL(i);
  }

  constexpr const Ordinal* nu()         const {
    return stencilWidths_->getNU();
  }
  constexpr const Ordinal nu(const int i)  const {
    return stencilWidths_->getNU(i);
  }

  constexpr const int bcl(const int dir) const {
    return getBCLocal()->getBCL(dir);
  }
  constexpr const int bcu(const int dir) const {
    return getBCLocal()->getBCU(dir);
  }

  /// \deprecated
  constexpr const Ordinal* sInd(const F fieldType) const {
    return indexSpace_->sInd(fieldType);
  }
  /// \deprecated
  constexpr const Ordinal* eInd( const F fieldType) const {
    return indexSpace_->eInd(fieldType);
  }

  /// \deprecated
  constexpr const Ordinal* sIndB(const F fieldType) const {
    return indexSpace_->sIndB(fieldType);
  }
  /// \deprecated
  constexpr const Ordinal* eIndB(const F fieldType) const {
    return indexSpace_->eIndB(fieldType);
  }


  constexpr const Ordinal si(const F fieldType, const int dir, const B withB=B::N) {
    return static_cast<bool>(withB)?
      indexSpace_->sIndB(fieldType, dir) : indexSpace_->sInd (fieldType, dir);
  }
  constexpr const Ordinal ei(const F fieldType, const int dir, const B withB=B::N) {
    return static_cast<bool>(withB)?
      indexSpace_->eIndB(fieldType, dir) : indexSpace_->eInd (fieldType, dir);
  }

  constexpr const Ordinal* ib() const {
    return procGrid_->getIB().getRawPtr();
  }

  constexpr const Ordinal* getShift()      const {
    return indexSpace_->getShift();
  }
  constexpr const Ordinal getShift(const int i) const {
    return indexSpace_->getShift(i);
  }

  constexpr const Ordinal* np()      const {
    return procGrid_->getNP().getRawPtr();
  }
  constexpr const Ordinal np(const int i) const {
    return procGrid_->getNP(i);
  }

  /// \}

  void print(std::ostream& out=std::cout) const {

    out << "\t---Grid: ---\n";

    stencilWidths_->print(out);

    out << "\t---GridSizeGlobal: ---\n";
    gridSizeGlobal_->print(out);

    out << "\t---GridSizeLocal: ---\n";
    gridSizeLocal_->print(out);

    indexSpace_->print(out);

    procGrid_->print(out);

    getBCGlobal()->print(out);

    getBCLocal()->print(out);

    //    coordGlobal_->print(out);
    //
    //		coordLocal_->print(out);
  }


  static Teuchos::RCP<const Teuchos::ParameterList>  getValidParameters()  {

    static Teuchos::RCP<const Teuchos::ParameterList> validPL;

    // Set all the valid parameters and their default values.
    if(is_null(validPL)) {

      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList("Space");
      pl->set<ST>("Re", 1., "Reynolds number");
      pl->set<ST>("alpha2", 1.,
                  "Womersley square alpha^2");
      // domain type
      pl->sublist("boundary conditions");
      pl->sublist("boundary conditions").set<int>("lower X", static_cast<int>(BC::Dirichlet));
      pl->sublist("boundary conditions").set<int>("upper X", static_cast<int>(BC::Dirichlet));
      pl->sublist("boundary conditions").set<int>("lower Y", static_cast<int>(BC::Dirichlet));
      pl->sublist("boundary conditions").set<int>("upper Y", static_cast<int>(BC::Dirichlet));
      pl->sublist("boundary conditions").set<int>("lower Z", static_cast<int>(BC::Dirichlet));
      pl->sublist("boundary conditions").set<int>("upper Z", static_cast<int>(BC::Dirichlet));

      //pl->set<int>("domain", 2,
      //"Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel");

      // domain size
      pl->set<ST>("lx", 1., "length in x-direction");
      pl->set<ST>("ly", 1., "length in y-direction");
      pl->set<ST>("lz", 1., "length in z-direction");

      pl->set<ST>("origin x", 0., "origin in x-direction");
      pl->set<ST>("origin y", 0., "origin in y-direction");
      pl->set<ST>("origin z", 0., "origin in z-direction");

      // grid size
      pl->set<bool>("spectral in time", false, "enables spectral time discretization");

      pl->set<OT>("nx", 33, "amount of grid points in x-direction: a*2**q+1");
      pl->set<OT>("ny", 33, "amount of grid points in y-direction: a*2**q+1");
      pl->set<OT>("nz", 33, "amount of grid points in z-direction: a*2**q+1");
      pl->set<OT>("nf", 4, "amount of grid points in f-direction");

      // grid stretching
      pl->sublist("Stretching in X");
      pl->sublist("Stretching in Y");
      pl->sublist("Stretching in Z");
      pl->sublist("Stretching in X").set<std::string>("Stretch Type", "none");
      pl->sublist("Stretching in Y").set<std::string>("Stretch Type", "none");
      pl->sublist("Stretching in Z").set<std::string>("Stretch Type", "none");

      // processor grid size
      pl->set<OT>("npx", 2, "amount of processors in x-direction");
      pl->set<OT>("npy", 2, "amount of processors in y-direction");
      pl->set<OT>("npz", 1, "amount of processors in z-direction");
      pl->set<OT>("npf", 1, "amount of processors in f-direction");

      validPL = pl;
    }
    return validPL;
  }


}; // end of class Grid



/// \relates Grid
template<class ST>
Teuchos::RCP<const ST>
create(Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList()) {

  return Teuchos::rcp(new ST(pl));
}



/// \relates Grid
/// \relates TransferOp
/// \relates CoarsenStrategy
/// \relates CoarsenStrategyGlobal
/// "same" grid with changing dimNC
template<class OGridT, class IGridT>
Teuchos::RCP<const OGridT >
createGrid(
  const Teuchos::RCP<const IGridT>& grid) {

  using Scalar    = typename OGridT::Scalar;
  using Ordinal   = typename OGridT::Ordinal;
  const int sdim  = OGridT::sdim;
  const int dim   = OGridT::dimension;
  const int dimNC = OGridT::dimNC;

  Teuchos::RCP<const StencilWidths<dim, dimNC > > stencilWidths =
    createStencilWidths<dim, dimNC >(grid->getStencilWidths()->spectralT());

  Teuchos::RCP<const DomainSize<Scalar, sdim> > domainSize = grid->getDomainSize();

  Teuchos::RCP<const BoundaryConditionsGlobal<dim> > boundaryConditionsGlobal = grid->getBCGlobal();

  Teuchos::RCP<const BoundaryConditionsLocal<dim> > boundaryConditionsLocal = grid->getBCLocal();

  Teuchos::RCP<const GridSizeGlobal<Ordinal, sdim> > gridSizeGlobal = grid->getGridSizeGlobal();

  Teuchos::RCP<const GridSizeLocal<Ordinal, sdim, dim> > gridSizeLocal = grid->getGridSizeLocal();

  Teuchos::RCP<const ProcGrid<Ordinal, dim> > procGrid = grid->getProcGrid();

  Teuchos::RCP<const IndexSpace<Ordinal, dim> > indexSpace =
    Pimpact::createIndexSpace<Ordinal, sdim, dim, dimNC>(
        stencilWidths,
        gridSizeLocal,
        boundaryConditionsLocal,
        procGrid);

  Teuchos::RCP<const CoordinatesGlobal<Scalar, Ordinal, dim> > coordGlobal = grid->getCoordinatesGlobal();

  Teuchos::RCP<const CoordinatesLocal<Scalar, Ordinal, dim, dimNC> >  coordLocal =
    Pimpact::createCoordinatesLocal(
        stencilWidths,
        domainSize,
        gridSizeGlobal,
        gridSizeLocal,
        boundaryConditionsGlobal,
        boundaryConditionsLocal,
        procGrid,
        coordGlobal);

  Teuchos::RCP<const InterpolateV2S<Scalar, Ordinal, sdim, dim, dimNC> > interV2S =
    Pimpact::createInterpolateV2S<Scalar, Ordinal, sdim, dim, dimNC>(
        indexSpace,
        gridSizeLocal,
        stencilWidths,
        domainSize,
        boundaryConditionsLocal,
        coordLocal);

  return Teuchos::rcp(
      new OGridT(
        stencilWidths,
        indexSpace,
        gridSizeGlobal,
        gridSizeLocal,
        procGrid,
        coordGlobal,
        coordLocal,
        domainSize,
        boundaryConditionsGlobal,
        boundaryConditionsLocal,
        interV2S));

} // end of createGrid



/// \brief creates coarse grid from fine and new GridSize
///
/// \tparam GridT type of grid
/// \param grid grid object
/// \param newGridSizeGlobal new grid size
///
/// \return new coarse grid
template<class GridT>
static Teuchos::RCP<const GridT > createGrid(
  const Teuchos::RCP<const GridT>& grid,
  const GridSizeGlobal<typename GridT::Ordinal, GridT::sdim>& newGridSizeGlobal) {

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;
  const int sdim = GridT::sdim;
  const int dim = GridT::dimension;
  const int dimNC = GridT::dimNC;

  Teuchos::RCP<const StencilWidths<dim, dimNC> > stencilWidths = grid->getStencilWidths();

  Teuchos::RCP<const DomainSize<Scalar, sdim> > domainSize = grid->getDomainSize();

  Teuchos::RCP<const BoundaryConditionsGlobal<dim> > boundaryConditionsGlobal = grid->getBCGlobal();
  Teuchos::RCP<const BoundaryConditionsLocal <dim> > boundaryConditionsLocal = grid->getBCLocal();


  Teuchos::RCP<const GridSizeGlobal<Ordinal, sdim> > gridSizeGlobal =
    createGridSizeGlobal<Ordinal, sdim>(
        newGridSizeGlobal);

  Teuchos::RCP<const ProcGrid<Ordinal, dim> > procGrid = grid->getProcGrid();

  Teuchos::RCP<const GridSizeLocal<Ordinal, sdim, dim> > gridSizeLocal =
    Pimpact::createGridSizeLocal<Ordinal, sdim, dim, dimNC>(
        gridSizeGlobal,
        procGrid,
        stencilWidths);

  Teuchos::RCP<const IndexSpace<Ordinal, dim> > indexSpace =
    Pimpact::createIndexSpace<Ordinal, sdim, dim, dimNC>(
        stencilWidths,
        gridSizeLocal,
        boundaryConditionsLocal,
        procGrid);

  Teuchos::RCP<const CoordinatesGlobal<Scalar, Ordinal, dim> > coordGlobal =
    Pimpact::createCoordinatesGlobal<Scalar, Ordinal, sdim, dim>(
        gridSizeGlobal,
        grid->getCoordinatesGlobal());

  Teuchos::RCP<const CoordinatesLocal<Scalar, Ordinal, dim, dimNC> > coordLocal =
    Pimpact::createCoordinatesLocal<Scalar, Ordinal, sdim, dim, dimNC>(
        stencilWidths,
        domainSize,
        gridSizeGlobal,
        gridSizeLocal,
        boundaryConditionsGlobal,
        boundaryConditionsLocal,
        procGrid,
        coordGlobal);

  Teuchos::RCP<const InterpolateV2S<Scalar, Ordinal, sdim, dim, dimNC> > interV2S =
    Pimpact::createInterpolateV2S<Scalar, Ordinal, sdim, dim, dimNC>(
        indexSpace,
        gridSizeLocal,
        stencilWidths,
        domainSize,
        boundaryConditionsLocal,
        coordLocal);

  return Teuchos::rcp(
      new GridT(
        stencilWidths,
        indexSpace,
        gridSizeGlobal,
        gridSizeLocal,
        procGrid,
        coordGlobal,
        coordLocal,
        domainSize,
        boundaryConditionsGlobal,
        boundaryConditionsLocal,
        interV2S));

} // end of createGrid



/// \note todo move ProcGrid createion to createProcGrid(...)
/// \note todo redo input GridSizeGlobal, ProcGridSize
/// \relates CoarsenStrategyGlobal
template<class GridT>
static Teuchos::RCP<const GridT > createGrid(
    const Teuchos::RCP<const GridT>& grid,
    const GridSizeGlobal<typename GridT::Ordinal, GridT::sdim>& newGridSizeGlobal,
    const Teuchos::Tuple<typename GridT::Ordinal, GridT::dimension>& npNew,
    Teuchos::Tuple<typename GridT::Ordinal, GridT::dimension>& stride,
    const Teuchos::Tuple<typename GridT::Ordinal, GridT::dimension>& npWorld,
    const Teuchos::Tuple<typename GridT::Ordinal, GridT::dimension>& ibWorld) {

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  const int sdim = GridT::sdim;
  const int dimension = GridT::dimension;
  const int dimNC = GridT::dimNC;

  using TO = typename Teuchos::Tuple<Ordinal, dimension>;

  Teuchos::RCP<const StencilWidths<dimension, dimNC> >
    stencilWidths = grid->getStencilWidths();

  Teuchos::RCP<const DomainSize<Scalar, sdim> > domainSize = grid->getDomainSize();
  Teuchos::RCP<const BoundaryConditionsGlobal<dimension> > boundaryConditionsGlobal = grid->getBCGlobal();


  // --- coarsen gridSizeGlobal ---
  auto gridSizeGlobal =
    createGridSizeGlobal<Ordinal, sdim>(newGridSizeGlobal);


  /// redo procGrid create new communicator, make procgrid from communicator
  TO np = grid->getProcGrid()->getNP();
  TO nGather;
  TO ib;
  for(Ordinal dir=0; dir<dimension; ++dir) {
    nGather[dir] = np[dir]/npNew[dir];
    stride[dir] *= nGather[dir];
    ib[dir]      = (ibWorld[dir] - 1)*npNew[dir]/npWorld[dir] + 1;
  }

  bool participating = false;
  MPI_Comm commWorld = grid->getProcGrid()->getCommWorld();
  MPI_Comm commSub ;//= grid->getProcGrid()->getCommS();
  //std::cout << "rank: " << grid->rankST() << "\tcomm: " << commWorld<< "\n";

  int rankWorld = grid->getProcGrid()->getRank();
  int rankSub = grid->getProcGrid()->getRank(); // necessary?

  Teuchos::Tuple<int, dimension> rankL;
  Teuchos::Tuple<int, dimension> rankU;
  for(int dir=0; dir<dimension; ++dir)
    MPI_Cart_shift(
        commWorld,      // communicator with Cartesian structure
        dir,            // coordinate dimension of shift
        stride[dir],    // displacement
        &rankL[dir],    // rank of source process
        &rankU[dir]);	// rank of destination process

  Ordinal gather_yes = 1;
  for(int i=0; i<dimension; ++i) {
    gather_yes *= nGather[i];
  }
  if(gather_yes>1) {
    int n = 1;
    for(int i=0; i<dimension; ++i)
      n *= npNew[i];
    int* newRanks = new int[n];

    TO rankCoord;

    for(int i=0; i<npNew[0]; ++i) {
      rankCoord[0] = (i*stride[0])%npWorld[0];
      for(int j=0; j<npNew[1]; ++j) {
        rankCoord[1] = (j*stride[1])%npWorld[1];
        for(int k=0; k<npNew[2]; ++k) {
          rankCoord[2] = (k*stride[2])%npWorld[2];
          if(4==dimension) {
            for(int l=0; l<npNew[3]; ++l) {
              rankCoord[3] = (l*stride[3])%npWorld[3];
              MPI_Cart_rank(
                  commWorld,									          // comm
                  rankCoord.getRawPtr(),                // processor coordinate
                  &newRanks[i
                  +j*npNew[0]
                  +k*npNew[0]*npNew[1]
                  +l*npNew[0]*npNew[1]*npNew[2] ]);   	// according rank to coordinate
              if(rankWorld==newRanks[i + j*npNew[0] + k*npNew[0]*npNew[1] + l*npNew[0]*npNew[1]*npNew[2] ])
                participating = true;
            }
          } else {
            MPI_Cart_rank(
                commWorld,									                    // comm
                rankCoord.getRawPtr(),                          // processor coordinate
                &newRanks[i+j*npNew[0]+k*npNew[0]*npNew[1] ]); // according rank to coordinate
            if(rankWorld==newRanks[i+j*npNew[0]+k*npNew[0]*npNew[1] ])
              participating = true;
          }
        }
      }
    }

    MPI_Comm commTemp;
    MPI_Group baseGroup, newGroup;

    MPI_Comm_group(commWorld, &baseGroup);
    MPI_Group_incl(baseGroup, n, newRanks, &newGroup);
    MPI_Comm_create(commWorld, newGroup, &commTemp);
    MPI_Group_free(&baseGroup);
    MPI_Group_free(&newGroup);

    Teuchos::Tuple<int, dimension> periodic = boundaryConditionsGlobal->periodic();

    if(participating) {

      MPI_Cart_create(
          commTemp,		          // communicator without Cartesian information
          dimension,            // number of dimensions
          npNew.getRawPtr(),    // number of processors in each dimension
          periodic.getRawPtr(),	// array for mpi to signal which dimension is periodic
          false,                // false means ranking is not reordered
          &commSub);           // new communicator with Cartesian information
      if(4==dimension) {
        MPI_Comm commTemp_;
        int temp[] = {1, 1, 1, 0};
        MPI_Cart_sub(commSub, temp, &commTemp_);
        MPI_Comm_free(&commSub);
        commSub = commTemp_;
      }
      MPI_Comm_free(&commTemp);
    } else
      commSub=MPI_COMM_NULL;

    delete[] newRanks;

  }

  if(commSub==MPI_COMM_NULL)
    rankSub = -1;
  else
    MPI_Comm_rank(commSub, &rankSub); // get rank

  Teuchos::RCP<const ProcGrid<Ordinal, dimension> > procGrid =
    Teuchos::rcp(
        new ProcGrid<Ordinal, dimension>(
          npNew,
          participating,
          commWorld,
          commSub,
          rankWorld,
          rankSub,
          ib,
          rankL,
          rankU));
  //procGrid->print();

  //std::cout << "rank: " << grid->rankST() << "\tcomm: " << commWorld<< "\n";

  Teuchos::RCP<const GridSizeLocal<Ordinal, sdim, dimension> > gridSizeLocal =
    Pimpact::createGridSizeLocal<Ordinal, sdim, dimension, dimNC>(
        gridSizeGlobal,
        procGrid,
        stencilWidths);

  Teuchos::RCP<const BoundaryConditionsLocal<dimension> > boundaryConditionsLocal =
    createBoudaryConditionsLocal<Ordinal, dimension>(
        boundaryConditionsGlobal,
        procGrid);

  Teuchos::RCP<const IndexSpace<Ordinal, dimension> > indexSpace =
    Pimpact::createIndexSpace<Ordinal, sdim, dimension, dimNC>(
        stencilWidths,
        gridSizeLocal,
        boundaryConditionsLocal,
        procGrid);

  Teuchos::RCP<const CoordinatesGlobal<Scalar, Ordinal, dimension> > coordGlobal =
    Pimpact::createCoordinatesGlobal<Scalar, Ordinal, sdim, dimension>(
        gridSizeGlobal,
        grid->getCoordinatesGlobal());

  Teuchos::RCP<const CoordinatesLocal<Scalar, Ordinal, dimension, dimNC> > coordLocal =
    Pimpact::createCoordinatesLocal<Scalar, Ordinal, sdim, dimension, dimNC>(
        stencilWidths,
        domainSize,
        gridSizeGlobal,
        gridSizeLocal,
        boundaryConditionsGlobal,
        boundaryConditionsLocal,
        procGrid,
        coordGlobal);

  Teuchos::RCP<const InterpolateV2S<Scalar, Ordinal, sdim, dimension, dimNC> > interV2S =
    Pimpact::createInterpolateV2S<Scalar, Ordinal, sdim, dimension, dimNC>(
        indexSpace,
        gridSizeLocal,
        stencilWidths,
        domainSize,
        boundaryConditionsLocal,
        coordLocal);

  return Teuchos::rcp(
      new GridT(
        stencilWidths,
        indexSpace,
        gridSizeGlobal,
        gridSizeLocal,
        procGrid,
        coordGlobal,
        coordLocal,
        domainSize,
        boundaryConditionsGlobal,
        boundaryConditionsLocal,
        interV2S));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SPACE_HPP
