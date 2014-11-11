#pragma once
#ifndef PIMPACT_COARSENSTRATEGY_HPP
#define PIMPACT_COARSENSTRATEGY_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {



/// \brief first model implementation, where one coarsens until every processor has 3 dofs.
///
/// \tparam FSpaceT space on finest level not necessary the same for the coarse grids( difference in \c dim_nc).
/// it should be that the coardinates are taken from the fine grid and are halved.
/// \todo compute coordinates correctly (grid stretching only on finest grid, afterwards simple coarsening),
/// cleanest version would be to use same grid stretching on every level, makes
/// interpolation and restriction slightly more complicated.
///\todo the \c CSpaceT parameter allows even smaller Grids nLoc=1 but therefore
///the exception in GridSilzeLocal has to be adapted on StencilWidths
template<class FSpaceT,class CSpaceT>
class CoarsenStrategy {

  typedef typename FSpaceT::Scalar  Scalar;
  typedef typename FSpaceT::Ordinal Ordinal;

  /// should be same for finest space and coarse spaces
  static const int dimension = FSpaceT::dimension;

  /// can be different for finest and coarse spaces
  static const int dimNCF = FSpaceT::dimNC;

  static const int dimNCC = CSpaceT::dimNC;

public:

  static std::vector<Teuchos::RCP<const CSpaceT> > getMultiSpace(
      const Teuchos::RCP<const FSpaceT> space,
      int maxGrids=10 ) {

    // FSpace>CSpace if necessary
//    Teuchos::RCP<const CSpaceT> tempSpace;
//    Teuchos::RCP<const Space<typename CSpaceT::Scalar,typename CSpaceT::Ordinal,CSpaceT::dimension,CSpaceT::dimNC> > tempSpace;

//    if( dimNCF==dimNCC )
      // doesnot work, would be nice feature makes only seens when one would
      // also work on the finest level with second order, then otherwise some
      // uncesary pointers are stored
//      tempSpace = space;
      //
//    else {
      auto tempSpace = createCoarseSpaceT( space );
//    }
    std::vector<Teuchos::RCP<const CSpaceT> > multiSpace( 1, tempSpace );

    const Ordinal* nLoc_ = space->nLoc();
    Ordinal nLoc[dimension];

    for( int i=0; i<dimension; ++i ) {
      nLoc[i] = nLoc_[i];
    }

    bool coarsen_yes;
    bool coarsen_dir[dimension];

    for( int i=1; i<maxGrids; ++i ) {
      coarsen_yes = false;

      for( int j=0; j<dimension; ++j ) {
        coarsen_dir[j] = false;
        if( j<3 ) {
          if( ( (nLoc[j]-1)%2 )==0 && nLoc[j]>=5 ) {
//          if( ( (nLoc[j]-1)%2 )==0 && nLoc[j]>=3 ) {
            nLoc[j] = (nLoc[j]-1)/2 + 1;
            coarsen_yes = true;
            coarsen_dir[j] = true;
          }
        }
        else
          if( ( (nLoc[j])%2 )==0 && nLoc[j]>1 ) {
            nLoc[j] = (nLoc[j])/2;
            coarsen_yes = true;
            coarsen_dir[j] = true;
          }
      }
      if( coarsen_yes )  multiSpace.push_back( createCoarseSpace( multiSpace.back(), coarsen_dir ) );
    }

    multiSpace.shrink_to_fit();

    return( multiSpace );

  }

protected:

  template<class SpaceT>
  static Teuchos::RCP< const SpaceT > createCoarseSpace(
      const Teuchos::RCP<const SpaceT>& space,
      bool coarsen_dir[dimension] ) {

//    Teuchos::RCP<const StencilWidths<dimension,Space::dimNC> > stencilWidths;
    auto stencilWidths = space->getStencilWidths();

    auto domain = space->getDomain();

    auto boundaryConditionsLocal = domain->getBCLocal();

    auto procGridSize = space->getProcGridSize();

    auto gridSizeGlobalTup = space->getGridSizeGlobal()->getTuple();

    for( int i=0; i<dimension; ++i )
      if( coarsen_dir[i] )
        gridSizeGlobalTup[i] = (gridSizeGlobalTup[i]-1)/2 +1;

    // which template paramters are necessary?
    auto gridSizeGlobal = createGridSizeGlobal<Ordinal,dimension>( gridSizeGlobalTup ); //coarsen

    auto gridSizeLocal = Pimpact::createGridSizeLocal<Ordinal,dimension>(
        gridSizeGlobal, procGridSize );

    auto procGrid = Pimpact::createProcGrid<Ordinal,dimension>(
        gridSizeLocal, domain->getBCGlobal(), procGridSize );


    auto indexSpace = Pimpact::createIndexSpace<Ordinal,dimension>(
        stencilWidths, gridSizeLocal, boundaryConditionsLocal, false );

    auto  coordGlobal = Pimpact::createGridCoordinatesGlobal<Scalar,Ordinal,dimension>(
        gridSizeGlobal,
        domain->getDomainSize() );

    auto  coordLocal = Pimpact::createGridCoordinatesLocal<Scalar,Ordinal,dimension>(
        stencilWidths,
        domain->getDomainSize(),
        gridSizeGlobal,
        gridSizeLocal,
        domain->getBCGlobal(),
        domain->getBCLocal(),
        procGrid,
        coordGlobal );

    auto interV2S =
        Pimpact::createInterpolateV2S<Scalar,Ordinal,dimension>(
            procGrid,
            gridSizeLocal,
            stencilWidths,
            domain,
            coordLocal );

    return(
         Teuchos::rcp(
             new SpaceT(
                 stencilWidths,
                 indexSpace,
                 gridSizeGlobal,
                 gridSizeLocal,
                 procGridSize,
                 procGrid,
                 coordGlobal,
                 coordLocal,
                 domain,
                 interV2S )
         )
    );

  }

  static Teuchos::RCP< const CSpaceT > createCoarseSpaceT(
      const Teuchos::RCP<const FSpaceT>& space ) {

    auto stencilWidths = createStencilWidths<dimension,dimNCC>();

    auto domain = space->getDomain();

    auto boundaryConditionsLocal = domain->getBCLocal();

    auto procGridSize = space->getProcGridSize();

    auto gridSizeGlobalTup = space->getGridSizeGlobal()->getTuple();

    auto gridSizeGlobal = space->getGridSizeGlobal();

    auto gridSizeLocal = space->getGridSizeLocal();

    auto procGrid = space->getProcGrid();

    auto indexSpace = space->getIndexSpace();

    auto  coordGlobal = space->getCoordinatesGlobal();

    auto  coordLocal = Pimpact::createGridCoordinatesLocal(
        stencilWidths,
        domain->getDomainSize(),
        gridSizeGlobal,
        gridSizeLocal,
        domain->getBCGlobal(),
        domain->getBCLocal(),
        procGrid,
        coordGlobal );

    auto interV2S =
        Pimpact::createInterpolateV2S(
            procGrid,
            gridSizeLocal,
            stencilWidths,
            domain,
            coordLocal );

    return(
         Teuchos::rcp(
             new CSpaceT(
                 stencilWidths,
                 indexSpace,
                 gridSizeGlobal,
                 gridSizeLocal,
                 procGridSize,
                 procGrid,
                 coordGlobal,
                 coordLocal,
                 domain,
                 interV2S )
         )
    );

  }

}; // end of class CoarsenStrategy



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COARSENSTRATEGY_HPP
