#pragma once
#ifndef PIMPACT_COARSENSTRATEGY_HPP
#define PIMPACT_COARSENSTRATEGY_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {


/// \brief first model implementation, where one coarsens until every processor has 3 dofs.
///
/// it should be that the coardinates are taken from the fine grid and are halved.
/// \todo compute coordinates correctly (grid stretching only on finest grid, afterwards simple coarsening),
/// cleanest version would be to use same grid stretching on every level, makes interpolation and restriction slightly more complicated.
template<class Field>
class CoarsenStrategy {

  typedef typename Field::Scalar  Scalar;
  typedef typename Field::Ordinal Ordinal;

  static const int dimension = Field::dimension;

  typedef typename Field::SpaceT SpaceT;

public:


  static std::vector<Teuchos::RCP<const SpaceT> > getMultiSpace(
      const Teuchos::RCP<const SpaceT> space,
      int maxGrids=10 ) {

    std::vector<Teuchos::RCP<const SpaceT> > multiSpace( 1, space );

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

  static Teuchos::RCP< const SpaceT > createCoarseSpace(
      const Teuchos::RCP<const SpaceT>& space,
      bool coarsen_dir[dimension] ) {

    auto domain = space->getDomain();

    auto boundaryConditionsLocal = domain->getBCLocal();

    auto procGridSize = space->getProcGridSize();

    auto gridSizeGlobalTup = space->getGridSizeGlobal()->getTuple();

    for( int i=0; i<dimension; ++i )
      if( coarsen_dir[i] )
        gridSizeGlobalTup[i] = (gridSizeGlobalTup[i]-1)/2 +1;

    auto gridSizeGlobal = createGridSizeGlobal<Ordinal,dimension>( gridSizeGlobalTup ); //coarsen

    auto gridSizeLocal = Pimpact::createGridSizeLocal<Ordinal,dimension>(
        gridSizeGlobal, procGridSize );

    auto procGrid = Pimpact::createProcGrid<Ordinal,dimension>(
        gridSizeLocal, domain->getBCGlobal(), procGridSize );

    auto fieldSpace = space->getFieldSpace();

    auto scalarIndexSpace = Pimpact::createScalarIndexSpace<Ordinal,dimension>(
        fieldSpace, gridSizeLocal, boundaryConditionsLocal, false );

    auto innerIndexSpace = Pimpact::createInnerFieldIndexSpaces<Ordinal,dimension>(
        fieldSpace, gridSizeLocal, boundaryConditionsLocal, false );

    auto  fullIndexSpace = Pimpact::createFullFieldIndexSpaces<Ordinal,dimension>(
        fieldSpace, gridSizeLocal, boundaryConditionsLocal, false );

    auto  coordGlobal = Pimpact::createGridCoordinatesGlobal<Scalar,Ordinal,dimension>(
        gridSizeGlobal,
        domain->getDomainSize() );

    auto  coordLocal = Pimpact::createGridCoordinatesLocal<Scalar,Ordinal,dimension>(
        fieldSpace,
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
            fieldSpace,
            domain,
            coordLocal );

    return(
         Teuchos::rcp(
             new SpaceT(
                 fieldSpace,
                 scalarIndexSpace,
                 innerIndexSpace,
                 fullIndexSpace,
                 gridSizeGlobal,
                 gridSizeLocal,
                 procGridSize,
                 procGrid,
                 coordGlobal,
                 coordLocal,
                 domain,
                 interV2S ) ) );


  }

}; // end of class CoarsenStrategy



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COARSENSTRATEGY_HPP
