#pragma once
#ifndef PIMPACT_COARSENSTRATEGY_HPP
#define PIMPACT_COARSENSTRATEGY_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {



/// \brief first model implementation, where one coarsens until every processor has 3 dofs.
///
/// \tparam FSpaceT space on finest level not necessary the same for the coarse grids(difference in \c dim_nc).
/// it should be that the coardinates are taken from the fine grid and are halved.
/// \note compute coordinates correctly (grid stretching only on finest grid, afterwards simple coarsening),
/// cleanest version would be to use same grid stretching on every level, makes
/// interpolation and restriction slightly more complicated.
///\note the \c CSpaceT parameter allows even smaller Grids nLoc=1 but therefore
///the exception in GridSilzeLocal has to be adapted on StencilWidths
/// \ingroup MG
template<class FSpaceT, class CST>
class CoarsenStrategy {

public:

  using SpaceT = FSpaceT;
  using CSpaceT = CST;

protected:

  using Scalar = typename FSpaceT::Scalar;
  using Ordinal = typename FSpaceT::Ordinal;

  /// should be same for finest space and coarse spaces
  static const int sdim = FSpaceT::sdim;
  static const int dimension = FSpaceT::dimension;

  using TO = typename Teuchos::Tuple<Ordinal, dimension>;

  /// can be different for finest and coarse spaces
  static const int dimNCF = FSpaceT::dimNC;

  static const int dimNCC = CSpaceT::dimNC;

public:

  static std::vector<Teuchos::RCP<const CSpaceT> > getMultiSpace(
    const Teuchos::RCP<const FSpaceT> space,
    int maxGrids=10) {

    // creating low order space
    Teuchos::RCP<const CSpaceT> tempSpace = createSpace<CSpaceT, FSpaceT>(space);

    std::vector<Teuchos::RCP<const CSpaceT> > multiSpace(1, tempSpace);

    Teuchos::Tuple<Ordinal, dimension> nLoc = *space->getGridSizeLocal();
    GridSizeGlobal<Ordinal, sdim> nGlo = *space->getGridSizeGlobal();

    for(Ordinal i=1; i<maxGrids; ++i) {
      bool coarsen_yes = false;
      for(Ordinal j=0; j<dimension; ++j) {
        if(j<3) {
          if(((nLoc[j]-1)%2)==0 && ((nLoc[j]-1)/2 + 1)%2!=0 && (nLoc[j]-1)/2 + 1>3) {
            nLoc[j] = (nLoc[j]-1)/2 + 1;
            nGlo[j] = (nGlo[j]-1)/2 + 1;
            coarsen_yes = true;
          }
        } else if(!space->getStencilWidths()->spectralT() && ((nLoc[j])%2)==0 && nLoc[j]>1) {
          nLoc[j] = (nLoc[j])/2;
          nGlo[j] = (nGlo[j])/2;
          coarsen_yes = true;
        }
      }
      if(coarsen_yes)
        multiSpace.push_back(createSpace(multiSpace.back(), nGlo));
    }


    // not working on brutus
    //multiSpace.shrink_to_fit();

    return multiSpace;
  }


}; // end of class CoarsenStrategy


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COARSENSTRATEGY_HPP
