#pragma once
#ifndef PIMPACT_MGSPACES_HPP
#define PIMPACT_MGSPACES_HPP


#include <vector>

#include "Pimpact_Grid.hpp"




namespace Pimpact {


/// \brief contains hierarchy of Grids
/// \ingroup MG
/// \tparam FGridT \c Grid type on finest level
/// \tparam CGridT \c Grid type on coarser levels allows to have tighter boundaries, less memory used, less comunication
template<class FST, class CGT>
class MGGrids {

public:

  using FGridT = FST;
  using CGridT = CGT;

protected:

  template<class CoarsenStrategyT>
  friend Teuchos::RCP<const MGGrids<typename CoarsenStrategyT::GridT,
         typename CoarsenStrategyT::CGridT> >
         createMGGrids(const Teuchos::RCP<const typename CoarsenStrategyT::GridT>&
                         grid, int nGridsMax);

  Teuchos::RCP<const FGridT> grid_;
  std::vector<Teuchos::RCP<const CGridT> > grids_;


  MGGrids(
    const Teuchos::RCP<const FGridT>& grid,
    const std::vector<Teuchos::RCP<const CGridT> >& grids):
    grid_(grid),
    grids_(grids) {}

public:

  /// \brief get number of grids
  ///
  /// \return  number of grids
  constexpr int getNGrids() const {
    return grids_.size();
  }

  constexpr const Teuchos::RCP<const FGridT>&  get()        const {
    return grid_;
  }

  /// \brief gets ith grid, similar to python i=-1 is gets you the coarses grid
  ///
  /// \param i index of grid level if negative it is counted from coarsest grid
  ///
  /// \return ith grid
  constexpr const Teuchos::RCP<const CGridT>&  get(int i) const {
    if(i<0)
      return grids_[ getNGrids()+i ];
    else
      return grids_[i];
  }


  constexpr bool participating(int i) const {
    return get(i)->getProcGrid()->participating();
  }

  void print(std::ostream& out=std::cout) const {

    for(int i=0; i<getNGrids(); ++i) {
      out << "-------------------------\n";
      out << "-------- grid : "<< i << "--------\n";
      out << "-------------------------\n";
      get(i)->print(out);
    }
  }

}; // end of class MGGrids



/// \relates MGGrids
template<class CoarsenStrategy>
Teuchos::RCP<const MGGrids<typename CoarsenStrategy::GridT, typename CoarsenStrategy::CGridT> >
createMGGrids(
  const Teuchos::RCP<const typename CoarsenStrategy::GridT>& grid,
  int maxGrids=10) {

  std::vector<Teuchos::RCP<const typename CoarsenStrategy::CGridT> > grids =
    CoarsenStrategy::getMultiGrid(grid, maxGrids);

  return Teuchos::rcp(new MGGrids<typename CoarsenStrategy::GridT, typename
      CoarsenStrategy::CGridT>(grid, grids));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MGSPACES_HPP
