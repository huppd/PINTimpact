/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_SPACEFACTORY_HPP
#define PIMPACT_SPACEFACTORY_HPP


//#include "Teuchos_Tuple.hpp"
//#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
//#include "Teuchos_ParameterList.hpp"
//#include "Teuchos_XMLParameterListCoreHelpers.hpp"




namespace Pimpact {



/// \brief creates non const template from const
template<template<class> class GridObjectT, class GridT>
Teuchos::RCP<GridObjectT<GridT> >
create(const Teuchos::RCP<const GridT>& grid) {

  return Teuchos::rcp(new GridObjectT<GridT>(grid));
}


/// \brief creates non const template from non const
template<template<class> class GridObjectT, class GridT>
Teuchos::RCP<GridObjectT<GridT> >
create(const Teuchos::RCP<GridT>& grid) {
  return Teuchos::rcp(new GridObjectT<GridT>(grid));
}


/// \brief creates non const from non const
template<class GridObjectT>
Teuchos::RCP<GridObjectT>
create(const Teuchos::RCP<const typename GridObjectT::GridT>& grid) {

  return Teuchos::rcp(new GridObjectT(grid));
}


/// \brief creates non const from non const
template<class GridObjectT, class GridT>
Teuchos::RCP<GridObjectT >
create(const Teuchos::RCP<GridT>& grid) {

  return Teuchos::rcp(new GridObjectT(grid));
}


template<template<class> class GridObjectT, class GridT>
Teuchos::RCP<const GridObjectT<GridT> >
createConst(const Teuchos::RCP<const GridT>& grid) {

  return Teuchos::rcp(new GridObjectT<GridT>(grid));
}


template<template<class> class GridObjectT, class GridT>
Teuchos::RCP<const GridObjectT<GridT> >
createConst(const Teuchos::RCP<GridT>& grid) {

  return Teuchos::rcp(new GridObjectT<GridT>(grid));
}


template<class GridObjectT, class GridT>
Teuchos::RCP<const GridObjectT >
createConst(const Teuchos::RCP<const GridT>& grid) {

  return Teuchos::rcp(new GridObjectT(grid));
}


template<class GridObjectT, class GridT>
Teuchos::RCP<const GridObjectT >
createConst(const Teuchos::RCP<GridT>& grid) {
  return Teuchos::rcp(new GridObjectT(grid));
}


/// \relates TransferOp
template<class OpT>
Teuchos::RCP<const OpT >
create(
  const Teuchos::RCP<const typename OpT::FGridT>& fGrid,
  const Teuchos::RCP<const typename OpT::CGridT>& cGrid) {

  return Teuchos::rcp(new OpT(fGrid, cGrid));
}




} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SPACEFACTORY_HPP
