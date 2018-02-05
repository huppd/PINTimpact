#pragma once
#ifndef PIMPACT_MGFIELDS_HPP
#define PIMPACT_MGFIELDS_HPP


#include <vector>

#include "Teuchos_RCP.hpp"

#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \ingroup MG
template<class MGGridsT, template<class> class FieldT>
class MGFields {

public:

  using FGridT = typename MGGridsT::FGridT;
  using CGridT = typename MGGridsT::CGridT;

  using FFieldT = FieldT<FGridT>;
  using CFieldT = FieldT<CGridT>;

protected:

//  template<template<class> class FieldTT, class MGGridsTT >
//  friend
//  Teuchos::RCP<MGFields<MGGridsTT, FieldTT> >
//  createMGFields(
//      const Teuchos::RCP<const MGGridsTT>& grid,
//      F type=F::S);

  Teuchos::RCP<const MGGridsT> mgGrids_;

  FFieldT              fField_;
  std::vector<CFieldT> cFields_;

public:

  MGFields(const Teuchos::RCP<const MGGridsT>& mgGrids):
    mgGrids_(mgGrids),
    fField_(mgGrids_->get()),
    cFields_() {

    for(int i=0; i<mgGrids_->getNGrids(); ++i)
      //			if(0==i || mgGrids_->participating(i-1))
      cFields_.push_back(CFieldT(mgGrids_->get(i)));

    // not working on brutus
    //cFields_.shrink_to_fit();
  }

public:

  constexpr const FFieldT& get() {
    return fField_;
  }
  FFieldT& get() {
    return fField_;
  }

  /// \brief gets ith operator, similar to python i=-1 is gets you the coarses grid
  constexpr const CFieldT& get(int i) const {
    if(i<0)
      return cFields_[mgGrids_->getNGrids()+i];
    else
      return cFields_[i];
  }

  /// \brief gets ith operator, similar to python i=-1 is gets you the coarses grid
  CFieldT& get(int i)  {
    if(i<0)
      return cFields_[mgGrids_->getNGrids()+i];
    else
      return cFields_[i];
  }


  //  void print(std::ostream& out=std::cout) const {
  //
  //  }

}; // end of class MGFields



/// \relates MGFields
template<template<class> class FieldT, class MGGridsT >
Teuchos::RCP<MGFields<MGGridsT, FieldT> >
createMGFields(const Teuchos::RCP<const MGGridsT>& mgGrids) {

  return Teuchos::rcp(new MGFields<MGGridsT, FieldT>(mgGrids));
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MGFIELDS_HPP
