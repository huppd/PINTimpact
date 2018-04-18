/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_MGTRANSFERS_HPP
#define PIMPACT_MGTRANSFERS_HPP


#include "Pimpact_InterpolationOp.hpp"
#include "Pimpact_MGFields.hpp"
#include "Pimpact_MGGrids.hpp"
#include "Pimpact_RestrictionSFOp.hpp"
#include "Pimpact_RestrictionVFOp.hpp"
#include "Pimpact_TransferOp.hpp"
#include "Pimpact_VectorFieldOpWrap.hpp"



namespace Pimpact {



/// \ingroup MG
template<class MGGridsT,
         template<class, class> class TransT,
         template<class> class RestrT,
         template<class> class InterT >
class MGTransfers {

public:

  using FGridT = typename MGGridsT::FGridT;
  using CGridT = typename MGGridsT::CGridT;

  using Scalar = typename FGridT::Scalar;
  using Ordinal = typename FGridT::Ordinal;

  static const int dimension = FGridT::dimension;

  static const int dimNCF = FGridT::dimNC;
  static const int dimNCC = CGridT::dimNC;


  using TransferOpT = TransT<FGridT, CGridT>;
  using RestrictionOpT = RestrT<CGridT>;
  using InterpolationOpT = InterT<CGridT>;

//  template<
//    class MGGridsTT,
//    template<class, class> class TransTT,
//    template<class> class RestrTT,
//    template<class> class InterTT >
//  friend
//  Teuchos::RCP<const MGTransfers<MGGridsTT, TransTT, RestrTT, InterTT> >
//  createMGTransfers(
//      const Teuchos::RCP<const MGGridsTT>& grid);

protected:

  Teuchos::RCP<const MGGridsT> mgGrids_;

  Teuchos::RCP<const TransferOpT> transferOp_;

  std::vector<Teuchos::RCP<const RestrictionOpT> >   restrictionOps_;
  std::vector<Teuchos::RCP<const InterpolationOpT> > interpolationOps_;

public:

  MGTransfers(const Teuchos::RCP<const MGGridsT>& mgGrids):
    mgGrids_(mgGrids),
    transferOp_(create<TransferOpT>(mgGrids_->get(), mgGrids_->get(0))),
    restrictionOps_(),
    interpolationOps_() {

    for(unsigned i=0; i <mgGrids_->getNGrids()-1; ++i) {
      restrictionOps_.push_back(
        Teuchos::rcp(
          new RestrT<CGridT>(
            mgGrids_->get(i),
            mgGrids_->get(i+1),
            mgGrids_->get()->getProcGrid()->getNP()
         )
       )
     );
      interpolationOps_.push_back(
        Teuchos::rcp(
          new InterT<CGridT>(
            mgGrids_->get(i+1),
            mgGrids_->get(i),
            mgGrids_->get()->getProcGrid()->getNP()
         )
       )
     );
    }
    // not working on brutus(intel)
    //interpolationOps_.shrink_to_fit();
  }

public:

  constexpr const Teuchos::RCP<const TransferOpT>&       getTransferOp     (     ) const {
    return transferOp_;
  }

  /// \brief gets ith RestrictionOp, similar to python i=-1 is gets you the coarses grid
  constexpr const Teuchos::RCP<const RestrictionOpT>&    getRestrictionOp  (int i) const {
    if(i<0)
      return restrictionOps_[mgGrids_->getNGrids()+i];
    else
      return restrictionOps_[i];
  }

  /// \brief gets ith InterpolationOp, similar to python i=-1 is gets you the coarses grid
  constexpr const Teuchos::RCP<const InterpolationOpT>&  getInterpolationOp(int i) const {
    if(i<0)
      return interpolationOps_[mgGrids_->getNGrids()+i];
    else
      return interpolationOps_[i];
  }

  void print(std::ostream& out=std::cout) const {

    transferOp_->print(out);

    for(int i = 0; i<restrictionOps_.size(); ++i) {
      if(mgGrids_->participating(i)) {
        out << "-------- restrictor: "<< i << "--------\n";
        restrictionOps_[i]->print(out);
      }
    }
    for(int i = 0; i<interpolationOps_.size(); ++i) {
      if(mgGrids_->participating(i)) {
        out << "-------- interpolator: "<< i << "--------\n";
        interpolationOps_[i]->print(out);
      }
    }
  }

  template<template<class> class FieldT>
  void interpolation(MGFields<MGGridsT, FieldT>& x) const {

    for(int i=-2; i>=-mgGrids_->getNGrids(); --i)
      if(mgGrids_->participating(i))
        getInterpolationOp(i)->apply(x.get(i+1), x.get(i));

    getTransferOp()->apply(x.get(0), x.get());
  }

  template<template<class> class FieldT>
  void restriction(MGFields<MGGridsT, FieldT>& x) const {

    getTransferOp()->apply(x.get(), x.get(0));
    for(int i=0; i<mgGrids_->getNGrids()-1; ++i)
      if(mgGrids_->participating(i))
        getRestrictionOp(i)->apply(x.get(i), x.get(i+1));
  }

}; // end of class MGTransfers



/// \relates MGTransfers
template<
  template<class, class> class TransT,
  template<class> class RestrT,
  template<class> class InterT,
  class MGGridsT >
Teuchos::RCP<const MGTransfers<MGGridsT, TransT, RestrT, InterT> >
createMGTransfers(
  const Teuchos::RCP<const MGGridsT>& mgGrids) {

  return Teuchos::rcp(new MGTransfers<MGGridsT, TransT, RestrT, InterT>(mgGrids));
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MGTRANSFERS_HPP
