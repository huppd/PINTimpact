#pragma once
#ifndef PIMPACT_MGOPERATORS_HPP
#define PIMPACT_MGOPERATORS_HPP


#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \ingroup MG
/// \note todo add constructror from FOperator
/// \note todo generalize MGContainer, specify to Fields and Operators
template<class MGGT, template<class> class FOT, template<class> class COT>
class MGOperators {

public:

  using MGGridsT = MGGT;

  using FGridT = typename MGGridsT::FGridT;
  using CGridT = typename MGGridsT::CGridT;

  using FOperatorT = FOT<FGridT>;
  using COperatorT = COT<CGridT>;

protected:

  // why?
  //  template<template<class> class FOTT, template<class> class COTT, class MGGridsTT >
  //  friend
  //  Teuchos::RCP<const MGOperators<MGGridsTT, FOTT, COTT> >
  //  createMGOperators(const Teuchos::RCP<const MGGridsTT>& grid);
  //  template<template<class> class FOTT, template<class> class COTT, class MGGTT >
  //  friend Teuchos::RCP<const MGOperators<MGGTT, FOTT, COTT> >
  //  createMGOperators(
  //      const Teuchos::RCP<const MGGTT>& mgGrids);

  Teuchos::RCP<const MGGridsT> mgGrids_;

  Teuchos::RCP<FOperatorT>                fOperator_;
  std::vector<Teuchos::RCP<COperatorT> >  cOperator_;

public:

  MGOperators(const Teuchos::RCP<const MGGridsT>& mgGrids,
      const Teuchos::RCP<FOperatorT>& fOperator):
    mgGrids_(mgGrids),
    fOperator_(fOperator),
    cOperator_(mgGrids_->getNGrids()) {

    for(int i=0; i<mgGrids_->getNGrids(); ++i)
      if(mgGrids_->participating(i))
        cOperator_[i] = Teuchos::rcp(new COperatorT(mgGrids_->get(i)));

    // not working on brutus
    //cOperator_.shrink_to_fit();
  }

  constexpr const Teuchos::RCP<const MGGridsT>& getMGGrids() const {
    return mgGrids_;
  }

  constexpr const Teuchos::RCP<FOperatorT>& get() const {
    return fOperator_;
  }

  /// \brief gets ith operator, similar to python i=-1 is gets you the coarses grid
  constexpr const Teuchos::RCP<COperatorT>& get(int i) const {
    if(i<0)
      return cOperator_[mgGrids_->getNGrids()+i];
    else
      return cOperator_[i];
  }

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) const {

    fOperator_->setParameter(para);

    for(int i=0; i<mgGrids_->getNGrids(); ++i)
      if(mgGrids_->participating(i))
        cOperator_[i]->setParameter(para);
  }

  void print(std::ostream& out=std::cout) const {

    fOperator_->print();
    for(int i=0; i<mgGrids_->getNGrids(); ++i) {
      std::cout << "\n\n --- level: " << i << " ---\n\n";
      get(i)->print();
    }
  }

}; // end of class MGOperators



/// \relates MGOperators
template<template<class> class FOperatorT, template<class> class COperatorT=FOperatorT, class MGGridsT >
Teuchos::RCP<const MGOperators<MGGridsT, FOperatorT, COperatorT> >
createMGOperators(
  const Teuchos::RCP<const MGGridsT>& mgGrids,
  const Teuchos::RCP<FOperatorT<typename MGGridsT::FGridT> >& fOperator) {

  return Teuchos::rcp(
      new MGOperators<MGGridsT, FOperatorT, COperatorT>(mgGrids, fOperator));
}



} // end of namespace Pimpact



#include "Pimpact_MGGrids.hpp"
#include "Pimpact_DivGradOp.hpp"
#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_NonlinearOp.hpp"
#include "Pimpact_ConvectionDiffusionSOp.hpp"
template<class T> using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;


#endif // end of #ifndef PIMPACT_MGOPERATORS_HPP
