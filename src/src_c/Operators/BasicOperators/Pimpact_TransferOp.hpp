#pragma once
#ifndef PIMPACT_TRANSFEROP_HPP
#define PIMPACT_TRANSFEROP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"




namespace Pimpact{


extern "C" {

void OP_Transfer(
    const int* const N,
    const int* const bLI,
    const int* const bUI,
    const int* const SI,
    const int* const NI,
    const int* const bLO,
    const int* const bUO,
    const int* const SO,
    const int* const NO,
    const double* const phiIN,
    double* const phiOUT );

}


/// \brief Transfers fields from "coarse" to "fine" spaces, necessary when \c Space::dimNC is  different.
///
/// Goes in both direction. If this is used a lot, it could be beneficial, to
/// seperate StencilWidths with data layout, so having same datalayout for all
/// stencile. Coping could be beneficial because Cash effects are bether
///
/// \tparam FSpaceT fine space type in the sense of stencil order.
/// \tparam CSpaceT coase space type
/// \ingroup BaseOperator
/// \todo add test if nloc is equal
template<class FST, class CST>
class TransferOp {

public:

  typedef FST SpaceT;

  typedef FST FSpaceT;
  typedef CST CSpaceT;

  typedef typename FSpaceT::Scalar Scalar;
  typedef typename FSpaceT::Ordinal Ordinal;

  typedef ScalarField<FSpaceT>  DomainFieldT;
  typedef ScalarField<CSpaceT>  RangeFieldT;

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  Teuchos::RCP<const FSpaceT> fSpace_;
  Teuchos::RCP<const CSpaceT> cSpace_;

public:

  TransferOp(
      const Teuchos::RCP<const FSpaceT>& fSpace,
      const Teuchos::RCP<const CSpaceT>& cSpace ):
        fSpace_(fSpace),cSpace_(cSpace) {}


  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    auto fType = x.fType_;

    TEUCHOS_TEST_FOR_EXCEPTION(
        fType != y.fType_,
        std::logic_error,
        "Pimpact::TransferOp!!!\n");

    OP_Transfer(
        fSpace_->nLoc(),
        fSpace_->bl(),
        fSpace_->bu(),
        fSpace_->sInd(fType),
        fSpace_->eInd(fType),
        cSpace_->bl(),
        cSpace_->bu(),
        cSpace_->sInd(fType),
        cSpace_->eInd(fType),
        x.s_,
        y.s_ );

    y.changed();
  }

  void apply(const RangeFieldT& x, DomainFieldT& y) const {
    auto fType = x.fType_;

    TEUCHOS_TEST_FOR_EXCEPTION(
        fType != y.fType_,
        std::logic_error,
        "Pimpact::TransferOp!!!\n");

    OP_Transfer(
        cSpace_->nLoc(),
        cSpace_->bl(),
        cSpace_->bu(),
        cSpace_->sInd(fType),
        cSpace_->eInd(fType),
        fSpace_->bl(),
        fSpace_->bu(),
        fSpace_->sInd(fType),
        fSpace_->eInd(fType),
        x.s_,
        y.s_ );

    y.changed();
  }

  void assignField( const RangeFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }


}; // end of class TransferOp



/// \relates TransferOp
template<class FSpaceT, class CSpaceT>
Teuchos::RCP<const TransferOp<FSpaceT,CSpaceT> >
createTransferOp(
    const Teuchos::RCP<const FSpaceT>& fSpace,
    const Teuchos::RCP<const CSpaceT>& cSpace ) {

  return(
      Teuchos::rcp(
          new TransferOp<FSpaceT,CSpaceT>(fSpace,cSpace)
      )
  );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRANSFEROP_HPP
