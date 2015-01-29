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


  template< class SP1T, class SP2T>
  void apply( const ScalarField<SP1T>& x, ScalarField<SP2T>& y ) const {

    auto fType = x.getType();

    TEUCHOS_TEST_FOR_EXCEPT( fType != y.getType() );

    for( int i=0; i<x.space()->dim(); ++i )
      TEUCHOS_TEST_FOR_EXCEPT( x.space()->nLoc(i) != y.space()->nLoc(i) );

    OP_Transfer(
        x.space()->nLoc(),
        x.space()->bl(),
        x.space()->bu(),
        x.space()->sInd(fType),
        x.space()->eInd(fType),
        y.space()->bl(),
        y.space()->bu(),
        y.space()->sInd(fType),
        y.space()->eInd(fType),
        x.getConstRawPtr(),
        y.getRawPtr() );

    y.changed();

  }

  void assignField( const RangeFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }


}; // end of class TransferOp



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRANSFEROP_HPP
