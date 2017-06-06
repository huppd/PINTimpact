#pragma once
#ifndef PIMPACT_TRANSFERCOMPOUNDOP_HPP
#define PIMPACT_TRANSFERCOMPOUNDOP_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_CompoundField.hpp"




namespace Pimpact {



template<class TransVT, class TransST>
class TransferCompoundOp {

public:

  using FSpaceT = typename TransVT::FSpaceT;
  using CSpaceT = typename TransVT::CSpaceT;

  using SpaceT = typename TransVT::SpaceT;

  using DomainFieldT = CompoundField<typename TransVT::DomainFieldT, typename TransST::DomainFieldT >;
  using RangeFieldT = CompoundField<typename TransVT::RangeFieldT,  typename TransST::RangeFieldT  >;

protected:

  Teuchos::RCP<TransVT> opV_;
  Teuchos::RCP<TransST> opS_;

public:

  TransferCompoundOp(
    const Teuchos::RCP<const FSpaceT>& spaceC,
    const Teuchos::RCP<const CSpaceT>& spaceF ):
    opV_( Teuchos::rcp( new TransVT( spaceC, spaceF ) ) ),
    opS_( Teuchos::rcp( new TransST( spaceC, spaceF ) ) ) {}

  TransferCompoundOp(
    const Teuchos::RCP<const FSpaceT>& spaceC,
    const Teuchos::RCP<const CSpaceT>& spaceF,
    const Teuchos::Tuple<int,SpaceT::dimension>& nb ):
    opV_( Teuchos::rcp( new TransVT( spaceC, spaceF ) ) ),
    opS_( Teuchos::rcp( new TransST( spaceC, spaceF ) ) ) {}


  template<class DT, class RT>
  void apply( const DT& x, RT& y ) const {

    opV_->apply( x.getVField(), y.getVField() );
    opS_->apply( x.getSField(), y.getSField() );
  }

//	void apply( const RangeFieldT& x, DomainFieldT& y ) const {
//
//    opV_->apply( x.getVField(), y.getVField() );
//    opS_->apply( x.getSField(), y.getSField() );
//
//	}


  void print(  std::ostream& out=std::cout ) const {

    out << "=== TransferCompoundOP ===\n";
    opV_->print( out );
    opS_->print( out );
  }

}; // end of class TransferCompoundOp



} // end of namespace Pimpact



#ifdef COMPILE_ETI
#include "Pimpact_Space.hpp"
#include "Pimpact_InterpolationOp.hpp"
#include "Pimpact_VectorFieldOpWrap.hpp"
//extern template class Pimpact::TransferCompoundOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > >;
//extern template class Pimpact::TransferCompoundOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > >;
//extern template class Pimpact::TransferCompoundOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > > >;
//extern template class Pimpact::TransferCompoundOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > > >;
#endif


#endif // end of #ifndef PIMPACT_TRANSFERCOMPOUNDOP_HPP
