#pragma once
#ifndef PIMPACT_INTRESCOMPOUNDOP_HPP
#define PIMPACT_INTRESCOMPOUNDOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Teuchos_TestForException.hpp"

#include "Pimpact_CompoundField.hpp"




namespace Pimpact {



template<class TransVT, class TransST>
class IntResCompoundOp {

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

  IntResCompoundOp(
    const Teuchos::RCP<const FSpaceT>& spaceC,
    const Teuchos::RCP<const CSpaceT>& spaceF ):
    opV_( Teuchos::rcp( new TransVT( spaceC, spaceF ) ) ),
    opS_( Teuchos::rcp( new TransST( spaceC, spaceF ) ) ) {}

  IntResCompoundOp(
    const Teuchos::RCP<const FSpaceT>& spaceC,
    const Teuchos::RCP<const CSpaceT>& spaceF,
    const Teuchos::Tuple<int,SpaceT::dimension>& nb ):
    opV_( Teuchos::rcp( new TransVT( spaceC, spaceF ) ) ),
    opS_( Teuchos::rcp( new TransST( spaceC, spaceF ) ) ) {}


  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    opV_->apply( x.getVField(), y.getVField() );
    opS_->apply( x.getSField(), y.getSField() );
  }


  void print( std::ostream& out=std::cout ) const {

    out << "=== IntResCompoundOP ===\n";
    opV_->print( out );
    opS_->print( out );
  }

}; // end of class IntResCompoundOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTRESCOMPOUNDOP_HPP
