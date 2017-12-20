#pragma once
#ifndef PIMPACT_TRANSFERMODEOP_HPP
#define PIMPACT_TRANSFERMODEOP_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_ModeField.hpp"




namespace Pimpact {



template<class InterT>
class TransferModeOp {

public:

  using FSpaceT = typename InterT::FSpaceT;
  using CSpaceT = typename InterT::CSpaceT;

  using SpaceT = typename InterT::SpaceT;

  using DomainFieldT = ModeField<typename InterT::DomainFieldT>;
  using RangeFieldT = ModeField<typename InterT::RangeFieldT>;

protected:

  using OT = typename SpaceT::Ordinal;

  Teuchos::RCP<InterT> op_;

public:

  TransferModeOp(
    const Teuchos::RCP<const FSpaceT>& spaceC,
    const Teuchos::RCP<const CSpaceT>& spaceF ):
    op_( Teuchos::rcp( new InterT( spaceC, spaceF ) ) ) {}

  TransferModeOp(
    const Teuchos::RCP<const FSpaceT>& spaceC,
    const Teuchos::RCP<const CSpaceT>& spaceF,
    const Teuchos::Tuple<int,SpaceT::dimension>& nb ):
    op_( Teuchos::rcp( new InterT( spaceC, spaceF, nb ) ) ) {}



  template<class DT, class RT>
  void apply( const DT& x, RT& y ) const {

    op_->apply( x.getCField(), y.getCField() );
    op_->apply( x.getSField(), y.getSField() );
  }


  void print( std::ostream& out=std::cout ) const {

    out << "=== TransferModeOP ===\n";
    op_->print( out );
  }


}; // end of class TransferModeOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRANSFEROP_HPP
