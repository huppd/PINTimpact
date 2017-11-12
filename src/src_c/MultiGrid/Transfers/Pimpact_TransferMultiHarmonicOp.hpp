#pragma once
#ifndef PIMPACT_TRANSFERMULTIHARMONICOP_HPP
#define PIMPACT_TRANSFERMULTIHARMONICOP_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_MultiHarmonicField.hpp"




namespace Pimpact {



template<class InterT>
class TransferMultiHarmonicOp {

public:

  using FSpaceT = typename InterT::FSpaceT;
  using CSpaceT = typename InterT::CSpaceT;

  using SpaceT = typename InterT::SpaceT;

  using DomainFieldT = MultiHarmonicField<typename InterT::DomainFieldT>;
  using RangeFieldT = MultiHarmonicField<typename InterT::RangeFieldT>;

protected:

  using OT = typename SpaceT::Ordinal;

  Teuchos::RCP<InterT> op_;

public:

  TransferMultiHarmonicOp(
    const Teuchos::RCP<const FSpaceT>& spaceC,
    const Teuchos::RCP<const CSpaceT>& spaceF ):
    op_( Teuchos::rcp( new InterT( spaceC, spaceF ) ) ) {}

  TransferMultiHarmonicOp(
    const Teuchos::RCP<const FSpaceT>& spaceC,
    const Teuchos::RCP<const CSpaceT>& spaceF,
    const Teuchos::Tuple<int,SpaceT::dimension>& nb ):
    op_( Teuchos::rcp( new InterT( spaceC, spaceF, nb ) ) ) {}



  template<class DT, class RT>
  void apply( const DT& x_ref, RT& y ) const {

    Teuchos::RCP<const DT> x;
    if( x_ref.global() )
      x = Teuchos::rcpFromRef( x_ref );
    else {
      Teuchos::RCP<DT> temp = Teuchos::rcp( new DT( x_ref.space(), true ) );
      *temp = x_ref;
      x = temp;
    }

    x->exchange();
    //std::cout << y.space()->si(F::U,3) << "\n";

    if( 0==y.space()->si(F::U,3) )
      y.get0Field() = x->get0Field();
      //op_->apply( x->get0Field(), y.get0Field() );

    OT iS = std::max( y.space()->si(F::U,3), 1 );
    OT iE = std::min( x->space()->nGlo(3), y.space()->ei(F::U,3) );

    for( OT i=iS; i<=iE; ++i ) {
      //op_->apply( x->getCField(i), y.getCField(i) );
      //op_->apply( x->getSField(i), y.getSField(i) );
      y.getCField(i) = x->getCField(i);
      y.getSField(i) = x->getSField(i);
    }

    iS = std::max(x->space()->nGlo(3)+1, y.space()->si(F::U,3));
    iE = y.space()->ei(F::U,3);
    for( OT i=iS; i<=iE; ++i ) {
      y.getCField(i).init();
      y.getSField(i).init();
    }
    y.changed();
  }


  void print(  std::ostream& out=std::cout ) const {

    out << "=== TransferMultiHarmonicOP ===\n";
    op_->print( out );
  }


}; // end of class TransferMultiHarmonicOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRANSFEROP_HPP
