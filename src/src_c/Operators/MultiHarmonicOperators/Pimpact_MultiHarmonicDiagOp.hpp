#pragma once
#ifndef PIMPACT_MULDTIHARMONICDIAGOP_HPP
#define PIMPACT_MULDTIHARMONICDIAGOP_HPP


#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_Utils.hpp"



namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template< class ZeroOpT, class ModeOpT >
class MultiHarmonicDiagOp {

public:

  using SpaceT = typename ZeroOpT::SpaceT;

  using DomainFieldT = MultiHarmonicField<typename ZeroOpT::DomainFieldT>;
  using RangeFieldT = MultiHarmonicField<typename ZeroOpT::RangeFieldT>;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  Teuchos::RCP<ZeroOpT> zeroOp_;
  Teuchos::RCP<ModeOpT> modeOp_;

public:

  MultiHarmonicDiagOp( const Teuchos::RCP<ZeroOpT>& zeroOp, const Teuchos::RCP<ModeOpT>& modeOp ):
    zeroOp_( zeroOp ),
    modeOp_( modeOp ) {};


  /// \todo either tangle this with wind of operator or make an exchange ZeroField (tangle preferred)
  void assignField( const DomainFieldT& y_ref ) {

    Teuchos::RCP<const DomainFieldT> y;

    if( y_ref.global()==DomainFieldT::Global::Y )
      y = Teuchos::rcpFromRef( y_ref );
    else {
      Teuchos::RCP<DomainFieldT> temp = Teuchos::rcp( new DomainFieldT( space(), DomainFieldT::Global::Y ) );
      *temp = y_ref;
      y = temp;
      //std::cout << "prec: y->global(): " << y->global() << "\n";
    }

    y->exchange();
    zeroOp_->assignField( y->get0Field() );
    //modeOp->assignField( y->get0Field() );
  };


  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    Scalar iRe = 1./space()->getDomainSize()->getRe();
    Scalar a2 = space()->getDomainSize()->getAlpha2()*iRe;

    Scalar mulI;

    // computing zero mode of z
    // set paramteters
    if( 0==space()->si(F::U,3) ) {

      // set parameters
      //
      auto para = Teuchos::parameterList();
      para->set<Scalar>( "mulI", 0. );
      para->set<Scalar>( "mulC", 1. );
      para->set<Scalar>( "mulL", iRe );
      zeroOp_->setParameter( para );

      zeroOp_->apply( x.get0Field(), y.get0Field() );
    }

    for( Ordinal i=std::max(space()->si(F::U,3),1); i<=space()->ei(F::U,3); ++i ) {

      // set parameters
      auto para = Teuchos::parameterList();
      para->set<Scalar>( "mulI", a2*i );
      para->set<Scalar>( "mulC", 1. );
      para->set<Scalar>( "mulL", iRe );
      modeOp_->setParameter( para );

      modeOp_->apply( x.getField(i), y.getField(i) );
    }
    y.changed();
  }

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return zeroOp_->space();
  };

  void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {
    zeroOp_->setParameter( para );
    modeOp_->setParameter( para );
  }

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "MultiHarmonicDiagOp ";
  };

  void print( std::ostream& out=std::cout ) const {
    out << getLabel() << ":\n";
    zeroOp_->print( out );
    modeOp_->print( out );
  }

}; // end of class MultiHarmonicDiagOp



/// \relates MultiHarmonicDiagOp
template<class ZeroOpT, class ModeOpT>
Teuchos::RCP<MultiHarmonicDiagOp<ZeroOpT,ModeOpT> >
createMultiHarmonicDiagOp(
  const Teuchos::RCP<ZeroOpT>& zeroOp,
  const Teuchos::RCP<ModeOpT>& modeOp ) {

  return Teuchos::rcp( new MultiHarmonicDiagOp<ZeroOpT,ModeOpT>( zeroOp, modeOp ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULDTIHARMONICDIAGOP_HPP
