#pragma once
#ifndef PIMPACT_MULDTIHARMONICDIAGOP_HPP
#define PIMPACT_MULDTIHARMONICDIAGOP_HPP


#include "Pimpact_EddyPrec.hpp"
#include "Pimpact_FieldFactory.hpp"
#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_Types.hpp"



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

  /// \todo get nf from grid
  MultiHarmonicDiagOp( const Teuchos::RCP<ZeroOpT>& zeroOp, const Teuchos::RCP<ModeOpT>& modeOp ):
		zeroOp_( zeroOp ),
		modeOp_( modeOp ) {};


  void assignField( const DomainFieldT& mv ) {
    zeroOp_->assignField( mv.getConst0Field() );
  };


  void apply( const DomainFieldT& x, RangeFieldT& y ) const {
		
		Scalar iRe = 1./zeroOp_->space()->getDomainSize()->getRe();
		Scalar a2 = zeroOp_->space()->getDomainSize()->getAlpha2()*iRe;

		Scalar mulI;

    // computing zero mode of z
		// set paramteters
		auto para = Teuchos::parameterList();
		para->set<Scalar>( "mulI", 0. );
		para->set<Scalar>( "mulC", 1. );
		para->set<Scalar>( "mulL", iRe );
		zeroOp_->setParameter( para );

		if( 0==space()->sInd(U,3) )
			zeroOp_->apply( x.getConst0Field(), y.get0Field() );

		for( Ordinal i=std::max(space()->sInd(U,3),1); i<=space()->eInd(U,3); ++i ) {
			// set parameters
			para->set<Scalar>( "mulI", a2*i );
			para->set<Scalar>( "mulC", 1. );
			para->set<Scalar>( "mulL", iRe );
			modeOp_->setParameter( para );
			modeOp_->apply( x.getConstField(i), y.getField(i) );
		}
		y.changed();
	}

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(zeroOp_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "MultiHarmonicDiagOp " ); };

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

  return( Teuchos::rcp( new MultiHarmonicDiagOp<ZeroOpT,ModeOpT>( zeroOp, modeOp ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULDTIHARMONICDIAGOP_HPP
