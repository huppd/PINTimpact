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

	using Ordinal = typename SpaceT::Ordinal;

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
	void apply( const DT& x, RT& y ) const {

		x.exchange();

		if( y.space()->begin(U,3)<0 )
			op_->apply( x.getConst0Field(), y.get0Field() );

		//    int m = std::min( x.space()->nGlo(3), y.space()->nGlo(3) );
		//    for( int i=0; i<m; ++i ) {

		Ordinal iS = std::max( y.space()->begin(U,3), 0 );
		Ordinal iE = std::min( x.space()->nGlo(3), y.space()->end(U,3) );

		for( Ordinal i=iS; i<iE; ++i ) {
			op_->apply( x.getConstCField(i), y.getCField(i) );
			op_->apply( x.getConstSField(i), y.getSField(i) );
		}

		iS = std::max(x.space()->nGlo(3),y.space()->begin(U,3));
		iE = y.space()->end(U,3);
		for( Ordinal i=iS; i<iE; ++i ) {
			y.getCField(i).initField();
			y.getSField(i).initField();
		}

		y.changed();

	}


  void print(  std::ostream& out=std::cout ) const {

		out << "=== TransferMultiHarmonicOP ===\n";
		op_->print( out );
  }


}; // end of class TransferMultiHarmonicOp



} // end of namespace Pimpact



#ifdef COMPILE_ETI
#include "Pimpact_Space.hpp"
#include "Pimpact_InterpolationOp.hpp"
#include "Pimpact_VectorFieldOpWrap.hpp"
extern template class Pimpact::TransferMultiHarmonicOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::TransferMultiHarmonicOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > > >;
extern template class Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > > >;
#endif


#endif // end of #ifndef PIMPACT_TRANSFEROP_HPP
