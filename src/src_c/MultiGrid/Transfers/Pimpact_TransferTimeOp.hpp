#pragma once
#ifndef PIMPACT_TRANSFERTIMEOP_HPP
#define PIMPACT_TRANSFERTIMEOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Teuchos_TestForException.hpp"



namespace Pimpact {


template<class TransVT>
class TransferTimeOp {

public:

	using FSpaceT = typename TransVT::FSpaceT;
	using CSpaceT = typename TransVT::CSpaceT;

	using SpaceT = typename TransVT::SpaceT;

	using DomainFieldT = TimeField<typename TransVT::DomainFieldT>;
	using RangeFieldT = TimeField<typename TransVT::RangeFieldT>;
	
protected:

	Teuchos::RCP<TransVT> opV_;

public:

	TransferTimeOp(
			const Teuchos::RCP<const FSpaceT>& spaceC,
			const Teuchos::RCP<const CSpaceT>& spaceF ):
		opV_( Teuchos::rcp( new TransVT( spaceC, spaceF ) ) ) {}


	template< class DT, class RT>
	void apply( const DT& x, RT& y ) const {

		for (int i=x.space()->begin(F::S,3); i<x.space()->end(F::S,3); ++i)
			opV_->apply( x(i), y(i) );

		y.changed();
	}

	
	void print(  std::ostream& out=std::cout ) const {
		out << "=== TransferTimeOP ===\n";
		opV_->print( out );
	}

}; // end of class TransferTime


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
