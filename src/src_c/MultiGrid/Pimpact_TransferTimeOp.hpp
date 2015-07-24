#pragma once
#ifndef PIMPACT_TRANSFERTIMEOP_HPP
#define PIMPACT_TRANSFERTIMEOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Teuchos_TestForException.hpp"

//#include "Pimpact_CompoundField.hpp"


namespace Pimpact {


template<class TransVT>
class TransferTimeOp {

public:

	typedef typename TransVT::FSpaceT FSpaceT;
	typedef typename TransVT::CSpaceT CSpaceT;

	typedef typename TransVT::SpaceT SpaceT;
  
	typedef TimeField<typename TransVT::DomainFieldT> DomainFieldT;
	typedef TimeField<typename TransVT::RangeFieldT> RangeFieldT;
	
protected:

	Teuchos::RCP<TransVT> opV_;
public:

	TransferTimeOp(
			const Teuchos::RCP<const SpaceT>& spaceC,
			const Teuchos::RCP<const SpaceT>& spaceF ):
		opV_( Teuchos::rcp( new TransVT( spaceC, spaceF ) ) ) {}


	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		for (int i=x.space()->sInd(S,3); i<x.space()->eInd(S,3); ++i)
	   		 opV_->apply( x.getConstVField(i), y.getVField(i) );

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
