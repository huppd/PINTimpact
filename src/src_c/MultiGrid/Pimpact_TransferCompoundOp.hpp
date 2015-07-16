#pragma once
#ifndef PIMPACT_TRANSFERCOMPOUNDOP_HPP
#define PIMPACT_TRANSFERCOMPOUNDOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Teuchos_TestForException.hpp"

#include "Pimpact_CompoundField.hpp"




namespace Pimpact {



template<class TransVT, class TransST>
class TransferCompoundOp {

public:

	typedef typename TransVT::FSpaceT FSpaceT;
	typedef typename TransVT::CSpaceT CSpaceT;

	typedef typename TransVT::SpaceT SpaceT;
  
	typedef CompoundField<typename TransVT::DomainFieldT, typename TransST::DomainFieldT > DomainFieldT;
	typedef CompoundField<typename TransVT::RangeFieldT,  typename TransST::RangeFieldT  > RangeFieldT;
	
protected:

	Teuchos::RCP<TransVT> opV_;
	Teuchos::RCP<TransST> opS_;

public:

	TransferCompoundOp(
			const Teuchos::RCP<const SpaceT>& spaceC,
			const Teuchos::RCP<const SpaceT>& spaceF ):
		opV_( Teuchos::rcp( new TransVT( spaceC, spaceF ) ) ),
		opS_( Teuchos::rcp( new TransST( spaceC, spaceF ) ) ) {}

	TransferCompoundOp(
			const Teuchos::RCP<const SpaceT>& spaceC,
			const Teuchos::RCP<const SpaceT>& spaceF,
			const Teuchos::Tuple<int,SpaceT::dimension>& nb ):
		opV_( Teuchos::rcp( new TransVT( spaceC, spaceF ) ) ),
		opS_( Teuchos::rcp( new TransST( spaceC, spaceF ) ) ) {}


	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    opV_->apply( x.getConstVField(), y.getVField() );
    opS_->apply( x.getConstSField(), y.getSField() );

	}


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
