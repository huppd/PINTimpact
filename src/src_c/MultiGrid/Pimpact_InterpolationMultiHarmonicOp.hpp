#pragma once
#ifndef PIMPACT_INTERPOLATIONMULTIHARMONICOP_HPP
#define PIMPACT_INTERPOLATIONMULTIHARMONICOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Teuchos_TestForException.hpp"

#include "Pimpact_Space.hpp"
#include "Pimpact_MultiHarmonicField.hpp"




namespace Pimpact {



template<class InterT>
class InterpolationMultiHarmonicOp {

public:

	typedef typename InterT::FSpaceT FSpaceT;
	typedef typename InterT::CSpaceT CSpaceT;

	typedef typename InterT::SpaceT SpaceT;
  
	typedef MultiHarmonicField<typename InterT::DomainFieldT> DomainFieldT;
	typedef MultiHarmonicField<typename InterT::RangeFieldT> RangeFieldT;
	
protected:

//	typedef typename SpaceT::Scalar Scalar;
//  typedef typename SpaceT::Ordinal Ordinal;

//  typedef ScalarField<SpaceT>  DomainFieldT;
//  typedef ScalarField<SpaceT>  RangeFieldT;

	Teuchos::RCP<InterT> op_;

public:

	InterpolationMultiHarmonicOp(
			const Teuchos::RCP<const SpaceT>& spaceC,
			const Teuchos::RCP<const SpaceT>& spaceF ):
		op_( Teuchos::rcp( new InterT( spaceC, spaceF ) ) ) {}

	InterpolationMultiHarmonicOp(
			const Teuchos::RCP<const SpaceT>& spaceC,
			const Teuchos::RCP<const SpaceT>& spaceF,
			const Teuchos::Tuple<int,SpaceT::dimension>& nb ):
		op_( Teuchos::rcp( new InterT( spaceC, spaceF, nb ) ) ) {}


	void apply( const DomainFieldT& x, RangeFieldT& y ) const {


    op_->apply( x.getConst0Field(), y.get0Field() );

    int m = std::min( x.getNumberModes(), y.getNumberModes() );

    for( int i=0; i<m; ++i ) {
      op_->apply( x.getConstCField(i), y.getCField(i) );
      op_->apply( x.getConstSField(i), y.getSField(i) );
    }

    for( int i=m; i<y.getNumberModes(); ++i ) {
      y.getCField(i).initField();
      y.getSField(i).initField();
    }

	}


  void print(  std::ostream& out=std::cout ) const {

		out << "=== InterpolationMultiHarmonicOP ===\n";
		op_->print( out );

  }

}; // end of class InterpolationMultiHarmonicOp



} // end of namespace Pimpact



#ifdef COMPILE_ETI
#include "Pimpact_InterpolationOp.hpp"
#include "Pimpact_VectorFieldOpWrap.hpp"
extern template class Pimpact::InterpolationMultiHarmonicOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::InterpolationMultiHarmonicOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::InterpolationMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > > >;
extern template class Pimpact::InterpolationMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > > >;
#endif



#endif // end of #ifndef PIMPACT_INTERPOLATIONOP_HPP
