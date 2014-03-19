#pragma once
#ifndef PIMPACT_MULTIDTHELMHOLTZOP_HPP
#define PIMPACT_MULTIDTHELMHOLTZOP_HPP


#include "Pimpact_Types.hpp"

//#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

#include "Pimpact_DtHelmholtzOp.hpp"



namespace Pimpact{


/// \ingroup MultiHarmonicOperator
/// for generalizing this(MultiHarmonicwrapper for ModeOperator), ModeOperator needs a non Mode representation
template<class Scalar,class Ordinal>
class MultiDtHelmholtz {

//  Scalar omega_;
  Teuchos::RCP<DtL<Scalar,Ordinal> > op_;

public:

  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  RangeFieldT;


  MultiDtHelmholtz( Scalar omega=1., Scalar mulI=0., Scalar mulL=1.):
    op_( createDtL<Scalar,Ordinal>( omega, mulI, mulL )) {};

  MultiDtHelmholtz( const Teuchos::RCP<DtL<Scalar,Ordinal> >& op ):
    op_(op) {};


  void apply(const DomainFieldT& x, RangeFieldT& y, int k=1 ) const {
    op_->getInnerOpPtr()->apply( x.getConst0Field(), y.get0Field() );

    for( int i=0; i<x.getNumberModes(); ++i ) {
      op_->apply( x.getConstField(i), y.getField(i) );
      op_->apply( x.getConstField(i), y.getField(i) );
    }
  }


  void assignField( const DomainFieldT& mv ) {};


  Teuchos::RCP< Helmholtz<Scalar,Ordinal> > getInnerOpPtr() {
    return( op_ );
  }


  bool hasApplyTranspose() const { return( false ); }
};



template< class Scalar, class Ordinal>
Teuchos::RCP< MultiDtHelmholtz<Scalar,Ordinal> > createMultiDtHelmholtz( Scalar omega=1., Scalar mulI=0., Scalar mulL=1. ) {
  return( Teuchos::rcp( new MultiDtHelmholtz<Scalar,Ordinal>( omega, mulI, mulL ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIDTHELMHOLTZOP_HPP
