#pragma once
#ifndef PIMPACT_MULTIDTHELMHOLTZOP_HPP
#define PIMPACT_MULTIDTHELMHOLTZOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"

#include "Pimpact_DtLapOp.hpp"


namespace Pimpact{


/// \ingroup MultiHarmonicOperator
/// for generalizing this(MultiHarmonicwrapper for ModeOperator), ModeOperator needs a non Mode representation
template<class Scalar,class Ordinal>
class MultiDtHelmholtz {

  Teuchos::RCP<DtLapOp<Scalar,Ordinal> > op_;

public:

  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  RangeFieldT;


  MultiDtHelmholtz( const Teuchos::RCP< const Space<Scalar,Ordinal,3> >& space, Scalar alpha2=1., Scalar iRe=1.):
    op_( createDtLapOp<Scalar,Ordinal>( space, alpha2, iRe) ) {};

  MultiDtHelmholtz( const Teuchos::RCP<DtLapOp<Scalar,Ordinal> >& op ):
    op_(op) {};


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {
    op_->getInnerOpPtr()->apply( x.getConst0Field(), y.get0Field() );

    for( int i=0; i<x.getNumberModes(); ++i ) {
      op_->apply( x.getConstField(i), y.getField(i), i+1 );
      op_->apply( x.getConstField(i), y.getField(i), i+1 );
    }
  }


  void assignField( const DomainFieldT& mv ) {};


  Teuchos::RCP< HelmholtzOp<Scalar,Ordinal> > getInnerOpPtr() {
    return( op_ );
  }


  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiDtHelmholtz



/// \relates MultiDtHelmholtz
template< class S, class O>
Teuchos::RCP< MultiDtHelmholtz<S,O> > createMultiDtHelmholtz(
    const Teuchos::RCP< const Space<S,O,3> >& space,
    S omega=1.,
    S mulL=1. ) {
  return( Teuchos::rcp( new MultiDtHelmholtz<S,O>( space, omega, mulL ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIDTHELMHOLTZOP_HPP
