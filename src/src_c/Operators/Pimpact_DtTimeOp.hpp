#pragma once
#ifndef PIMPACT_DTTIMEOP_HPP
#define PIMPACT_DTTIMEOP_HPP


#include <cmath>

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_TimeField.hpp"



namespace Pimpact{


/// \ingroup TimeOperator
template<class Scalar,class Ordinal>
class DtTimeOp {

  Scalar alpha2_;

public:

  DtTimeOp():alpha2_(1.) {};
  DtTimeOp(Scalar alpha2):alpha2_(alpha2) {};

  typedef TimeField< VectorField<Scalar,Ordinal,4> >  DomainFieldT;
  typedef TimeField< VectorField<Scalar,Ordinal,4> >  RangeFieldT;


  void apply( const DomainFieldT& x, RangeFieldT& y ) const {
    const_cast<DomainFieldT&>(x).exchange();

    Scalar pi = 4.*std::atan(1.);

    Scalar idt = ((Scalar)y.getSpace()->nGlo()[3])/2./pi;

    typename RangeFieldT::Iter j = x.beginI_;
    for( typename DomainFieldT::Iter i=y.beginI_; i<y.endI_; ++i ) {
       (*i)->add( idt*alpha2_, **j, -idt*alpha2_, **( j-1 )  );
       ++j;
    }
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DtTimeOp



/// \relates DtTimeOp
template< class Scalar=double, class Ordinal=int>
Teuchos::RCP< DtTimeOp<Scalar,Ordinal> > createDtTimeOp( Scalar alpha2 = 1. ) {
  return( Teuchos::rcp( new DtTimeOp<Scalar,Ordinal>( alpha2 ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DTTIMEOP_HPP
