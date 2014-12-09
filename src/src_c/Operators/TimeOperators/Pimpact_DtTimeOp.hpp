#pragma once
#ifndef PIMPACT_DTTIMEOP_HPP
#define PIMPACT_DTTIMEOP_HPP


#include <cmath>

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_TimeField.hpp"



namespace Pimpact{


/// \ingroup TimeOperator
template<class ST>
class DtTimeOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;

  typedef typename SpaceT::Ordinal Ordinal;

  typedef TimeField< VectorField<SpaceT> >  DomainFieldT;
  typedef TimeField< VectorField<SpaceT> >  RangeFieldT;

protected:

  Teuchos::RCP<const SpaceT> space_;

  Scalar alpha2_;

public:

  DtTimeOp( const Teuchos::RCP<const SpaceT>& space ):space_(space) {
    Scalar pi = 4.*std::atan(1.);
    Scalar idt = ((Scalar)space_->nGlo()[3])/2./pi;
    alpha2_ = space_->getDomain()->getDomainSize()->getAlpha2()*idt/space_->getDomain()->getDomainSize()->getRe();
  };


  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    x.exchange();

//    typename RangeFieldT::Iter j = x.sInd_;
//    for( typename DomainFieldT::Iter i=y.sInd_; i<y.eInd_; ++i ) {

    for( Ordinal i=space_->sInd(S,3); i<space_->eInd(S,3); ++i ) {
       y.getFieldPtr(i)->add( alpha2_, x.getConstField(i), -alpha2_, x.getConstField(i-1) );
//       (*i)->add( alpha2_, **j, -alpha2_, **( j-1 )  );
//       ++j;
    }

    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DtTimeOp



///// \relates DtTimeOp
//template< class SpaceT>
//Teuchos::RCP< DtTimeOp<SpaceT> > createDtTimeOp( typename SpaceT::Scalar alpha2 = 1. ) {
//  return( Teuchos::rcp( new DtTimeOp<SpaceT>( alpha2 ) ) );
//}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DTTIMEOP_HPP
