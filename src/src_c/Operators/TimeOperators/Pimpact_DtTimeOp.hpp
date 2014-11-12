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

protected:

  Scalar alpha2_;

public:

  DtTimeOp():alpha2_(1.) {};
  DtTimeOp(Scalar alpha2):alpha2_(alpha2) {};

  typedef TimeField< VectorField<SpaceT> >  DomainFieldT;
  typedef TimeField< VectorField<SpaceT> >  RangeFieldT;


  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    x.exchange();

    typename RangeFieldT::Iter j = x.sInd_;
    for( typename DomainFieldT::Iter i=y.sInd_; i<y.eInd_; ++i ) {
       (*i)->add( alpha2_, **j, -alpha2_, **( j-1 )  );
       ++j;
    }

    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DtTimeOp



/// \relates DtTimeOp
template< class SpaceT>
Teuchos::RCP< DtTimeOp<SpaceT> > createDtTimeOp( typename SpaceT::Scalar alpha2 = 1. ) {
  return( Teuchos::rcp( new DtTimeOp<SpaceT>( alpha2 ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DTTIMEOP_HPP
