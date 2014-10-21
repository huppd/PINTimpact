#pragma once
#ifndef PIMPACT_DTMODEOP_HPP
#define PIMPACT_DTMODEOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"



namespace Pimpact{


/// \ingroup ModeOperator
template<class Scalar,class Ordinal>
class DtModeOp {

  Scalar alpha2_;

public:

  DtModeOp():alpha2_(1.) {};
  DtModeOp(Scalar alpha2):alpha2_(alpha2) {};

  typedef ModeField<VectorField<Scalar,Ordinal> >  DomainFieldT;
  typedef ModeField<VectorField<Scalar,Ordinal> >  RangeFieldT;


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {
    y.getCFieldPtr()->add(      0.,  x.getConstCField(), alpha2_, x.getConstSField() );
    y.getSFieldPtr()->add( -alpha2_,  x.getConstCField(),     0., x.getConstSField() );
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DtModeOp



/// \relates DtModeOp
template< class Scalar, class Ordinal>
Teuchos::RCP< DtModeOp<Scalar,Ordinal> > createDtModeOp( Scalar alpha2 = 1. ) {
  return( Teuchos::rcp( new DtModeOp<Scalar,Ordinal>( alpha2 ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DTMODEOP_HPP
