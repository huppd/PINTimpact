#pragma once
#ifndef PIMPACT_DTMODEOP_HPP
#define PIMPACT_DTMODEOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"



namespace Pimpact{


/// \ingroup ModeOperator
template<class ST>
class DtModeOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

protected:

  Scalar alpha2_;

public:

  DtModeOp():alpha2_(1.) {};
  DtModeOp(Scalar alpha2):alpha2_(alpha2) {};

  typedef ModeField<VectorField<SpaceT> >  DomainFieldT;
  typedef ModeField<VectorField<SpaceT> >  RangeFieldT;


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {
    y.getCFieldPtr()->add(      0.,  x.getConstCField(), alpha2_, x.getConstSField() );
    y.getSFieldPtr()->add( -alpha2_,  x.getConstCField(),     0., x.getConstSField() );
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DtModeOp



/// \relates DtModeOp
template<class SpaceT>
Teuchos::RCP< DtModeOp<SpaceT> > createDtModeOp( typename SpaceT::Scalar alpha2 = 1. ) {
  return( Teuchos::rcp( new DtModeOp<SpaceT>( alpha2 ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DTMODEOP_HPP
