#pragma once
#ifndef PIMPACT_DTOP_HPP
#define PIMPACT_DTOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"



namespace Pimpact{


/// \ingroup ModeOperator
template<class Scalar,class Ordinal>
class Dt {
  Scalar omega_;
public:
  Dt():omega_(1.) {};
  Dt(Scalar omega):omega_(omega) {};

  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;
  typedef ModeOp OpType;


  void apply(const ModeField<DomainFieldT>& x, ModeField<RangeFieldT>& y ) const {
    y.getCFieldPtr()->add(      0.,  x.getConstCField(), omega_, x.getConstSField() );
    y.getSFieldPtr()->add( -omega_,  x.getConstCField(),     0., x.getConstSField() );
//    y.getFieldC()->assign( x.getConstFieldS() );
//    y.getFieldC()->scale( -omega_ );
//    y.getFieldS()->assign( x.getConstFieldC() );
//    y.getFieldS()->scale( omega_ );
  }
  bool hasApplyTranspose() const { return( false ); }
};


template< class Scalar, class Ordinal>
Teuchos::RCP<OperatorMV< Dt<Scalar,Ordinal> > > createDt( Scalar omega = 1. ) {
  return( Teuchos::rcp( new OperatorMV<Dt<Scalar,Ordinal> >( Teuchos::rcp( new Dt<Scalar,Ordinal>( omega ) ) ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DTDOP_HPP
