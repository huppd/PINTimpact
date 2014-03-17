#pragma once
#ifndef PIMPACT_DTHELMHOLTZOP_HPP
#define PIMPACT_DTHELMHOLTZOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"



namespace Pimpact{


/// \ingroup ModeOperator
template<class Scalar,class Ordinal>
class DtL {
  Scalar omega_;
  Teuchos::RCP<Helmholtz<Scalar,Ordinal> > L_;
public:

  DtL():omega_(1.),L_(Teuchos::rcp(new Helmholtz<Scalar,Ordinal>( 0., 1. )) ) {};

  DtL( Scalar omega, Scalar mulI, Scalar mulL):omega_(omega),L_(Teuchos::rcp( new Helmholtz<Scalar,Ordinal>(mulI,mulL) ) ) {};

  DtL( Scalar omega, Teuchos::RCP<Helmholtz<Scalar,Ordinal> > L ):omega_(omega),L_(L) {};

  typedef ModeField<VectorField<Scalar,Ordinal> >  DomainFieldT;
  typedef ModeField<VectorField<Scalar,Ordinal> >  RangeFieldT;
  typedef ModeOp OpType;


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {
    L_->apply( x.getConstCField(), y.getCField() );
    y.getCFieldPtr()->add( 1., y.getConstCField(), omega_, x.getConstSField() );
    L_->apply( x.getConstSField(), y.getSField() );
    y.getSFieldPtr()->add( -omega_, x.getConstCField(), 1., y.getConstSField() );
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }
};



template< class Scalar, class Ordinal>
Teuchos::RCP< DtL<Scalar,Ordinal> > createDtL( Scalar omega=1., Scalar mulI=0., Scalar mulL=1. ) {
  return( Teuchos::rcp( new DtL<Scalar,Ordinal>( omega, mulI, mulL ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DTHELMHOLTZOP_HPP
