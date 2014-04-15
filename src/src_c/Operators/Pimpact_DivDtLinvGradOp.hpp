#pragma once
#ifndef PIMPACT_DIVDTLINVGRADOP_HPP
#define PIMPACT_DIVDTLINVGRADOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_ModeField.hpp"

#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_DtHelmholtzOp.hpp"

//#include "Pimpact_OperatorMV.hpp"
//#include "Pimpact_ModeOpWrap.hpp"
//#include "Pimpact_MultiOpWrap.hpp"

#include "Pimpact_LinearProblem.hpp"




namespace Pimpact{


/// \brief periodic Operator for Schur complement.
/// \ingroup ModeOperator
/// \depracted
template< class Scalar, class Ordinal >
class DivDtLinvGrad {

public:

  typedef ModeField< ScalarField<Scalar,Ordinal> > DomainFieldT;
  typedef ModeField< ScalarField<Scalar,Ordinal> > RangeFieldT;
  typedef MultiField<ModeField<VectorField<Scalar,Ordinal> > > MVF;
//  typedef ModeOp OpType;

private:
  Teuchos::RCP< MVF > temp0_;
  Teuchos::RCP< MVF > temp1_;
  Teuchos::RCP< Div<Scalar,Ordinal> >  div_;
  Teuchos::RCP< Grad<Scalar,Ordinal> > grad_;
  Teuchos::RCP< LinearProblem<MVF> > H_;


public:
  DivDtLinvGrad():
        temp0_(Teuchos::null), temp1_(Teuchos::null),
        div_(Teuchos::rcp( new Div<Scalar,Ordinal> ) ),
        grad_(Teuchos::rcp( new Grad<Scalar,Ordinal> ) ),
        H_(Teuchos::null) {};

  DivDtLinvGrad( const Teuchos::RCP<MVF>& temp,
      const Teuchos::RCP< LinearProblem<MVF> >& H ):
        temp0_(temp->clone(1)), temp1_(temp->clone(1)),
        div_(Teuchos::rcp( new Div<Scalar,Ordinal> ) ),
        grad_(Teuchos::rcp( new Grad<Scalar,Ordinal> ) ),
        H_(H) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    grad_->apply( x.getConstCField(), temp0_->getFieldPtr(0)->getCField() );
    grad_->apply( x.getConstSField(), temp0_->getFieldPtr(0)->getSField() );
    H_->solve( temp1_, temp0_);
    div_->apply( temp1_->getFieldPtr(0)->getConstCField(), y.getCField() );
    div_->apply( temp1_->getFieldPtr(0)->getConstSField(), y.getSField() );
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

};


template< class Scalar, class Ordinal>
Teuchos::RCP< DivDtLinvGrad<Scalar,Ordinal> > createDivDtLinvGrad(
    const Teuchos::RCP<MultiField<ModeField<VectorField<Scalar,Ordinal> > > >& temp=Teuchos::null,
    const Teuchos::RCP< LinearProblem<MultiField<ModeField<VectorField<Scalar,Ordinal> > > > >& H=Teuchos::null ) {
  return(
      Teuchos::rcp( new DivDtLinvGrad<Scalar,Ordinal>( temp, H ) )
  );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVDTLINVGRADOP_HPP
