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
#include "Pimpact_OperatorMV.hpp"

#include "Pimpact_LinearProblem.hpp"

namespace Pimpact{



/// \brief periodic Operator for Schur complement.
/// \ingroup ModeOperator
/// \deprecated use flexicble \c DivOpGrad operator
template< class Scalar, class Ordinal >
class Div_DtLinv_Grad {
public:
  typedef ScalarField<Scalar,Ordinal>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
  typedef MultiField<ModeField<VectorField<Scalar,Ordinal> > > MVF;
  typedef OperatorMV<DtL<Scalar,Ordinal> > HType;
  typedef ModeOp OpType;
private:
  Teuchos::RCP< MVF > temp0_;
  Teuchos::RCP< MVF > temp1_;
  Teuchos::RCP< Div<Scalar,Ordinal> > div_;
  Teuchos::RCP< Grad<Scalar,Ordinal> > grad_;
  Teuchos::RCP< LinearProblem<Scalar, MVF, HType > > H_;

public:
  Div_DtLinv_Grad( Teuchos::RCP<MVF> temp,
//        Teuchos::RCP<Div<Scalar,Ordinal> > div,
//        Teuchos::RCP<Grad<Scalar,Ordinal> > grad,
      Teuchos::RCP< LinearProblem<Scalar, MVF, HType > > H ):
        temp0_(temp->clone(1)), temp1_(temp->clone(1)),
        div_(Teuchos::rcp( new Div<Scalar,Ordinal> ) ),
        grad_(Teuchos::rcp( new Grad<Scalar,Ordinal> ) ),
        H_(H) {};

  void apply(const ModeField<DomainFieldT>& x, ModeField<RangeFieldT>& y) const {
    grad_->apply( x.getConstCField(), temp0_->getFieldPtr(0)->getCField() );
    grad_->apply( x.getConstSField(), temp0_->getFieldPtr(0)->getSField() );
    H_->solve( temp1_, temp0_);
    div_->apply( temp1_->getFieldPtr(0)->getConstCField(), y.getCField() );
    div_->apply( temp1_->getFieldPtr(0)->getConstSField(), y.getSField() );
  }

  bool hasApplyTranspose() const { return( false ); }

};


template< class Scalar, class Ordinal>
Teuchos::RCP<OperatorMV< Div_DtLinv_Grad<Scalar,Ordinal> > > createDivDtLinvGrad(
    Teuchos::RCP<MultiField<ModeField<VectorField<Scalar,Ordinal> > > > temp,
    Teuchos::RCP< LinearProblem<Scalar, MultiField<ModeField<VectorField<Scalar,Ordinal> > >,
    OperatorMV<DtL<Scalar,Ordinal> > > > H ) {

  return(
      Teuchos::rcp(
          new OperatorMV<Div_DtLinv_Grad<Scalar,Ordinal> >(
              Teuchos::rcp( new Div_DtLinv_Grad<Scalar,Ordinal>( temp, H) ) ) )
  );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVDTLINVGRADOP_HPP
