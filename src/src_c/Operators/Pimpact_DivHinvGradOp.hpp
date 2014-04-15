#pragma once
#ifndef PIMPACT_DIVHINVGRADOP_HPP
#define PIMPACT_DIVHINVGRADOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_OperatorMV.hpp"

#include "Pimpact_LinearProblem.hpp"

namespace Pimpact{


//extern "C" {
//  void OP_div_grad( const bool& corner_yes, double* phi, double* lap );
//}

/// \brief stationary Operator for Schur complement.
/// \ingroup BaseOperator
template< class Scalar, class Ordinal >
class DivHinvGrad {
public:
  typedef ScalarField<Scalar,Ordinal>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
  typedef MultiField<VectorField<Scalar,Ordinal> > MVF;
  typedef OperatorBase<MVF> HType;
//  typedef OperatorMV<Helmholtz<Scalar,Ordinal> > HType;
//  typedef NonModeOp OpType;
private:
  Teuchos::RCP< MVF > temp0_;
  Teuchos::RCP< MVF > temp1_;
  Teuchos::RCP< Div<Scalar,Ordinal> > div_;
  Teuchos::RCP< Grad<Scalar,Ordinal> > grad_;
  Teuchos::RCP< LinearProblem<MVF> > H_;

public:
  DivHinvGrad(
      Teuchos::RCP<MVF> temp=Teuchos::null,
      Teuchos::RCP< LinearProblem<MVF > > H=Teuchos::null ):
        temp0_(temp->clone(1)), temp1_(temp->clone(1)),
        div_(Teuchos::rcp( new Div<Scalar,Ordinal> ) ),
        grad_(Teuchos::rcp( new Grad<Scalar,Ordinal> ) ),
        H_(H) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    grad_->apply(x,temp0_->getField(0) );
    H_->solve( temp1_, temp0_);
    div_->apply(temp1_->getField(0),y);
    return;
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DivHinvGrad



/// \relates DivHinvGrad
template< class Scalar, class Ordinal>
Teuchos::RCP<DivHinvGrad<Scalar,Ordinal> > createDivHinvGrad(
    Teuchos::RCP<MultiField<VectorField<Scalar,Ordinal> > > temp,
    Teuchos::RCP< LinearProblem<MultiField<VectorField<Scalar,Ordinal> > > > H ) {

  return(
      Teuchos::rcp( new DivHinvGrad<Scalar,Ordinal>( temp, H ) ) );

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVHINVGRADOP_HPP
