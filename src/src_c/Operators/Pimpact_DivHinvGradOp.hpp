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

namespace Pimpact{


//extern "C" {
//  void OP_div_grad( const bool& corner_yes, double* phi, double* lap );
//}

/// \brief stationary Operator for Schur complement.
/// \ingroup BaseOperator
/// \deprecated use flexicble \c DivOpGrad operator
template< class Scalar, class Ordinal >
class Div_Hinv_Grad {
public:
  typedef ScalarField<Scalar,Ordinal>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
  typedef MultiField<VectorField<Scalar,Ordinal> > MVF;
  typedef OperatorMV<Helmholtz<Scalar,Ordinal> > HType;
  typedef NonModeOp OpType;
private:
  Teuchos::RCP< MVF > temp0_;
  Teuchos::RCP< MVF > temp1_;
  Teuchos::RCP< Div<Scalar,Ordinal> > div_;
  Teuchos::RCP< Grad<Scalar,Ordinal> > grad_;
  Teuchos::RCP< LinearProblem<Scalar, MVF, HType > > H_;

public:
  Div_Hinv_Grad( Teuchos::RCP<MVF> temp,
//        Teuchos::RCP<Div<Scalar,Ordinal> > div,
//        Teuchos::RCP<Grad<Scalar,Ordinal> > grad,
      Teuchos::RCP< LinearProblem<Scalar, MVF, HType > > H ):
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
  bool hasApplyTranspose() const { return( false ); }

};


template< class Scalar, class Ordinal>
Teuchos::RCP<OperatorMV< Div_Hinv_Grad<Scalar,Ordinal> > > createDivHinvGrad(
    Teuchos::RCP<MultiField<VectorField<Scalar,Ordinal> > > temp,
    Teuchos::RCP< LinearProblem<Scalar, MultiField<VectorField<Scalar,Ordinal> >, OperatorMV<Helmholtz<Scalar,Ordinal> > > > H ) {

  return(
      Teuchos::rcp( new OperatorMV<Div_Hinv_Grad<Scalar,Ordinal> >( Teuchos::rcp( new Div_Hinv_Grad<Scalar,Ordinal>( temp, H ) ) ) )
  );

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVHINVGRADOP_HPP
