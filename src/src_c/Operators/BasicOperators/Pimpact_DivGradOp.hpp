#pragma once
#ifndef PIMPACT_DIVGRADOP_HPP
#define PIMPACT_DIVGRADOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"




namespace Pimpact{


extern "C" {

  void SF_level( double* const phi );

}


/// \brief "laplace" for pressure.
/// \ingroup BaseOperator
/// \todo not workin properly?
/// \todo add temporary variable
/// \warning does not hold test.
template<class Scalar,class Ordinal, int dimension=3>
class DivGradOp {

  Teuchos::RCP<VectorField<Scalar,Ordinal,dimension> > temp_;
  Teuchos::RCP<DivOp<Scalar,Ordinal,dimension> > div_;
  Teuchos::RCP<GradOp<Scalar,Ordinal,dimension> > grad_;

public:

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  RangeFieldT;


  DivGradOp(
      const Teuchos::RCP<VectorField<Scalar,Ordinal,dimension> >& temp,
      const Teuchos::RCP< DivOp<Scalar,Ordinal,dimension> >& div,
      const Teuchos::RCP< GradOp<Scalar,Ordinal,dimension> >& grad ):
    temp_(temp),
    div_ (div),
    grad_(grad) {};

  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {
    grad_->apply( x, *temp_ );
    div_->apply( *temp_, y );

//    OP_div_grad( true, x.s_, y.s_ );
//    SF_level( y.s_ );
//    y.change();

  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of DivGradOp



/// \relates DivGradOp
template<class S, class O, int d=3>
Teuchos::RCP< DivGradOp<S,O,d> > createDivGradOp(
    const Teuchos::RCP<VectorField<S,O,d> >& temp,
    const Teuchos::RCP< DivOp<S,O,d> >& div,
    const Teuchos::RCP< GradOp<S,O,d> >& grad ) {
  return(
      Teuchos::rcp( new DivGradOp<S,O,d>( temp, div, grad ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVGRADOP_HPP
