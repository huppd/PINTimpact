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
template<class SpaceT>
class DivGradOp {

  Teuchos::RCP<VectorField<SpaceT> > temp_;
  Teuchos::RCP<DivOp<SpaceT> > div_;
  Teuchos::RCP<GradOp<SpaceT> > grad_;

public:

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;


  DivGradOp(
      const Teuchos::RCP<VectorField<SpaceT> >& temp,
      const Teuchos::RCP< DivOp<SpaceT> >& div,
      const Teuchos::RCP< GradOp<SpaceT> >& grad ):
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
template<class SpaceT>
Teuchos::RCP< DivGradOp<SpaceT> > createDivGradOp(
    const Teuchos::RCP<VectorField<SpaceT> >& temp,
    const Teuchos::RCP< DivOp<SpaceT> >& div,
    const Teuchos::RCP< GradOp<SpaceT> >& grad ) {
  return(
      Teuchos::rcp( new DivGradOp<SpaceT>( temp, div, grad ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVGRADOP_HPP
