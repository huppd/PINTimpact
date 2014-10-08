#pragma once
#ifndef PIMPACT_COMPOUNDSTOKESOP_HPP
#define PIMPACT_COMPOUNDSTOKESOP_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_CompoundField.hpp"

#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_OperatorMV.hpp"




namespace Pimpact {


/// \ingroup CompoundOperator
template<class Scalar,class Ordinal>
class CompoundStokes {
  Scalar omega_;
  Teuchos::RCP<HelmholtzOp<Scalar,Ordinal> > L_;
  Teuchos::RCP<DivOp<Scalar,Ordinal> >  div_;
  Teuchos::RCP<GradOp<Scalar,Ordinal> > grad_;


  typedef ScalarField<Scalar,Ordinal>  SF;
  typedef VectorField<Scalar,Ordinal>  VF;
  typedef ModeField<SF>                MSF;
  typedef ModeField<VF>                MVF;

  Teuchos::RCP<VF> temp_;

public:

  typedef CompoundField<MVF,MSF>  DomainFieldT;
  typedef CompoundField<MVF,MSF>  RangeFieldT;

//  CompoundStokes():
//    omega_(1.),
//    L_(Teuchos::rcp(new HelmholtzOp<Scalar,Ordinal>( 0., 1. )) ),
//    div_(Teuchos::rcp(new DivOp<Scalar,Ordinal>())),
//    grad_(Teuchos::rcp(new GradOp<Scalar,Ordinal>() )),
//    temp_(Teuchos::null) {};

  CompoundStokes(
      const Teuchos::RCP<const Space<Scalar,Ordinal> >& space,
      Scalar omega,
      Scalar mulI,
      Scalar mulL,
      Teuchos::RCP<VF> temp ):
    omega_(omega),
    L_(Teuchos::rcp( new HelmholtzOp<Scalar,Ordinal>(space, mulI,mulL) ) ),
    div_(Teuchos::rcp(new DivOp<Scalar,Ordinal>(space))),
    grad_(Teuchos::rcp(new GradOp<Scalar,Ordinal>(space) )),
    temp_(temp) {};


  void apply(const DomainFieldT& x_, RangeFieldT& y_,
      Belos::ETrans trans=Belos::NOTRANS  ) const {
    auto x = x_.getConstVFieldPtr();
    auto y = y_.getVFieldPtr();
    auto xp = x_.getConstSFieldPtr();
    auto yp = y_.getSFieldPtr();
    // H-blockz
    L_->apply( x->getConstCField(), y->getCField() );
    y->getCFieldPtr()->add( 1., y->getConstCField(), omega_, x->getConstSField() );
    L_->apply( x->getConstSField(), y->getSField() );
    y->getSFieldPtr()->add( -omega_, x->getConstCField(), 1., y->getConstSField() );
    // div
    div_->apply( x->getConstCField(), yp->getCField() );
    div_->apply( x->getConstSField(), yp->getSField() );
    // grad pressure
    grad_->apply( xp->getConstCField(), *temp_ );
    y->getCFieldPtr()->add( 1., y->getCField(), 1., *temp_ );
    grad_->apply( xp->getConstSField(), *temp_ );
    y->getSFieldPtr()->add( 1., y->getSField(), 1., *temp_ );
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class CompoundStokes



/// \relates CompoundStokes
template< class Scalar, class Ordinal>
Teuchos::RCP<CompoundStokes<Scalar,Ordinal> > createCompoundStokes(
    const Teuchos::RCP<const Space<Scalar,Ordinal> >& space,
    Scalar omega, Scalar mulI, Scalar mulL, Teuchos::RCP<VectorField<Scalar,Ordinal> > temp ) {

  return(
      Teuchos::rcp( new CompoundStokes<Scalar,Ordinal>( space, omega,mulI,mulL,temp) )
      );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDSTOKESOP_HPP
