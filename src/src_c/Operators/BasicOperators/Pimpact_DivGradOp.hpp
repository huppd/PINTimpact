#pragma once
#ifndef PIMPACT_DIVGRADOP_HPP
#define PIMPACT_DIVGRADOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"




namespace Pimpact{


//extern "C" {
//  void SF_level( double* const phi );
//}


/// \brief "laplace" for pressure.
/// \ingroup BaseOperator
/// \todo not workin properly?
/// \warning does not hold test.
template<class ST>
class DivGradOp {

public:

  typedef ST SpaceT;

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  Teuchos::RCP<VectorField<SpaceT> > temp_;
  Teuchos::RCP<DivOp<SpaceT> > div_;
  Teuchos::RCP<GradOp<SpaceT> > grad_;

public:

  DivGradOp( const Teuchos::RCP<const SpaceT>& space ):
    temp_( create<VectorField>( space ) ),
    div_ ( create<DivOp>( space ) ),
    grad_( create<GradOp>( space ) ) {};

  DivGradOp(
      const Teuchos::RCP< DivOp<SpaceT> >& div,
      const Teuchos::RCP< GradOp<SpaceT> >& grad ):
    temp_( create<VectorField>( div->space() ) ),
    div_ (div),
    grad_(grad) {};

  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {

    grad_->apply( x, *temp_ );
    div_->apply( *temp_, y );

//		OP_SetBCZero(
//				space_->nLoc(),
//				space_->bl(),
//				space_->bu(),
//				space_->getDomain()->getBCLocal()->getBCL(),
//				space_->getDomain()->getBCLocal()->getBCU(),
//				space_->sIndB(i),
//				space_->eIndB(i),
//				y.getRawPtr(i) );
		SF_handle_corner(
				space()->nLoc(),
				space()->bl(),
				space()->bu(),
				space()->getDomain()->getBCLocal()->getBCL(),
				space()->getDomain()->getBCLocal()->getBCU(),
				y.getRawPtr() );
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(div_->space() ); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

}; // end of DivGradOp



/// \relates DivGradOp
template<class SpaceT>
Teuchos::RCP< DivGradOp<SpaceT> > createDivGradOp(
    const Teuchos::RCP< DivOp<SpaceT> >& div,
    const Teuchos::RCP< GradOp<SpaceT> >& grad ) {

  return(
      Teuchos::rcp( new DivGradOp<SpaceT>( div, grad ) )
  );

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVGRADOP_HPP
