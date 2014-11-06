#pragma once
#ifndef PIMPACT_MULTICONVECTIONJACOBIANOP_HPP
#define PIMPACT_MULTICONVECTIONJACOBIANOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_MultiConvectionOp.hpp"




namespace Pimpact {


/// \ingroup MultiHarmonicOperator
template<class SpaceT>
class MultiHarmonicConvectionJacobianOp : private MultiHarmonicConvectionOp<SpaceT> {

public:

  typedef typename MultiHarmonicConvectionOp<SpaceT>::DomainFieldT DomainFieldT;
  typedef typename MultiHarmonicConvectionOp<SpaceT>::RangeFieldT RangeFieldT;

protected:

  Teuchos::RCP<DomainFieldT> u_;

  const bool isNewton_;

public:

  MultiHarmonicConvectionJacobianOp(
      const Teuchos::RCP<const SpaceT>& space,
      const Teuchos::RCP<DomainFieldT>& u=Teuchos::null,
      const bool& isNewton=true ):
        MultiHarmonicConvectionOp<SpaceT>(space),
        u_(u),
        isNewton_(isNewton) {};

  void assignField( const DomainFieldT& mv ) {
    if( Teuchos::is_null( u_ ) )
      u_ = mv.clone();
    else
      u_->assign( mv );
  };

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    MultiHarmonicConvectionOp<SpaceT>::apply( *u_,  x,  y, true );

    if( isNewton_ )
      MultiHarmonicConvectionOp<SpaceT>::apply(  x,  *u_, y, false );

  }

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiHarmonicConvectionJacobianOp



/// \relates MultiHarmonicConvectionJacobianOp
template<class SpaceT>
Teuchos::RCP<MultiHarmonicConvectionJacobianOp<SpaceT> > createMultiHarmonicConvectionJacobianOp(
    const Teuchos::RCP<const SpaceT>& space,
    const Teuchos::RCP<typename MultiHarmonicConvectionJacobianOp<SpaceT>::DomainFieldT>& u = Teuchos::null,
    const bool& isNewton=true ) {

  return( Teuchos::rcp( new MultiHarmonicConvectionJacobianOp<SpaceT>( space, u, isNewton ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTICONVECTIONJACOBIANOP_HPP
