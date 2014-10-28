#pragma once
#ifndef PIMPACT_MULTICONVECTIONJACOBIANOP_HPP
#define PIMPACT_MULTICONVECTIONJACOBIANOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_MultiConvectionOp.hpp"




namespace Pimpact {


/// \ingroup MultiHarmonicOperator
template<class S,class O>
class MultiHarmonicConvectionJacobianOp : private MultiHarmonicConvectionOp<S,O> {

public:

  typedef typename MultiHarmonicConvectionOp<S,O>::DomainFieldT DomainFieldT;
  typedef typename MultiHarmonicConvectionOp<S,O>::RangeFieldT RangeFieldT;

protected:

  Teuchos::RCP<DomainFieldT> u_;

  const bool isNewton_;

public:

  MultiHarmonicConvectionJacobianOp(
      const Teuchos::RCP<const Space<S,O,3> >& space,
      const Teuchos::RCP<DomainFieldT>& u=Teuchos::null,
      const bool& isNewton=true ):
        MultiHarmonicConvectionOp<S,O>(space),
        u_(u),
        isNewton_(isNewton) {};

  void assignField( const DomainFieldT& mv ) {
    if( Teuchos::is_null( u_ ) )
      u_ = mv.clone();
    else
      u_->assign( mv );
  };

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    MultiHarmonicConvectionOp<S,O>::apply( *u_,  x,  y, true );

    if( isNewton_ )
      MultiHarmonicConvectionOp<S,O>::apply(  x,  *u_, y, false );

  }

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiHarmonicConvectionJacobianOp



/// \relates MultiHarmonicConvectionJacobianOp
template< class S=double , class O=int >
Teuchos::RCP<MultiHarmonicConvectionJacobianOp<S,O> > createMultiHarmonicConvectionJacobianOp(
    const Teuchos::RCP<const Space<S,O,3> >& space,
    const Teuchos::RCP<typename MultiHarmonicConvectionJacobianOp<S,O>::DomainFieldT>& u = Teuchos::null,
    const bool& isNewton=true ) {

  return( Teuchos::rcp( new MultiHarmonicConvectionJacobianOp<S,O>( space, u, isNewton ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTICONVECTIONJACOBIANOP_HPP
