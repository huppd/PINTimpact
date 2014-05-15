#pragma once
#ifndef PIMPACT_MULTINONLINEARJACOBIANOP_HPP
#define PIMPACT_MULTINONLINEARJACOBIANOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_MultiNonlinearOp.hpp"




namespace Pimpact {


/// \ingroup MultiHarmonicOperator
/// \note u_ has to contain appropriate BC, temp_ and y doesnot matter, x should have zero BC
template<class S,class O>
class MultiHarmonicNonlinearJacobian : private MultiHarmonicNonlinear<S,O> {

public:

  typedef typename MultiHarmonicNonlinear<S,O>::DomainFieldT DomainFieldT;
  typedef typename MultiHarmonicNonlinear<S,O>::RangeFieldT RangeFieldT;

protected:

  Teuchos::RCP<DomainFieldT> u_;

  const bool isNewton_;

public:

  MultiHarmonicNonlinearJacobian(
//      const Teuchos::RCP<VectorField<S,O> >& temp=Teuchos::null,
      const Teuchos::RCP<DomainFieldT>& u=Teuchos::null,
      const bool& isNewton=true ):
        MultiHarmonicNonlinear<S,O>(/*temp*/),
        u_(u),
        isNewton_(isNewton) {};

  void assignField( const DomainFieldT& mv ) {
    if( Teuchos::is_null( u_ ) )
      u_ = mv.clone();
    else
      u_->assign( mv );
  };

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    MultiHarmonicNonlinear<S,O>::apply( *u_,  x,  y, true );

    if( isNewton_ )
      MultiHarmonicNonlinear<S,O>::apply(  x,  *u_, y, false );

  }

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiHarmonicNonlinearJacobian



/// \relates MultiHarmonicNonlinearJacobian
template< class S, class O>
Teuchos::RCP<MultiHarmonicNonlinearJacobian<S,O> > createMultiHarmonicNonlinearJacobian(
//    const Teuchos::RCP< VectorField<S,O> >& temp = Teuchos::null,
    const Teuchos::RCP< typename MultiHarmonicNonlinearJacobian<S,O>::DomainFieldT>& u = Teuchos::null,
    const bool& isNewton=true ) {

  return( Teuchos::rcp( new MultiHarmonicNonlinearJacobian<S,O>( /*temp,*/ u, isNewton ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTINONLINEARJACOBIANOP_HPP
