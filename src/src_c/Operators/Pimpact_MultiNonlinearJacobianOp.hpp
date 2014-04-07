#pragma once
#ifndef PIMPACT_MULTINONLINEARJACOBIANOP_HPP
#define PIMPACT_MULTINONLINEARJACOBIANOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_MultiNonlinearOp.hpp"



//extern "C" {
//  void OP_nonlinear( const bool& exch_yes,
//      double* phi1U, double* phi1V, double* phi1W,
//      double* phi2U, double* phi2V, double* phi2W,
//      double* nl1,   double* nl2,   double* nl3 );
//}



namespace Pimpact {


/// \ingroup BaseOperator
/// \note u_ has to contain appropriate BC, temp_ and y doesnot matter, x should have zero BC
template<class S,class O>
class MultiHarmonicNonlinearJacobian : private MultiHarmonicNonlinear<S,O> {

public:

  typedef typename MultiHarmonicNonlinear<S,O>::DomainFieldT DomainFieldT;
  typedef typename MultiHarmonicNonlinear<S,O>::RangeFieldT RangeFieldT;

protected:

  Teuchos::RCP<DomainFieldT> u_;
  Teuchos::RCP<DomainFieldT> temp_;

  const bool isNewton_;

public:

  MultiHarmonicNonlinearJacobian(
      const Teuchos::RCP<VectorField<S,O> >& temp=Teuchos::null,
      const Teuchos::RCP<DomainFieldT>& u=Teuchos::null,
      const bool& isNewton=true ):
        MultiHarmonicNonlinear<S,O>(temp),u_(u),temp_(u->clone()),isNewton_(isNewton) {};

  void assignField( const DomainFieldT& mv ) {
    if( Teuchos::is_null( u_ ) )
      u_ = mv.clone();
    else
      u_->assign( mv );
  };

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    if( isNewton_ ) {
      MultiHarmonicNonlinear<S,O>::apply( *u_,  x,  *temp_ );
      MultiHarmonicNonlinear<S,O>::apply(  x,  *u_, y );
      y.add( 1., *temp_, 1., y );
    }
    else {
      MultiHarmonicNonlinear<S,O>::apply( *u_,  x, y );
    }
  }

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiHarmonicNonlinearJacobian



/// \relates MultiHarmonicNonlinearJacobian
template< class S, class O>
Teuchos::RCP<MultiHarmonicNonlinearJacobian<S,O> > createMultiHarmonicNonlinearJacobian(
    const Teuchos::RCP< VectorField<S,O> >& temp = Teuchos::null,
    const Teuchos::RCP< typename MultiHarmonicNonlinearJacobian<S,O>::DomainFieldT>& u = Teuchos::null,
    const bool& isNewton=true ) {

  return( Teuchos::rcp( new MultiHarmonicNonlinearJacobian<S,O>( temp, u, isNewton ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTINONLINEARJACOBIANOP_HPP
