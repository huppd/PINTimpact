#pragma once
#ifndef PIMPACT_MULTIDIAGNONLINEARJACOBIANOP_HPP
#define PIMPACT_MULTIDIAGNONLINEARJACOBIANOP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"

#include "Pimpact_ConvectionVOp.hpp"



namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class Scalar,class Ordinal>
class MultiHarmonicDiagNonlinearJacobian {

  typedef Scalar S;
  typedef Ordinal O;

public:

  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  RangeFieldT;

protected:

  Teuchos::RCP<DomainFieldT> u_;
  //  Teuchos::RCP<DomainFieldT> mtemp_;
  //  Teuchos::RCP<VectorField<S,O> > temp_;
  Teuchos::RCP<ConvectionOp<S,O> > op_;

  const bool isNewton_;

public:

  MultiHarmonicDiagNonlinearJacobian( const bool& isNewton=true ):
    u_(Teuchos::null),
    //    mtemp_(Teuchos::null),
    //    temp_(Teuchos::null),
    op_( Teuchos::rcp( new ConvectionOp<S,O>() )),
    isNewton_(isNewton) {};

  MultiHarmonicDiagNonlinearJacobian(
      const Teuchos::RCP<const Space<S,O,3> >& space,
      const Teuchos::RCP<DomainFieldT>& temp,
      const bool& isNewton=true ):
        u_(temp->clone()),
        //    mtemp_(temp->clone()),
        //    temp_(temp->get0FieldPtr()->clone()),
        op_( Teuchos::rcp( new ConvectionOp<S,O>( space ) )),
        isNewton_(isNewton) {};

  void assignField( const DomainFieldT& mv ) {
    if( Teuchos::is_null( u_ ) )
      u_ = mv.clone();
    else
      u_->assign( mv );
  };

protected:

  void apply(const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {
    int Nf = x.getNumberModes();
    if( init_yes )
      z.init( 0. );

    // computing zero mode of y
    op_->apply( x.getConst0Field(), y.getConst0Field(), z.get0Field(), 1. );
    //    z.get0Field().add( 1., z.getConst0Field(), 1., *temp_ );


    for( int i=1; i<=Nf; ++i ) {
      // computing cos mode of y
      op_->apply( x.getConst0Field(), y.getConstCField(i-1), z.getCField(i-1), 1. );
      //      z.getCField(i-1).add( 1., z.getConstCField(i-1), 1., *temp_ );

      // computing sin mode of y
      op_->apply( x.getConst0Field(), y.getConstSField(i-1), z.getSField(i-1), 1. );
      //      z.getSField(i-1).add( 1., z.getConstSField(i-1), 1., *temp_ );

      if( 2*i<=Nf ) {

        // computing cos mode of y
        op_->apply( x.getConstCField(i+i-1), y.getConstCField(i-1), z.getCField(i-1), 0.5 );
        //        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );

        op_->apply( x.getConstSField(i+i-1), y.getConstSField(i-1), z.getCField(i-1), 0.5 );
        //        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );

        // computing sin mode of y
        op_->apply( x.getConstCField(i+i-1), y.getConstSField(i-1), z.getSField(i-1), -0.5 );
        //        z.getSField(i-1).add( 1., z.getConstSField(i-1), -0.5, *temp_ );

        op_->apply( x.getConstSField(i+i-1), y.getConstCField(i-1), z.getSField(i-1), 0.5 );
        //        z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );

      }

    }
  }

public:

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    apply( *u_,  x,  y, true );

    if( isNewton_ )
      apply(  x,  *u_, y, false );

  }


  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiHarmonicDiagNonlinearJacobianOp



/// \relates MultiHarmonicDiagNonlinearJacobianOp
template< class S, class O>
Teuchos::RCP<MultiHarmonicDiagNonlinearJacobian<S,O> > createMultiHarmonicDiagNonlinearJacobian(
    const Teuchos::RCP<const Space<S,O,3> >& space,
    const Teuchos::RCP<typename MultiHarmonicDiagNonlinearJacobian<S,O>::DomainFieldT>& u = Teuchos::null,
    const bool& isNewton=true ) {
  return( Teuchos::rcp( new MultiHarmonicDiagNonlinearJacobian<S,O>( space, u, isNewton ) ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTIDIAGNONLINEARJACOBIANOP_HPP
