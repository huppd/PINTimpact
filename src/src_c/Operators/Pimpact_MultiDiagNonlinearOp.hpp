#pragma once
#ifndef PIMPACT_MULTIDIAGNONLINEARJACOBIANOP_HPP
#define PIMPACT_MULTIDIAGNONLINEARJACOBIANOP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"

#include "Pimpact_NonlinearOp.hpp"



namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class Scalar,class Ordinal>
class MultiHarmonicDiagNonlinearJacobian {

  typedef Scalar S;
  typedef Ordinal O;

public:

  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  RangeFieldT;
  typedef NonModeOp OpType;

protected:

  Teuchos::RCP<DomainFieldT> u_;
  Teuchos::RCP<DomainFieldT> mtemp_;
  Teuchos::RCP<VectorField<S,O> > temp_;
  Teuchos::RCP<Nonlinear<S,O> > op_;

  const bool isNewton_;

public:

  MultiHarmonicDiagNonlinearJacobian( const bool& isNewton=true ):
    u_(Teuchos::null),mtemp_(Teuchos::null),temp_(Teuchos::null),
    op_( Teuchos::rcp( new Nonlinear<S,O>() )),
    isNewton_(isNewton) {};

  MultiHarmonicDiagNonlinearJacobian( const Teuchos::RCP<DomainFieldT>& temp, const bool& isNewton=true ):
    u_(temp->clone()),mtemp_(temp->clone()),temp_(temp->get0FieldPtr()->clone()),
    op_( Teuchos::rcp( new Nonlinear<S,O>() )),
    isNewton_(isNewton) {};

  void assignField( const DomainFieldT& mv ) {
    if( Teuchos::is_null( u_ ) )
       u_ = mv.clone();
     else
       u_->assign( mv );
  };

protected:

  void apply(const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z) const {
    int Nf = x.getNumberModes();
    z.init( 0. );

    // computing zero mode of y
    op_->apply( x.getConst0Field(), y.getConst0Field(), *temp_ );
    z.get0Field().add( 1., z.getConst0Field(), 1., *temp_ );

//    for( int i=1; i<=Nf; ++i ) {
//      op_->apply( x.getConstCField(i-1), y.getConstCField(i-1), *temp_ );
//      z.get0Field().add( 1., z.getConst0Field(), 0.5, *temp_ );
//    }
//
//    for( int i=1; i<=Nf; ++i ) {
//      op_->apply( x.getConstSField(i-1), y.getConstSField(i-1), *temp_ );
//      z.get0Field().add( 1., z.getConst0Field(), 0.5, *temp_ );
//    }


//    // computing cos mode of y
    for( int i=1; i<=Nf; ++i ) {
      op_->apply( x.getConst0Field(), y.getConstCField(i-1), *temp_ );
      z.getCField(i-1).add( 1., z.getConstCField(i-1), 1., *temp_ );

//      op_->apply( x.getConstCField(i-1), y.getConst0Field(), *temp_ );
//      z.getCField(i-1).add( 1., z.getConstCField(i-1), 1., *temp_ );
//
//      for( int k=1; k+i<=Nf; ++k ) {
      if( 2*i<=Nf ) {
        op_->apply( x.getConstCField(i+i-1), y.getConstCField(i-1), *temp_ );
        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );
      }
//
//        op_->apply( x.getConstCField(k-1), y.getConstCField(k+i-1), *temp_ );
//        z.getCField(i-1).add( 1., z.getCField(i-1), 0.5, *temp_ );
//
      if( 2*i<=Nf ) {
        op_->apply( x.getConstSField(i+i-1), y.getConstSField(i-1), *temp_ );
        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );
      }
//
//        op_->apply( x.getConstSField(k-1), y.getConstSField(k+i-1), *temp_ );
//        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );
//      }
    }
//
//    // computing sin mode of y
    for( int i=1; i<=Nf; ++i ) {
      op_->apply( x.getConst0Field(), y.getConstSField(i-1), *temp_ );
      z.getSField(i-1).add( 1., z.getConstSField(i-1), 1., *temp_ );

//      op_->apply( x.getConstSField(i-1), y.getConst0Field(), *temp_ );
//      z.getSField(i-1).add( 1., z.getConstSField(i-1), 1., *temp_ );
//
//      for( int k=1; k+i<=Nf; ++k ) {
      if( 2*i<=Nf ) {
        op_->apply( x.getConstCField(i+i-1), y.getConstSField(i-1), *temp_ );
        z.getSField(i-1).add( 1., z.getConstSField(i-1), -0.5, *temp_ );
      }
//
//        op_->apply( x.getConstCField(k-1), y.getConstSField(k+i-1), *temp_ );
//        z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );
//
      if( 2*i<=Nf ) {
        op_->apply( x.getConstSField(i+i-1), y.getConstCField(i-1), *temp_ );
        z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );
      }
//
//        op_->apply( x.getConstSField(k-1), y.getConstCField(k+i-1), *temp_ );
//        z.getSField(i-1).add( 1., z.getConstSField(i-1), -0.5, *temp_ );
//      }
    }

    // strange terms
//    int i;
//    for( int k=1; k<=Nf; ++k ) {
//      for( int l=1; l<=Nf; ++l ) {
//        i = k+l;
//        if( i<=Nf ) {
//          op_->apply( x.getConstCField(k-1), y.getConstCField(l-1), *temp_ );
//          z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );
//          op_->apply( x.getConstSField(k-1), y.getConstSField(l-1), *temp_ );
//          z.getCField(i-1).add( 1., z.getConstCField(i-1), -0.5, *temp_ );
//
//          op_->apply( x.getConstCField(k-1), y.getConstSField(l-1), *temp_ );
//          z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );
//          op_->apply( x.getConstSField(k-1), y.getConstCField(l-1), *temp_ );
//          z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );
//        }
//      }
//    for( int k=1; k<=Nf; ++k ) {
//        i = 2*k;
//        if( i<=Nf ) {
//          op_->apply( x.getConstCField(i-1), y.getConstCField(i-1), *temp_ );
//          z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );
//          op_->apply( x.getConstSField(i-1), y.getConstSField(i-1), *temp_ );
//          z.getCField(i-1).add( 1., z.getConstCField(i-1), -0.5, *temp_ );
//
//          op_->apply( x.getConstCField(i-1), y.getConstSField(i-1), *temp_ );
//          z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );
//          op_->apply( x.getConstSField(i-1), y.getConstCField(i-1), *temp_ );
//          z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );
//        }
//      }
//    }
  }

public:

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    if( isNewton_ ) {
      apply( *u_,  x,  *mtemp_ );
      apply(  x,  *u_, y );
      y.add( 1., *mtemp_, 1., y );
    }
    else {
      apply( *u_,  x, y );
    }
  }


  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiHarmonicDiagNonlinearJacobianOp



/// \relates MultiHarmonicDiagNonlinearJacobianOp
template< class S, class O>
Teuchos::RCP<MultiHarmonicDiagNonlinearJacobian<S,O> > createMultiHarmonicDiagNonlinearJacobian(
    const Teuchos::RCP<typename MultiHarmonicDiagNonlinearJacobian<S,O>::DomainFieldT>& u = Teuchos::null, const bool& isNewton=true ) {
  return( Teuchos::rcp( new MultiHarmonicDiagNonlinearJacobian<S,O>( u, isNewton ) ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTIDIAGNONLINEARJACOBIANOP_HPP