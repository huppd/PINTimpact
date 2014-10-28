#pragma once
#ifndef PIMPACT_TIMENONLINEARJACOBIANOP_HPP
#define PIMPACT_TIMENONLINEARJACOBIANOP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_TimeField.hpp"

#include "Pimpact_ConvectionVOp.hpp"



namespace Pimpact {


/// \ingroup TimeOperator
/// \todo better design would be to wrap ConvectionJacobianOp(no because then you would store stencil mutliple times, so one would have to abstract Jacobian)
/// \note u_ has to contain appropriate BC, temp_ and y doesnot matter, x should have zero BC
template<class Scalar,class Ordinal,bool CrankNicolsonYes=false >
class TimeNonlinearJacobian {

public:

  typedef TimeField< VectorField<Scalar,Ordinal,4> > DomainFieldT;
  typedef TimeField< VectorField<Scalar,Ordinal,4> > RangeFieldT;

protected:

  Teuchos::RCP<DomainFieldT> u_;

  const bool isNewton_;

  Teuchos::RCP< VectorField<Scalar,Ordinal,4> > temp_;

  Teuchos::RCP<const ConvectionVOp<Scalar,Ordinal,4> > op_;

public:

  TimeNonlinearJacobian(
      const Teuchos::RCP<const Space<Scalar,Ordinal,4> >& space,
      const bool& isNewton=true ):
    u_( Teuchos::null ),
    isNewton_(isNewton),
    temp_( createVectorField<Scalar,Ordinal,4>(space) ),
    op_( createConvectionVOp<Scalar,Ordinal,4>(space) )
  {};

  void assignField( const DomainFieldT& mv ) {

    if( u_.is_null() )
      u_ = mv.clone();
    else
      u_->assign( mv );

    u_->exchange();

    for( typename DomainFieldT::Iter i = u_->mfs_.begin(); i<u_->mfs_.end(); ++i )
      (*i)->exchange();

  };

  void apply( const DomainFieldT& x, RangeFieldT& y, Scalar mul=0. ) const {

    x.exchange();

    if( CrankNicolsonYes ) {

      typename DomainFieldT::Iter k = u_->mfs_.begin();
      typename RangeFieldT::Iter  j = y.  mfs_.begin();
//      typename RangeFieldT::Iter  j = y.sInd_-1;
      for( typename DomainFieldT::Iter i = const_cast<DomainFieldT&>(x).mfs_.begin(); i<x.eInd_; ++i ) {

        (*i)->exchange();
        temp_->init(0);

        op_->apply( **k, **i, *temp_, 1. );
//        OP_nonlinear(
//            (*k)->vec(0), (*k)->vec(1), (*k)->vec(2),
//            (*i)->vec(0), (*i)->vec(1), (*i)->vec(2),
//            temp_->vec(0), temp_->vec(1), temp_->vec(2),
//            1. );

        if( isNewton_ ) {
          op_->apply( **i, **k, *temp_, 1. );
//          OP_nonlinear(
//              (*i)->vec(0), (*i)->vec(1), (*i)->vec(2),
//              (*k)->vec(0), (*k)->vec(1), (*k)->vec(2),
//              temp_->vec(0), temp_->vec(1), temp_->vec(2),
//              1. );
        }
        if( j>=y.sInd_ )
          (*j)->add( 1., **j, 0.5, *temp_ );
        ++j;
        if( j<y.eInd_ ) {
          (*j)->add( 0., **j, 0.5, *temp_ );
          (*j)->changed();
        }
//        ++j;
        ++k;
      }
    }
    else {
      y.init(0.);

      typename DomainFieldT::Iter k = u_->sInd_;
      typename RangeFieldT::Iter  j = y.sInd_;
      for( typename DomainFieldT::Iter i = const_cast<DomainFieldT&>(x).sInd_; i<x.eInd_; ++i ) {
        (*i)->exchange();
        op_->apply( **k, **i, **j, 1. );
//        OP_nonlinear(
//            (*k)->vec(0), (*k)->vec(1), (*k)->vec(2),
//            (*i)->vec(0), (*i)->vec(1), (*i)->vec(2),
//            (*j)->vec(0), (*j)->vec(1), (*j)->vec(2),
//            1. );

        if( isNewton_ ) {
          op_->apply( **i, **k, **j, 1. );
//          OP_nonlinear(
//              (*i)->vec(0), (*i)->vec(1), (*i)->vec(2),
//              (*k)->vec(0), (*k)->vec(1), (*k)->vec(2),
//              (*j)->vec(0), (*j)->vec(1), (*j)->vec(2),
//              1. );
        }
        (*j)->changed();
        ++j; ++k;
      }
    }

    y.changed();

  }

  bool hasApplyTranspose() const { return( false ); }


}; // end of class TimeNonlinearJacobian



/// \relates TimeNonlinearJacobian
template< class S=double, class O=int, bool CNY=false >
Teuchos::RCP<TimeNonlinearJacobian<S,O,CNY> > createTimeNonlinearJacobian(
    const Teuchos::RCP<const Space<S,O,4> >& space,
    const bool& isNewton=true ) {

    return( Teuchos::rcp( new TimeNonlinearJacobian<S,O,CNY>( space, isNewton ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_TIMENONLINEARJACOBIANOP_HPP
