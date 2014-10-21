#pragma once
#ifndef PIMPACT_TIMENONLINEARJACOBIANOP_HPP
#define PIMPACT_TIMENONLINEARJACOBIANOP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_TimeField.hpp"

//#include "Pimpact_FieldFactory.hpp"




extern "C" {

void OP_nonlinear(
      double* const phi1U, double* const phi1V, double* const phi1W,
      double* const phi2U, double* const phi2V, double* const phi2W,
      double* const nl1,   double* const nl2,   double* const nl3,
      const double& mul );

}



namespace Pimpact {


/// \ingroup TimeOperator
/// \todo better design would be to wrap NonlinearJacobian
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

public:

  TimeNonlinearJacobian( const bool& isNewton=true ):
    u_(Teuchos::null),
    isNewton_(isNewton) {};

  TimeNonlinearJacobian(
      const Teuchos::RCP<DomainFieldT>& u,
      const bool& isNewton=true,
      const Teuchos::RCP< VectorField<Scalar,Ordinal,4> >& temp=Teuchos::null ):
    u_( u ),
    isNewton_(isNewton),
    temp_(temp) {};

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

    const_cast<DomainFieldT&>(x).exchange();


    if( CrankNicolsonYes ) {

      typename DomainFieldT::Iter k = u_->mfs_.begin();
      typename RangeFieldT::Iter  j = y.  mfs_.begin();
//      typename RangeFieldT::Iter  j = y.beginI_-1;
      for( typename DomainFieldT::Iter i = const_cast<DomainFieldT&>(x).mfs_.begin(); i<x.endI_; ++i ) {
        (*i)->exchange();
        temp_->init(0);

        OP_nonlinear(
            (*k)->vec(0), (*k)->vec(1), (*k)->vec(2),
            (*i)->vec(0), (*i)->vec(1), (*i)->vec(2),
            temp_->vec(0), temp_->vec(1), temp_->vec(2),
            1. );

        if( isNewton_ ) {
          OP_nonlinear(
              (*i)->vec(0), (*i)->vec(1), (*i)->vec(2),
              (*k)->vec(0), (*k)->vec(1), (*k)->vec(2),
              temp_->vec(0), temp_->vec(1), temp_->vec(2),
              1. );
        }
        if( j>=y.beginI_ )
          (*j)->add( 1., **j, 0.5, *temp_ );
        ++j;
        if( j<y.endI_ ) {
          (*j)->add( 0., **j, 0.5, *temp_ );
          (*j)->changed();
        }
//        ++j;
        ++k;
      }
    }
    else {
      y.init(0.);

      typename DomainFieldT::Iter k = u_->beginI_;
      typename RangeFieldT::Iter  j = y.beginI_;
      for( typename DomainFieldT::Iter i = const_cast<DomainFieldT&>(x).beginI_; i<x.endI_; ++i ) {
        (*i)->exchange();
        OP_nonlinear(
            (*k)->vec(0), (*k)->vec(1), (*k)->vec(2),
            (*i)->vec(0), (*i)->vec(1), (*i)->vec(2),
            (*j)->vec(0), (*j)->vec(1), (*j)->vec(2),
            1. );

        if( isNewton_ ) {
          OP_nonlinear(
              (*i)->vec(0), (*i)->vec(1), (*i)->vec(2),
              (*k)->vec(0), (*k)->vec(1), (*k)->vec(2),
              (*j)->vec(0), (*j)->vec(1), (*j)->vec(2),
              1. );
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
    const Teuchos::RCP<typename TimeNonlinearJacobian<S,O>::DomainFieldT>& u = Teuchos::null,
    const bool& isNewton=true,
    const Teuchos::RCP< VectorField<S,O,4> >& temp=Teuchos::null ) {

//  if( Teuchos::is_null(u) )
//    return( Teuchos::rcp( new TimeNonlinearJacobian<S,O,CNY>( isNewton, temp ) ) );
//  else
    return( Teuchos::rcp( new TimeNonlinearJacobian<S,O,CNY>( u, isNewton, temp ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_TIMENONLINEARJACOBIANOP_HPP
