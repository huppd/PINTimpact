#pragma once
#ifndef PIMPACT_TIMENONLINEARJACOBIANOP_HPP
#define PIMPACT_TIMENONLINEARJACOBIANOP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_NonlinearOp.hpp"
#include "Pimpact_TimeField.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {


/// \ingroup TimeOperator
/// \todo better design would be to wrap ConvectionJacobianOp(no because then you would store stencil mutliple times, so one would have to abstract Jacobian)
/// \note u_ has to contain appropriate BC, temp_ and y doesnot matter, x should have zero BC
/// \todo add expcetion for space::dimension!=4
/// \deprecated
template<class ST,bool CrankNicolsonYes=false >
class TimeNonlinearJacobian {

public:

  using SpaceT = ST;

  using Scalar = typename SpaceT::Scalar;

  using DomainFieldT = TimeField< VectorField<SpaceT> >;
  using RangeFieldT = TimeField< VectorField<SpaceT> >;

protected:

  Teuchos::RCP<DomainFieldT> u_;

  const bool isNewton_;

  Teuchos::RCP< VectorField<SpaceT> > temp_;

  Teuchos::RCP<const NonlinearOp<SpaceT> > op_;

public:

  TimeNonlinearJacobian(
      const Teuchos::RCP<const SpaceT>& space,
      const bool& isNewton=true ):
    u_( Teuchos::null ),
    isNewton_(isNewton),
    temp_( create<Pimpact::VectorField>(space) ),
    op_( createNonlinearOp<SpaceT>(space) )
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

        if( isNewton_ ) {
          op_->apply( **i, **k, *temp_, 1. );
        }
        if( j>=y.sInd_ )
          (*j)->add( 1., **j, 0.5, *temp_ );
        ++j;
        if( j<y.eInd_ ) {
          (*j)->add( 0., **j, 0.5, *temp_ );
          (*j)->changed();
        }
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

        if( isNewton_ ) {
          op_->apply( **i, **k, **j, 1. );
        }
        (*j)->changed();
        ++j; ++k;
      }
    }

    y.changed();

  }

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "TimeNonlinearJacobian " ); };

}; // end of class TimeNonlinearJacobian



/// \relates TimeNonlinearJacobian
template< class SpaceT, bool CNY=false >
Teuchos::RCP<TimeNonlinearJacobian<SpaceT,CNY> >
createTimeNonlinearJacobian(
    const Teuchos::RCP<const SpaceT>& space,
    const bool& isNewton=true ) {

    return( Teuchos::rcp( new TimeNonlinearJacobian<SpaceT,CNY>( space, isNewton ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_TIMENONLINEARJACOBIANOP_HPP
