#pragma once
#ifndef PIMPACT_NONLINEARJACOBIANOP_HPP
#define PIMPACT_NONLINEARJACOBIANOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_FieldFactory.hpp"




extern "C" {

void OP_nonlinear(
      double* const phi1U, double* const phi1V, double* const phi1W,
      double* const phi2U, double* const phi2V, double* const phi2W,
      double* const nl1,   double* const nl2,   double* const nl3,
      const double& mul );

}



namespace Pimpact {



/// \ingroup BaseOperator
/// \note u_ has to contain appropriate BC, temp_ and y doesnot matter, x should have zero BC
template<class Scalar,class Ordinal, int dimension=3>
class NonlinearJacobian {

public:

  typedef VectorField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal,dimension>  RangeFieldT;

protected:

  Teuchos::RCP<DomainFieldT> u_;

  const bool isNewton_;

public:

  NonlinearJacobian( const bool& isNewton=true ):
    u_(Teuchos::null),
    isNewton_(isNewton) {};

  NonlinearJacobian( const Teuchos::RCP<DomainFieldT>& u, const bool& isNewton=true ):
    u_( u ),
    isNewton_(isNewton) {};

  void assignField( const DomainFieldT& mv ) {

    if( u_.is_null() )
      u_ = mv.clone();
    else
      u_->assign( mv );

    u_->exchange();

  };

  void apply( const DomainFieldT& x, RangeFieldT& y, Scalar mul=0. ) const {

    x.exchange();

    if( std::abs(mul) < 1.e-12 ) {
      y.init(0.);
      mul = 1.;
    }

    OP_nonlinear(
          u_->vec(0), u_->vec(1), u_->vec(2),
          x. vecC(0), x. vecC(1), x. vecC(2),
          y.  vec(0), y.  vec(1), y.  vec(2),
//        u_->vec_[0], u_->vec_[1], u_->vec_[2],
//        x.  vec_[0], x.  vec_[1], x.  vec_[2],
//        y.  vec_[0], y.  vec_[1], y.  vec_[2],
        mul );

    if( isNewton_ ) {
      OP_nonlinear(
          x. vecC(0), x. vecC(1), x. vecC(2),
          u_->vec(0), u_->vec(1), u_->vec(2),
          y.  vec(0), y.  vec(1), y.  vec(2),
          mul );
    }

    y.changed();

  }

  bool hasApplyTranspose() const { return( false ); }


}; // end of class NonlinearJacobian



/// \relates NonlinearJacobian
template< class S=double, class O=int, int d=3>
Teuchos::RCP<NonlinearJacobian<S,O,d> > createNonlinearJacobian(
    const Teuchos::RCP<typename NonlinearJacobian<S,O,d>::DomainFieldT>& u = Teuchos::null,
    const bool& isNewton=true ) {

  if( Teuchos::is_null(u) )
    return( Teuchos::rcp( new NonlinearJacobian<S,O,d>( isNewton ) ) );
  else
    return( Teuchos::rcp( new NonlinearJacobian<S,O,d>( u, isNewton ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_NONLINEARJACOBIANOP_HPP
