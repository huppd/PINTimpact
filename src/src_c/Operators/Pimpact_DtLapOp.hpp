#pragma once
#ifndef PIMPACT_DTLAPOP_HPP
#define PIMPACT_DTLAPOP_HPP


#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"

#include "Pimpact_HelmholtzOp.hpp"



namespace Pimpact{


extern "C" {
  void OP_DtHelmholtz(
      const int& dimens,
      const int* const N,
      const int* const bL,
      const int* const bU,
      const int& m,
      const double& mulI,
      const double& multL,
      double* const phiI,
      double* const phiL,
      double* const Lap);
}


/// \ingroup ModeOperator
template<class Scalar,class Ordinal>
class DtLapOp {

  Scalar alpha2_;
  Scalar iRe_;
  Teuchos::RCP<Helmholtz<Scalar,Ordinal> > L_;

public:

  DtLapOp( Scalar alpha2=1., Scalar iRe=1. ):
    alpha2_(alpha2),
    iRe_(iRe),
    L_(Teuchos::rcp( new Helmholtz<Scalar,Ordinal>(0,iRe) ) ) {};

  DtLapOp( Scalar alpha2, const Teuchos::RCP<Helmholtz<Scalar,Ordinal> >& L ):
    alpha2_(alpha2),
    iRe_(L->getMulL()),
    L_(L) {};

  typedef ModeField<VectorField<Scalar,Ordinal> >  DomainFieldT;
  typedef ModeField<VectorField<Scalar,Ordinal> >  RangeFieldT;


  /// \f[ \begin{bmatrix} \hat{\mathbf{y}}^c \\ \hat{\mathbf{y}}^s  \end{bmatrix}
  /// = \begin{bmatrix} mulL \mathbf{L} & \alpha^2 k \mathbf{I} \\ -\alpha^2 k \mathbf{I} & mulL \mathbf{L} \end{bmatrix}
  /// \begin{bmatrix} \hat{\mathbf{x}}^c \\ \hat{\mathbf{x}}^s  \end{bmatrix} \f]
  void apply(const DomainFieldT& x, RangeFieldT& y, int k=1 ) const {
//    L_->apply( x.getConstCField(), y.getCField() );
//    y.getCFieldPtr()->add( 1., y.getConstCField(), alpha2_*k, x.getConstSField() );
//    L_->apply( x.getConstSField(), y.getSField() );
//    y.getSFieldPtr()->add( -alpha2_*k, x.getConstCField(), 1., y.getConstSField() );

    const Ordinal& dim = y.getConstCField().dim();

    for( int vel_dir=0; vel_dir<dim; ++vel_dir ) {

         for( int dir=0; dir<dim; ++dir ) {
           if( !x.getConstCField().is_exchanged(vel_dir,dir) ) x.getConstCField().exchange( vel_dir, dir );
           if( !x.getConstSField().is_exchanged(vel_dir,dir) ) x.getConstSField().exchange( vel_dir, dir );
         }

         OP_DtHelmholtz(
             dim,
             y.getConstCField().nLoc(),
             y.getConstCField().bl(),
             y.getConstCField().bu(),
             vel_dir+1,
             alpha2_*k,
             iRe_,
             x.getConstSField().vec_[vel_dir],
             x.getConstCField().vec_[vel_dir],
             y.getCField().vec_[vel_dir] ) ;

         OP_DtHelmholtz(
             dim,
             y.getConstCField().nLoc(),
             y.getConstCField().bl(),
             y.getConstCField().bu(),
             vel_dir+1,
             -alpha2_*k,
             iRe_,
             x.getConstCField().vec_[vel_dir],
             x.getConstSField().vec_[vel_dir],
             y.getSField().vec_[vel_dir] ) ;
       }
       y.getCField().changed();
       y.getSField().changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  Teuchos::RCP< Helmholtz<Scalar,Ordinal> > getInnerOpPtr() {
    return( L_ );
  }

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DtL



/// \relates DtL
template< class Scalar, class Ordinal>
Teuchos::RCP< DtLapOp<Scalar,Ordinal> > createDtLapOp( Scalar alpha2=1., Scalar iRe=1. ) {
  return( Teuchos::rcp( new DtLapOp<Scalar,Ordinal>( alpha2, iRe ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DTLAPOP_HPP
