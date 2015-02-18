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
    const int* const ss,
    const int* const nn,
    const double* const c11,
    const double* const c22,
    const double* const c33,
    const double& mulI,
    const double& multL,
    const double* const phiI,
    const double* const phiL,
    	    double* const Lap );
}


/// \ingroup ModeOperator
/// \todo move c_ from Helmholtz
template<class ST>
class DtLapOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

protected:

  Scalar alpha2_;
  Scalar iRe_;
  Teuchos::RCP<HelmholtzOp<SpaceT> > L_; // necesary for 0 mode

public:

  typedef ModeField<VectorField<SpaceT> >  DomainFieldT;
  typedef ModeField<VectorField<SpaceT> >  RangeFieldT;

  DtLapOp(
      const Teuchos::RCP<const SpaceT>& space,
      Scalar alpha2=1.,
      Scalar iRe=1.,
      const Teuchos::RCP<HelmholtzOp<SpaceT> >& L=Teuchos::null ):
        alpha2_(alpha2),
        iRe_(iRe),
        L_(L) {
    if( L_.is_null() )
      L_ = create<HelmholtzOp>( space );
  };




  /// \f[ \begin{bmatrix} \hat{\mathbf{y}}^c \\ \hat{\mathbf{y}}^s  \end{bmatrix}
  /// = \begin{bmatrix} mulL \mathbf{L} & \alpha^2 k \mathbf{I} \\ -\alpha^2 k \mathbf{I} & mulL \mathbf{L} \end{bmatrix}
  /// \begin{bmatrix} \hat{\mathbf{x}}^c \\ \hat{\mathbf{x}}^s  \end{bmatrix} \f]
  void apply(const DomainFieldT& x, RangeFieldT& y, int k=1 ) const {

    const Ordinal& dim = space()->dim();

    for( int i=0; i<dim; ++i ) {

      for( int dir=0; dir<dim; ++dir ) {
        if( !x.getConstCField().is_exchanged(i,dir) ) x.getConstCField().exchange( i, dir );
        if( !x.getConstSField().is_exchanged(i,dir) ) x.getConstSField().exchange( i, dir );
      }

      EField fType = (EField)i;

      OP_DtHelmholtz(
          dim,
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType),
          space()->eInd(fType),
          L_->getC(X,fType),
          L_->getC(Y,fType),
          L_->getC(Z,fType),
          alpha2_*k,
          iRe_,
          x.getConstSField().vecC(i),
          x.getConstCField().vecC(i),
          y.getCField().getRawPtr(i) ) ;

      OP_DtHelmholtz(
          dim,
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType),
          space()->eInd(fType),
          L_->getC(X,fType),
          L_->getC(Y,fType),
          L_->getC(Z,fType),
          -alpha2_*k,
          iRe_,
          x.getConstCField().vecC(i),
          x.getConstSField().vecC(i),
          y.getSField().getRawPtr(i) ) ;
    }

    y.getCField().changed();
    y.getSField().changed();

  }

  void assignField( const DomainFieldT& mv ) {};

	Teuchos::RCP<const SpaceT> space() const { return(L_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  Teuchos::RCP< HelmholtzOp<SpaceT> > getInnerOpPtr() {
    return( L_ );
  }

  bool hasApplyTranspose() const { return( false ); }

}; // end of class DtL



/// \relates DtLapOp
template< class SpaceT>
Teuchos::RCP< DtLapOp<SpaceT> >
createDtLapOp(
    const Teuchos::RCP<const SpaceT>& space,
    typename SpaceT::Scalar alpha2=1.,
    typename SpaceT::Scalar iRe=1. ) {
  return( Teuchos::rcp( new DtLapOp<SpaceT>( space, alpha2, iRe ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DTLAPOP_HPP
