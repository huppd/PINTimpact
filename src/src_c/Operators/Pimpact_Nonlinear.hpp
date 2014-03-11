#pragma once
#ifndef PIMPACT_NONLINEAROP_HPP
#define PIMPACT_NONLINEAROP_HPP


#include "Pimpact_Types.hpp"

extern "C" {
//  void OP_grad( const int& m, double* phi, double *grad );
//  void OP_div( double* phiU, double* phiV, double* phiW, double* div);
//  void OP_helmholtz( const int& m, const bool& exch_yes, const double& mulI, const double& multL, double* phi, double* Lap);
//  void OP_div_grad( const bool& corner_yes, double* phi, double* lap );
  void OP_nonlinear( const bool& exch_yes,
      double* phi1U, double* phi1V, double* phi1W,
      double* phi2U, double* phi2V, double* phi2W,
      double* nl1,   double* nl2,   double* nl3 );
}


namespace Pimpact {


template<class Scalar,class Ordinal>
class Nonlinear {
public:
  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;
  typedef NonModeOp OpType;
private:
  Teuchos::RCP<DomainFieldT> u_;
public:
  Nonlinear():u_(Teuchos::null) {};

  void setU(const Teuchos::RCP<DomainFieldT>& u) { u_=u; }
  void assignU(Teuchos::RCP<DomainFieldT> u) { u_->assign(u); }

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
//    u_.assing( x );
    if( Teuchos::is_null(u_) )
      OP_nonlinear( true,
          x.vec_[0],x.vec_[1],x.vec_[2],
          x.vec_[0],x.vec_[1],x.vec_[2],
          y.vec_[0],y.vec_[1],y.vec_[2] );
    else
      OP_nonlinear( true,
          u_->vec_[0],u_->vec_[1],u_->vec_[2],
          x.vec_[0],x.vec_[1],x.vec_[2],
          y.vec_[0],y.vec_[1],y.vec_[2] );
    return;
  }
  bool hasApplyTranspose() const { return( false ); }
};


template< class S, class O >
Teuchos::RCP<Nonlinear<S,O> > createNonlinear() {
  return( Teuchos::rcp( new Nonlinear<S,O>() ) );
}

template<class Scalar,class Ordinal>
class NonlinearJacobian {
public:
  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;
  typedef NonModeOp OpType;
private:
  Teuchos::RCP<DomainFieldT> u_;
  Teuchos::RCP<DomainFieldT> temp_;
public:
  NonlinearJacobian():u_(Teuchos::null) {};

  void setU(Teuchos::RCP<DomainFieldT> u) {
    u_=u;
    temp_=u->clone(ShallowCopy);
  }
  void assignU(Teuchos::RCP<DomainFieldT> u) {
    u_->add( 0., u_, 1., u );
//    u_->assign(u);
  }

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
//    u_.assing( x );
    OP_nonlinear( true,
        u_->vec_[0],u_->vec_[1],u_->vec_[2],
        x.vec_[0],x.vec_[1],x.vec_[2],
        temp_->vec_[0],temp_->vec_[1],temp_->vec_[2] );

    OP_nonlinear( true,
        x.vec_[0],x.vec_[1],x.vec_[2],
        u_->vec_[0],u_->vec_[1],u_->vec_[2],
        y.vec_[0],y.vec_[1],y.vec_[2] );
    y.add( 1., *temp_, 1., y );

    return;
  }
  bool hasApplyTranspose() const { return( false ); }
}; // end of class NonlinearJacobian


template< class S, class O>
Teuchos::RCP<NonlinearJacobian<S,O> > createNonlinearJacobian() {
  return( Teuchos::rcp( new NonlinearJacobian<S,O>() ) );
}

//template<class Scalar,class Ordinal>
//class nonlinear_Udx {
//public:
//  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
//  typedef VectorField<Scalar,Ordinal>  RangeFieldT;
//  typedef NonModeOp OpType;
//private:
//  Teuchos::RCP<DomainFieldT> u_;
//public:
//
//  void set_U(Teuchos::RCP<DomainFieldT> u) { u_=u; }
//
//  void apply(const DomainFieldT& x, RangeFieldT& y) const { return; }
//  bool hasApplyTranspose() const { return( true ); }
//};


//template<class Scalar,class Ordinal>
//class nonlinear_dxU {
//public:
//  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
//  typedef VectorField<Scalar,Ordinal>  RangeFieldT;
//  typedef NonModeOp OpType;
//private:
//  Teuchos::RCP<DomainFieldT> u_;
//public:
//
//  void set_U(Teuchos::RCP<DomainFieldT> u) { u_=u; }
//
//  void apply(const DomainFieldT& x, RangeFieldT& y) const { return; }
//  bool hasApplyTranspose() const { return( true ); }
//};

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_NONLINEAROP_HPP
