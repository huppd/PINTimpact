#pragma once
#ifndef PIMPACT_NONLINEAROP_HPP
#define PIMPACT_NONLINEAROP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"

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


/// \ingroup BaseOperator
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

//  void setU(const Teuchos::RCP<DomainFieldT>& u) { u_=u; }
//  void assignU(Teuchos::RCP<DomainFieldT> u) { u_->assign(u); }

  void assignField( const DomainFieldT& mv ) {
    u_->assign( mv );
  };

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

  void apply( const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z) const {
//    u_.assing( x );
    if( Teuchos::is_null(u_) )
      OP_nonlinear( true,
          x.vec_[0],x.vec_[1],x.vec_[2],
          y.vec_[0],y.vec_[1],y.vec_[2],
          z.vec_[0],z.vec_[1],z.vec_[2] );
    return;
  }
  bool hasApplyTranspose() const { return( false ); }

}; // end of class Nonlinear



/// \relates Nonlinear
template< class S, class O >
Teuchos::RCP<Nonlinear<S,O> > createNonlinear() {
  return( Teuchos::rcp( new Nonlinear<S,O>() ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_NONLINEAROP_HPP
