#pragma once
#ifndef PIMPACT_NONLINEAROP_HPP
#define PIMPACT_NONLINEAROP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"



namespace Pimpact {


extern "C" {

  void OP_nonlinear( //const bool& exch_yes,
      double* phi1U, double* phi1V, double* phi1W,
      double* phi2U, double* phi2V, double* phi2W,
      double* nl1,   double* nl2,   double* nl3 );

}


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
//    u_->assign( mv );
//    u_->exchange();
  };

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    if( Teuchos::is_null(u_) )
      apply( x, x, y);
    else
      apply( *u_, x, y );

  }

  void apply( const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z) const {

    for( int vel_dir=0; vel_dir<x.dim(); ++vel_dir )
      x.exchange( vel_dir, vel_dir );
    y.exchange();

    OP_nonlinear( //false,
        x.vec_[0],x.vec_[1],x.vec_[2],
        y.vec_[0],y.vec_[1],y.vec_[2],
        z.vec_[0],z.vec_[1],z.vec_[2] );

    z.changed();

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
