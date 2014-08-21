#pragma once
#ifndef PIMPACT_NONLINEAROP_HPP
#define PIMPACT_NONLINEAROP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"



namespace Pimpact {


extern "C" {

void OP_nonlinear(
    double* const phi1U, double* const phi1V, double* const phi1W,
    double* const phi2U, double* const phi2V, double* const phi2W,
    double* const nl1,   double* const nl2,   double* const nl3,
    const double& mul );

}


/// \ingroup BaseOperator
template<class Scalar,class Ordinal, int dimension=3>
class Nonlinear {

public:

  typedef VectorField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal,dimension>  RangeFieldT;

private:

  Teuchos::RCP<DomainFieldT> u_;

public:

  Nonlinear():u_(Teuchos::null) {};

  void assignField( const DomainFieldT& mv ) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    //    if( Teuchos::is_null(u_) )
    apply( x, x, y);
    //    else
    //      apply( *u_, x, y );

  }

  void apply( const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=0. ) const {

    for( int vel_dir=0; vel_dir<x.dim(); ++vel_dir )
      x.exchange( vel_dir, vel_dir );
    y.exchange();

    if( std::abs(mul) < 1.e-12 ) {
      z.init( 0. );
      mul = 1.;
    }
    OP_nonlinear(
//        x.vec_[0],x.vec_[1],x.vec_[2],
        x.sFields_[0]->getRawPtr(),x.sFields_[1]->getRawPtr(),x.sFields_[2]->getRawPtr(),
        y.sFields_[0]->getRawPtr(),y.sFields_[1]->getRawPtr(),y.sFields_[2]->getRawPtr(),
        z.sFields_[0]->getRawPtr(),z.sFields_[1]->getRawPtr(),z.sFields_[2]->getRawPtr(), mul );

    z.changed();

  }


  bool hasApplyTranspose() const { return( false ); }


}; // end of class Nonlinear



/// \relates Nonlinear
template< class S=double, class O=int, int d=3 >
Teuchos::RCP<Nonlinear<S,O,d> > createNonlinear() {
  return( Teuchos::rcp( new Nonlinear<S,O,d>() ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_NONLINEAROP_HPP
