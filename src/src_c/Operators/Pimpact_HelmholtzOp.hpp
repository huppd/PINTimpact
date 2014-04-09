#pragma once
#ifndef PIMPACT_HELMHOLTZOP_HPP
#define PIMPACT_HELMHOLTDOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

#include "Pimpact_OperatorMV.hpp"

namespace Pimpact{


extern "C" {
//  void OP_helmholtz( const int& m, const bool& exch_yes, const double& mulI, const double& multL, double* phi, double* Lap);
  void OP_helmholtz( const int& m, const double& mulI, const double& multL, double* phi, double* Lap);
}

/// \brief Helmholtz operator
/// \ingroup BaseOperator
template<class Scalar,class Ordinal>
class Helmholtz {

protected:

  Scalar mulI_;
  Scalar mulL_;

public:

  Helmholtz():mulI_(1.),mulL_(1.) {};
  Helmholtz(Scalar mulI, Scalar mulL):mulI_(mulI),mulL_(mulL) {};

  void setMulI(Scalar mulI){ mulI_ = mulI;};
  void setMulL(Scalar mulL){ mulL_ = mulL;};

  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;
  typedef NonModeOp OpType;

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    for( int vel_dir=0; vel_dir<x.dim(); ++vel_dir ) {

      for( int dir=0; dir<x.dim(); ++dir )
        if( !x.is_exchanged(vel_dir,dir) ) x.exchange( vel_dir, dir );

      OP_helmholtz( vel_dir+1, mulI_, mulL_, x.vec_[vel_dir], y.vec_[vel_dir] ) ;

    }
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class Helmholtz




/// \relates Helmholtz
template<class Scalar,class Ordinal>
Teuchos::RCP<Helmholtz<Scalar,Ordinal> > createHelmholtz( Scalar mulI=0., Scalar mulL=1. ) {
  return(
      Teuchos::rcp( new Helmholtz<Scalar,Ordinal>(mulI, mulL) )
  );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_HELMHOLTZOP_HPP
