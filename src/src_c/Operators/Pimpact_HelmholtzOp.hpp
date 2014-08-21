#pragma once
#ifndef PIMPACT_HELMHOLTZOP_HPP
#define PIMPACT_HELMHOLTDOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

#include "Pimpact_OperatorMV.hpp"



namespace Pimpact{


extern "C" {

void OP_helmholtz(
    const int& m,
    const double& mulI,
    const double& multL,
    double* const phi,
    double* const Lap );

}

/// \brief HelmholtzOp operator
/// \ingroup BaseOperator
template<class Scalar,class Ordinal,int dimension=3>
class HelmholtzOp {

protected:

  Scalar mulI_;
  Scalar mulL_;

public:

  typedef VectorField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal,dimension>  RangeFieldT;

  HelmholtzOp( Scalar mulI=1., Scalar mulL=1. ):
    mulI_(mulI),mulL_(mulL) {};

  void setMulI(Scalar mulI){ mulI_ = mulI;};
  void setMulL(Scalar mulL){ mulL_ = mulL;};

  Scalar getMulI() const { return(mulI_); };
  Scalar getMulL() const { return(mulL_); };


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    for( int vel_dir=0; vel_dir<x.dim(); ++vel_dir ) {

      for( int dir=0; dir<x.dim(); ++dir )
        if( !x.is_exchanged(vel_dir,dir) ) x.exchange( vel_dir, dir );

      OP_helmholtz( vel_dir+1, mulI_, mulL_, x.sFields_[vel_dir]->getRawPtr(), y.sFields_[vel_dir]->getRawPtr() ) ;

    }
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class HelmholtzOp




/// \relates HelmholtzOp
template<class S=double,class O=int,int d=3>
Teuchos::RCP<HelmholtzOp<S,O,d> > createHelmholtzOp( S mulI=0., S mulL=1. ) {
  return(
      Teuchos::rcp( new HelmholtzOp<S,O,d>(mulI, mulL) )
  );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_HELMHOLTZOP_HPP
