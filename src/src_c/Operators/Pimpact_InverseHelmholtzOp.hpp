#pragma once
#define PIMPACT_INVERSEHELMHOLTZOP_HPP
#ifndef PIMPACT_INVERSEHELMHOLTZOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"



namespace Pimpact{


extern "C" {

void OP_SolveHelmholtz(
    const int& m,
    const double& mulL,
    const double& epsU,
    const int& n_it_max,
//    const bool& init_yes,
    const double* const bb,
    double* const phi,
    const bool& quiet_yes1,
    const bool& quiet_yes2 );

}

/// \brief InverseHelmholtz operator
///
/// uses bicgstab, and multigrid from Impact
/// \ingroup BaseOperator
/// \deprecated
template<class Scalar,class Ordinal>
class InverseHelmholtzOp {

protected:

  Scalar mulI_;
  Scalar mulL_;

  Scalar tol_;
  int nItMax_;

  bool initYes_;

  bool quietYes1_;
  bool quietYes2_;

public:

  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;

  InverseHelmholtzOp(
      Scalar mulI=0.,
      Scalar mulL=1.,
      Scalar tol=1.e-3,
      int nItMax=100,
      bool initYes=true,
      bool quietYes1=true,
      bool quietYes2=true ):
    mulI_(mulI),
    mulL_(mulL),
    tol_(tol),
    nItMax_(nItMax),
    initYes_(initYes),
    quietYes1_(quietYes1),
    quietYes2_(quietYes2) {};

  void setMulI(Scalar mulI){ mulI_ = mulI;};
  void setMulL(Scalar mulL){ mulL_ = mulL;};

  Scalar getMulI() const { return(mulI_); };
  Scalar getMulL() const { return(mulL_); };


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    if( initYes_ )
      y.init( 0. );

    x.exchange();
    y.exchange();

    for( int vel_dir=0; vel_dir<x.dim(); ++vel_dir ) {

//      for( int dir=0; dir<x.dim(); ++dir )
//        if( !x.is_exchanged(vel_dir,dir) ) x.exchange( vel_dir, dir );

      OP_SolveHelmholtz(
          vel_dir+1,
          mulL_/mulI_,
          tol_,
          nItMax_,
//          initYes_,
          x.vecC(vel_dir),
          y.vec (vel_dir),
          quietYes1_,
          quietYes2_ );
//      OP_InverseHelmholtzOp( vel_dir+1, mulI_, mulL_, x.vec_[vel_dir], y.vec_[vel_dir] ) ;

    }
    y.scale( mulI_ );
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class InverseHelmholtzOp




/// \relates InverseHelmholtzOp
template<class Scalar=double,class Ordinal=int>
Teuchos::RCP<InverseHelmholtzOp<Scalar,Ordinal> > createInverseHelmholtzOp(
    Scalar mulI=0.,
    Scalar mulL=1.,
    Scalar tol=1.e-3,
    int nItMax=100,
    bool initYes=true,
    bool quietYes1=true,
    bool quietYes2=true ) {

  return(
      Teuchos::rcp( new InverseHelmholtzOp<Scalar,Ordinal>(mulI, mulL, tol, nItMax, initYes, quietYes1, quietYes2 ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVERSEHELMHOLTZOP_HPP
