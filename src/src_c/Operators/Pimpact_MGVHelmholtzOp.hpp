#pragma once
#define PIMPACT_MGVHELMHOLTZOP_HPP
#ifndef PIMPACT_MGVHELMHOLTZOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"



namespace Pimpact{


extern "C" {

void OP_MGV(
    const int& m,
    const double& mulL,
    const bool& init_yes,
    const int& gstart,
    double* const bb,
    double* const phi,
    const int& problem_type );

}

/// \brief Multigrid V cycle operator
///
/// uses bicgstab, and multigrid from Impact
/// \ingroup BaseOperator
/// \deprecated
template<class Scalar,class Ordinal, int dimension=3>
class MGVHelmholtzOp {

protected:

  Scalar mulI_;
  Scalar mulL_;

  bool initYes_;

public:

  typedef VectorField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal,dimension>  RangeFieldT;

  MGVHelmholtzOp(
      Scalar mulI=0.,
      Scalar mulL=1.,
      bool initYes=true
      ):
    mulI_(mulI),
    mulL_(mulL),
    initYes_(initYes)
  {};

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

      for( int i=0; i<1; ++i )
        OP_MGV(
            vel_dir+1,
            mulL_/mulI_,
            false,
            1,
            x.vecC(vel_dir),
            y.vec (vel_dir),
            1 );
    }
    y.scale( mulI_ );
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MGVHelmholtzOp




/// \relates MGVHelmholtzOp
template<class S,class O, int d=3>
Teuchos::RCP<MGVHelmholtzOp<S,O,d> > createMGVHelmholtzOp(
    S mulI=0.,
    S mulL=1.,
//    Scalar tol=1.e-3,
//    int nItMax=100,
    bool initYes=true
//    bool quietYes1=true,
//    bool quietYes2=true
    ) {

  return(
      Teuchos::rcp( new MGVHelmholtzOp<S,O,d>(mulI, mulL, initYes ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MGVHELMHOLTZOP_HPP
