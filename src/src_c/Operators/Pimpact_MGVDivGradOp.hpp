#pragma once
#ifndef PIMPACT_MGVDIVGRADOP_HPP
#define PIMPACT_MGVDIVGRADOP_HPP

#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"



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

/// \brief Multigrid V-cylce operator
///
/// uses bicgstab, and multigrid from Impact
/// \ingroup BaseOperator
template<class Scalar,class Ordinal>
class MGVDivGradOp {

protected:

  Scalar mulI_;
  Scalar mulL_;

  bool initYes_;

public:

  typedef ScalarField<Scalar,Ordinal>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal>  RangeFieldT;

  MGVDivGradOp(
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

//    for( int i=0; i<10; ++i )
      OP_MGV(
          0,
          mulL_/mulI_,
          false,
          1,
          x.s_,
          y.s_,
          3 ); /// not sure what is the difference between 2,3,5
//    y.scale( mulI_ );
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MGVDivGradOp




/// \relates MGVDivGradOp
template<class Scalar,class Ordinal>
Teuchos::RCP<MGVDivGradOp<Scalar,Ordinal> > createMGVDivGradOp(
    Scalar mulI=0.,
    Scalar mulL=1.,
    bool initYes=true
) {

  return(
      Teuchos::rcp( new MGVDivGradOp<Scalar,Ordinal>(mulI, mulL, initYes ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MGVDIVGRADOP_HPP
