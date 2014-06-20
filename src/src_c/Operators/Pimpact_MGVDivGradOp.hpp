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
template<class Scalar=double,class Ordinal=int, int dimension=3>
class MGVDivGradOp {

protected:

  bool initYes_;

public:

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  RangeFieldT;

  MGVDivGradOp(
      bool initYes=true
  ):
    initYes_(initYes)
  {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    if( initYes_ )
      y.init( 0. );
    else
      y.exchange();

    x.exchange();

    for( int i=0; i<1; ++i )
    {
      OP_MGV(
          0,
          1.,
          false,
          1,
          x.s_,
          y.s_,
          2 ); /// not sure what is the difference between 2,4,5
      SF_level( y.s_ );
    }

    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MGVDivGradOp




/// \relates MGVDivGradOp
template<class S=double,class O=int,int d=3>
Teuchos::RCP<MGVDivGradOp<S,O,d> > createMGVDivGradOp(
    bool initYes=true ) {

  return(
      Teuchos::rcp( new MGVDivGradOp<S,O,d>( initYes ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MGVDIVGRADOP_HPP
