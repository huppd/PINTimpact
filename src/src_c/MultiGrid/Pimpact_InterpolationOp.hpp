#pragma once
#ifndef PIMPACT_INTERPOLATIONOP_HPP
#define PIMPACT_INTERPOLATIONOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Teuchos_TestForException.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {


extern "C" {

void MG_getCIS(
    const int& N,
    const int& bL,
    const int& bU,
    const double* const xs,
    double* const cI );

void MG_getCIV(
    const int& Nc,
    const int& bLc,
    const int& bUc,
    const int& SSc,
    const int& NNc,
    const int& BC_L,
    const int& BC_U,
    const int& Nf,
    const int& bLf,
    const int& bUf,
    const int& SSf,
//    const int& NNf,
    const double* const xc,
    const double* const xf,
    double* const cIV );

void MG_interpolate(
    const int& dimens,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const BCL,
    const int* const BCU,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const double* const cI1,
    const double* const cI2,
    const double* const cI3,
    const double* const phic,
    double* const phif );

void MG_interpolateV(
    const int& dimens,
    const int& dir,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const SSf,
    const int* const NNf,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const SSc,
    const int* const NNc,
    const int* const BCL,
    const int* const BCU,
    const double* const cIV,
    const double* const cI1,
    const double* const cI2,
    const double* const cI3,
    const double* const phif,
    double* const phic );

}



template<class SpaceT>
class InterpolationOp {

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

//  typedef Space<Scalar,Ordinal,dimension> SpaceT;

  Teuchos::RCP<const SpaceT> spaceC_;
  Teuchos::RCP<const SpaceT> spaceF_;

  Teuchos::Tuple<Scalar*,3> cIS_;

  Teuchos::Tuple<Scalar*,3> cIV_;

public:

  InterpolationOp(
      const Teuchos::RCP<const SpaceT>& spaceC,
      const Teuchos::RCP<const SpaceT>& spaceF ):
        spaceC_(spaceC),spaceF_(spaceF) {

    for( int i=0; i<3; ++i ) {

      cIS_[i] = new Scalar[ 2*( spaceC_->nLoc(i)-1+1 ) ];
      MG_getCIS(
          spaceC_->nLoc(i),
          spaceC_->bl(i),
          spaceC_->bu(i),
          spaceC_->getCoordinatesLocal()->getX( i, EField::S ),
          cIS_[i] );

      cIV_[i] = new Scalar[ 2*( spaceF_->nLoc(i)-0+1 ) ];
//      if( i<spaceC_->dim() )
        MG_getCIV(
            spaceC_->nLoc(i),
            spaceC_->bl(i),
            spaceC_->bu(i),
            spaceC_->sInd(i)[i],
            spaceC_->eInd(i)[i],
            spaceC_->getDomain()->getBCLocal()->getBCL(i),
            spaceC_->getDomain()->getBCLocal()->getBCU(i),
            spaceF_->nLoc(i),
            spaceF_->bl(i),
            spaceF_->bu(i),
            spaceF_->sInd(i)[i],
//            spaceF_->eIndB(i)[i],
            spaceC_->getCoordinatesLocal()->getX( i, i ),
            spaceF_->getCoordinatesLocal()->getX( i, i ),
            cIV_[i] );
    }

  }
  ~InterpolationOp() {
    for( int i=0; i<3; ++i ) {
      delete[] cIS_[i];
      delete[] cIV_[i];
    }
  }

  void apply( const DomainFieldT& x, RangeFieldT& y ) {

    EField fType = x.fType_;

    TEUCHOS_TEST_FOR_EXCEPTION( fType!=y.fType_, std::logic_error, "Error!!! has to be of same FieldType!!!");

    x.exchange();

    if( EField::S==fType ) {
      MG_interpolate(
          spaceC_->dim(),
          spaceC_->nLoc(),
          spaceC_->bl(),
          spaceC_->bu(),
          spaceC_->getDomain()->getBCLocal()->getBCL(),
          spaceC_->getDomain()->getBCLocal()->getBCU(),
          spaceF_->nLoc(),
          spaceF_->bl(),
          spaceF_->bu(),
          cIS_[0],
          cIS_[1],
          cIS_[2],
          x.s_,
          y.s_ );
    }
    else {
      int dir = fType;
      MG_interpolateV(
          spaceC_->dim(),
          dir+1,
          spaceC_->nLoc(),
          spaceC_->bl(),
          spaceC_->bu(),
          spaceC_->sIndB(fType),
          spaceC_->eIndB(fType),
          spaceF_->nLoc(),
          spaceF_->bl(),
          spaceF_->bu(),
          spaceF_->sIndB(fType),
          spaceF_->eIndB(fType),
          spaceF_->getDomain()->getBCGlobal()->getBCL(),
          spaceF_->getDomain()->getBCGlobal()->getBCU(),
          cIV_[dir],
          cIS_[0],
          cIS_[1],
          cIS_[2],
          x.s_,
          y.s_ );
    }
    y.changed();

  }

  void print(  std::ostream& out=std::cout ) const {

    for( int j=0; j<3; ++j ) {
      out << "\n Scalar dir: " << j << ":\n";
      for( int i=0; i<2*( spaceC_->eInd(EField::S)[j]-spaceC_->sInd(EField::S)[j]+1 ); ++i)
        out << cIS_[j][i] << "\t";
    }
    out << "\n";
    for( int j=0; j<3; ++j ) {
      out << "\n Vector dir: " << j << ":\n";
//      for( int i=0; i<4*( spaceC_->eIndB(j)[j]-spaceC_->sIndB(j)[j]+1 ); ++i)
      for( int i=0; i<( spaceF_->nLoc(j)-0+1 ); ++i) {
        for( int k=0; k<2; ++k ) {
          out << cIV_[j][i*2+k] << ", ";
        }
        out << "\n";
      }
    }
    out << "\n";

  }

}; // end of class InterpolationOp



template<class SpaceT>
Teuchos::RCP< InterpolationOp<SpaceT> > createInterpolationOp(
    const Teuchos::RCP<const SpaceT>& spaceC,
    const Teuchos::RCP<const SpaceT>& spaceF ) {

  return( Teuchos::rcp( new InterpolationOp<SpaceT>(spaceC,spaceF) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTERPOLATIONOP_HPP
