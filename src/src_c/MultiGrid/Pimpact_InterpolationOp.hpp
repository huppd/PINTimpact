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
    const int& SS,
    const int& NN,
    const double* const xs,
    double* const cI );

void MG_getCIV(
    const int& N,
    const int& bL,
    const int& bU,
    const int& SS,
    const int& NN,
    const int& BC_L,
    const int& BC_U,
    const double* const xs,
    const double* const xv,
    double* const cIV );

void MG_interpolate(
    const int& dimens,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const SSc,
    const int* const NNc,
    const int* const BCL,
    const int* const BCU,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const SSf,
    const int* const NNf,
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
    const double* const cIV,
    const double* const phif,
    double* const phic );

}



template< class Scalar=double, class Ordinal=int, int dimension=3 >
class InterpolationOp {

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  RangeFieldT;

  typedef Space<Scalar,Ordinal,dimension> SpaceT;

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
//      cIS_[i] = new Scalar[ 2*(spaceC_->eIndB(EField::S)[i]-spaceC_->sIndB(EField::S)[i]+1) ];
      cIS_[i] = new Scalar[ 2*(spaceC_->nLoc(i)-spaceC_->sIndB(EField::S)[i]+1) ];
      MG_getCIS(
          spaceC_->nLoc(i),
          spaceC_->bl(i),
          spaceC_->bu(i),
          spaceC_->sIndB(EField::S)[i],
//          spaceC_->eIndB(EField::S)[i],
          spaceC_->nLoc(i),
          spaceC_->getCoordinatesLocal()->getX( i, EField::S ),
          cIS_[i] );
      cIV_[i] = new Scalar[ 2*( spaceC_->eIndB(i)[i]-spaceC_->sIndB(i)[i]+1 ) ];
//      MG_getCIV(
//          spaceC_->getGridSizeLocal()->get(i),
//          spaceC_->bl(i),
//          spaceC_->bu(i),
//          spaceC_->sIndB(i)[i],
//          spaceC_->eIndB(i)[i],
//          spaceC_->getDomain()->getBCLocal()->getBCL(i),
//          spaceC_->getDomain()->getBCLocal()->getBCU(i),
//          spaceC_->getCoordinatesLocal()->getX( i, EField::S ),
//          spaceC_->getCoordinatesLocal()->getX( i, i ),
//          cIV_[i] );
    }

  }
  ~InterpolationOp() {
    for( int i=0; i<3; ++i ) {
      delete[] cIS_[i];
      delete[] cIV_[i];
    }
  }

  void apply( const DomainFieldT& x, RangeFieldT& y ) {

    TEUCHOS_TEST_FOR_EXCEPTION( x.fType_!=y.fType_     , std::logic_error, "Error!!! has to be of same FieldType!!!");

    if( EField::S==x.fType_ ) {
      x.exchange();
      MG_interpolate(
          x.dim(),
          x.nLoc(),
          x.bl(),
          x.bu(),
          x.sInd(),
          x.eInd(),
          x.bcL(),
          x.bcU(),
          y.nLoc(),
          y.bl(),
          y.bu(),
          y.sInd(),
          y.eInd(),
          cIS_[0],
          cIS_[1],
          cIS_[2],
          x.s_,
          y.s_ );
      y.changed();
    }
    else {
      int dir = x.fType_;
      x.exchange( dir );
//      for( int i=0; i<3; ++i ) {
//          std::cout << "SSc["<<i<<"]: " <<  y.sIndB()[i] << "\n";
//          std::cout << "NNc["<<i<<"]: " <<  y.eIndB()[i] << "\n";
//      }
//      MG_interpolateV(
//          x.dim(),
//          dir+1,
//          x.nLoc(),
//          x.bl(),
//          x.bu(),
//          x.sIndB(),
//          x.eIndB(),
//          y.nLoc(),
//          y.bl(),
//          y.bu(),
//          y.sIndB(),
//          y.eIndB(),
//          cIV_[dir],
//          x.s_,
//          y.s_ );
//      y.changed();
    }
  }

  void print(  std::ostream& out=std::cout ) const {
    for( int j=0; j<3; ++j ) {
      out << "\n Scalar dir: " << j << ":\n";
      for( int i=0; i<2*( spaceC_->eInd(EField::S)[j]-spaceC_->sInd(EField::S)[j]+1 ); ++i)
        out << cIS_[j][i] << "\t";
    }
    out << "\n";
//    for( int j=0; j<3; ++j ) {
//      out << "\n Vector dir: " << j << ":\n";
//      for( int i=0; i<2*( spaceC_->eIndB(j)[j]-spaceC_->sIndB(j)[j]+1 ); ++i)
//        out << cIV_[j][i] << "\t";
//    }
//    out << "\n";
  }

}; // end of class InterpolationOp



template<class S=double, class O=int, int d=3>
Teuchos::RCP< InterpolationOp<S,O,d> > createInterpolationOp(
    const Teuchos::RCP<const Space<S,O,d> >& spaceC,
    const Teuchos::RCP<const Space<S,O,d> >& spaceF ) {

  return( Teuchos::rcp( new InterpolationOp<S,O,d>(spaceC,spaceF) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTERPOLATIONOP_HPP
