#pragma once
#ifndef PIMPACT_RESTRICTIONOP_HPP
#define PIMPACT_RESTRICTIONOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Teuchos_TestForException.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {


extern "C" {

void MG_getCRS(
    const int& S,
    const int& N,
    const int& BC_L,
    const int& BC_U,
    double* const cR );

void MG_getCRV(
    const int& N,
    const int& bL,
    const int& bU,
    const int& SS,
    const int& NN,
    const int& BC_L,
    const int& BC_U,
    const double* const xs,
    const double* const xv,
    double* const cRV );

void MG_restrict(
    const int& dimens,
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
    const double* const cR1,
    const double* const cR2,
    const double* const cR3,
    const double* const phif,
    double* const phic );

void MG_restrictV(
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
    const double* const cRV,
    const double* const phif,
    double* const phic );

}



template< class Scalar=double, class Ordinal=int, int dimension=3 >
class RestrictionOp {

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  RangeFieldT;

  typedef Space<Scalar,Ordinal,dimension> SpaceT;

  Teuchos::RCP<const SpaceT> spaceF_;
  Teuchos::RCP<const SpaceT> spaceC_;

  Teuchos::Tuple<Scalar*,3> cRS_;

  Teuchos::Tuple<Scalar*,3> cRV_;

public:

  RestrictionOp(
      const Teuchos::RCP<const SpaceT>& spaceF,
      const Teuchos::RCP<const SpaceT>& spaceC ):
        spaceF_(spaceF),spaceC_(spaceC) {

    for( int i=0; i<3; ++i ) {
      cRS_[i] = new Scalar[ 3*(spaceC_->eInd(EField::S)[i]-spaceC_->sInd(EField::S)[i]+1) ];
      MG_getCRS(
          spaceC_->sInd(EField::S)[i],
          spaceC_->eInd(EField::S)[i],
          spaceC_->getDomain()->getBCLocal()->getBCL(i),
          spaceC_->getDomain()->getBCLocal()->getBCU(i),
          cRS_[i] );
      cRV_[i] = new Scalar[ 2*( spaceC_->eIndB(i)[i]-spaceC_->sIndB(i)[i]+1 ) ];
      MG_getCRV(
          spaceC_->getGridSizeLocal()->get(i),
          spaceC_->bl(i),
          spaceC_->bu(i),
          spaceC_->sIndB(i)[i],
          spaceC_->eIndB(i)[i],
          spaceC_->getDomain()->getBCLocal()->getBCL(i),
          spaceC_->getDomain()->getBCLocal()->getBCU(i),
          spaceC_->getCoordinatesLocal()->getX( i, EField::S ),
          spaceC_->getCoordinatesLocal()->getX( i, i ),
          cRV_[i] );
    }

  }
  ~RestrictionOp() {
    for( int i=0; i<3; ++i ) {
      delete[] cRS_[i];
      delete[] cRV_[i];
    }
  }

  void apply( const DomainFieldT& x, RangeFieldT& y ) {

    TEUCHOS_TEST_FOR_EXCEPTION( x.fType_!=y.fType_     , std::logic_error, "Error!!! has to be of same FieldType!!!");

    if( EField::S==x.fType_ ) {
      x.exchange();
      MG_restrict(
          x.dim(),
          x.nLoc(),
          x.bl(),
          x.bu(),
          x.sInd(),
          x.eInd(),
          y.nLoc(),
          y.bl(),
          y.bu(),
          y.sInd(),
          y.eInd(),
          cRS_[0],
          cRS_[1],
          cRS_[2],
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
      MG_restrictV(
          x.dim(),
          dir+1,
          x.nLoc(),
          x.bl(),
          x.bu(),
          x.sIndB(),
          x.eIndB(),
          y.nLoc(),
          y.bl(),
          y.bu(),
          y.sIndB(),
          y.eIndB(),
          cRV_[dir],
          x.s_,
          y.s_ );
      y.changed();
    }
  }

  void print(  std::ostream& out=std::cout ) const {
    for( int j=0; j<3; ++j ) {
      out << "\n Scalar dir: " << j << ":\n";
      for( int i=0; i<3*( spaceC_->eInd(EField::S)[j]-spaceC_->sInd(EField::S)[j]+1 ); ++i)
        out << cRS_[j][i] << "\t";
    }
    out << "\n";
    for( int j=0; j<3; ++j ) {
      out << "\n Vector dir: " << j << ":\n";
      for( int i=0; i<2*( spaceC_->eIndB(j)[j]-spaceC_->sIndB(j)[j]+1 ); ++i)
        out << cRV_[j][i] << "\t";
    }
    out << "\n";
  }

}; // end of class RestrictionOp



template<class S=double, class O=int, int d=3>
Teuchos::RCP< RestrictionOp<S,O,d> > createRestrictionOp(
    const Teuchos::RCP<const Space<S,O,d> >& spaceF,
    const Teuchos::RCP<const Space<S,O,d> >& spaceC ) {

  return( Teuchos::rcp( new RestrictionOp<S,O,d>(spaceF,spaceC) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_RESTRICTIONOP_HPP
