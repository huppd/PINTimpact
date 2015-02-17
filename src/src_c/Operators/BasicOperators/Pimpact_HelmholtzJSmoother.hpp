#pragma once
#ifndef PIMPACT_HELMHOLTZOP_HPP
#define PIMPACT_HELMHOLTZOP_HPP

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"



namespace Pimpact{


extern "C" {
void OP_helmholtz(
    const int& dimens,
    const int* const N,
    const int* const bl,
    const int* const bu,
    const int* const ss,
    const int* const nn,
    const double* const c11,
    const double* const c22,
    const double* const c33,
    const double& mulI,
    const double& multL,
    const double* const phi,
    double* const lap );
}



/// \brief Helmholtz operator
///
/// computes \f$ y = ( mulI_ I - mulL_ \Delta) x
/// \ingroup BaseOperator
template<class ST>
class HelmholtzOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef const Teuchos::Tuple<Scalar*,3> TO;

protected:

  Teuchos::RCP<const SpaceT> space_;

  Scalar mulI_;
  Scalar mulL_;

  TO cS_;
  TO cV_;

public:

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

  HelmholtzOp(
      const Teuchos::RCP<const SpaceT>& space,
      Scalar mulI=1.,
      Scalar mulL=1. ):
        space_(space),
        mulI_(mulI),
        mulL_(mulL) {

    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->bu(i) - space_->bl(i) + 1);

      cS_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->bl(i),
            space_->bu(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            int(EField::S)+1,
            i+1,
            2,
            0,
            true,
            space_->getStencilWidths()->getDimNcbC(i),
            space_->getStencilWidths()->getNcbC(i),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            cS_[i] );

      cV_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->bl(i),
            space_->bu(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            1,
            i+1,
            2,
            0,
            true,
            space_->getStencilWidths()->getDimNcbC(i),
            space_->getStencilWidths()->getNcbC(i),
            space_->getCoordinatesLocal()->getX( i, i ),
            space_->getCoordinatesLocal()->getX( i, i ),
            cV_[i] );
    }

  };

  ~HelmholtzOp() {
    for( int i=0; i<3; ++i ) {
      delete[] cS_[i];
      delete[] cV_[i];
    }
  }

//  void setMulI(Scalar mulI){ mulI_ = mulI;};
//  void setMulL(Scalar mulL){ mulL_ = mulL;};
//
//  Scalar getMulI() const { return(mulI_); };
//  Scalar getMulL() const { return(mulL_); };


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    x.exchange();

    OP_helmholtz(
        space_->dim(),
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->sInd(U),
        space_->eInd(U),
        cV_[0],
        cS_[1],
        cS_[2],
        mulI_,
        mulL_,
        x.vecC(U),
        y.vec (U) );
    OP_helmholtz(
        space_->dim(),
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->sInd(V),
        space_->eInd(V),
        cS_[0],
        cV_[1],
        cS_[2],
        mulI_,
        mulL_,
        x.vecC(V),
        y.vec (V) );
    if( 3==space_->dim() )
      OP_helmholtz(
          space_->dim(),
          space_->nLoc(),
          space_->bl(),
          space_->bu(),
          space_->sInd(W),
          space_->eInd(W),
          cS_[0],
          cS_[1],
          cV_[2],
          mulI_,
          mulL_,
          x.vecC(W),
          y.vec (W) );
    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return( space_; ) };

  void print( std::ostream& out=std::cout ) const {
    out << " --- scalar stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ndir: " << i << "\n";
      Ordinal nTemp1 = ( space_->nLoc(i) + 1 );
      Ordinal nTemp2 = ( space_->bu(i) - space_->bl(i) + 1 );
      for( int j=0; j<nTemp1; ++j ) {
        out << "\ni: " << j << "\t(";
        for( int k=0; k<nTemp2; ++k ) {
          out << cS_[i][k+nTemp2*j] <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }
    out << " --- velocity stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ndir: " << i << "\n";
      Ordinal nTemp1 = ( space_->nLoc(i) + 1 );
      Ordinal nTemp2 = ( space_->bu(i) - space_->bl(i) + 1 );
      for( int j=0; j<nTemp1; ++j ) {
        out << "\ni: " << j << "\t(";
        for( int k=0; k<nTemp2; ++k ) {
          out << cV_[i][k+nTemp2*j] <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }
  }

}; // end of class HelmholtzOp




} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_HELMHOLTZOP_HPP
