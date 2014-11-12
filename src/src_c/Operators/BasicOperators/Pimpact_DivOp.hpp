#pragma once
#ifndef PIMPACT_DIVOP_HPP
#define PIMPACT_DIVOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

namespace Pimpact{


extern "C" {

  void OP_div(
      const int& dimens,
      const int* const N,
      const int* const bl,
      const int* const bu,
      const int* const dl,
      const int* const du,
      const int* const ss,
      const int* const nn,
      const double* const c1,
      const double* const c2,
      const double* const c3,
      const double* const phiU,
      const double* const phiV,
      const double* const phiW,
      double* const lap );

}


/// \brief Divergence operator.
/// \ingroup BaseOperator
template<class ST>
class DivOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  Teuchos::RCP<const SpaceT> space_;

  TO c_;

public:

  DivOp( const Teuchos::RCP<const SpaceT>& space ):
    space_(space) {

    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->du(i) - space_->dl(i) + 1);

      c_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->dl(i),
            space_->du(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            3,
            i+1,
            1,
            0,
            true,
            space_->getStencilWidths()->getDimNcbD(i),
            space_->getStencilWidths()->getNcbD(i),
            space_->getCoordinatesLocal()->getX( i, i ),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            c_[i] );
    }
  };


  ~DivOp() {
    for( int i=0; i<3; ++i ) {
      delete[] c_[i];
    }
  }


  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {

    for( int dir=0; dir<space_->dim(); ++dir )
      x.exchange( dir, dir );

    OP_div(
        space_->dim(),
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->dl(),
        space_->du(),
        space_->sInd(S),
        space_->eInd(S),
        c_[0],
        c_[1],
        c_[2],
        x.vecC(U),
        x.vecC(V),
        x.vecC(W),
        y.s_ );

    y.changed();

  }

  void assignField( const RangeFieldT& mv ) const {};
  void assignField( const DomainFieldT& mv ) const {};

  bool hasApplyTranspose() const { return( false ); }


  void print( std::ostream& out=std::cout ) const {
    out << " --- stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ndir: " << i << "\n";
      Ordinal nTemp1 = ( space_->nLoc(i) + 1 );
      Ordinal nTemp2 = ( space_->du(i) - space_->dl(i) + 1);
      for( int j=0; j<nTemp1; ++j ) {
        out << "\ni: " << j << "\t(";
        for( int k=0; k<nTemp2; ++k ) {
          out << c_[i][k+nTemp2*j] <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }
  }

};



/// \relates DivOp
template<class SpaceT>
Teuchos::RCP< DivOp<SpaceT> > createDivOp(
    const Teuchos::RCP<const SpaceT>& space ) {
  return( Teuchos::rcp( new DivOp<SpaceT>( space ) ) );
}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVOP_HPP
