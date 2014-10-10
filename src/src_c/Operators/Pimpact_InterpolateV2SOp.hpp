#pragma once
#ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
#define PIMPACT_INTERPOLATEVTOSOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

namespace Pimpact{


extern "C" {

  void OP_interpolateV2S(
      const int& m,
      const int* const N,
      const int* const bl,
      const int* const bu,
      const int* const dl,
      const int* const du,
      const int* const ss,
      const int* const nn,
      const double* const c,
      const double* const phi,
      double* const inter );

}


/// \brief Interpolation operator.
/// \ingroup BaseOperator
/// relates io+ nonlinear
template<class Scalar,class Ordinal,int dimension=3>
class InterpolateV2S {

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  Teuchos::RCP< const Space<Scalar,Ordinal,dimension> > space_;

  TO c_;

public:

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  RangeFieldT;

  InterpolateV2S( const Teuchos::RCP<const Space<Scalar,Ordinal,dimension> >& space ):
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
            0,
            0,
            true,
            space_->getFieldSpace()->getDimNcbD(i),
            space_->getFieldSpace()->getNcbD(i),
            space_->getCoordinatesLocal()->getX( i, i ),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            c_[i] );
    }
  };


  ~InterpolateV2S() {
    for( int i=0; i<3; ++i ) {
      delete[] c_[i];
    }
  }


  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {

    TEUCHOS_TEST_FOR_EXCEPTION(
           x.fType_ == S,
           std::logic_error,
           "Pimpact::InterpolateV2S:: can only interpolate from VectorField!!!\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
           y.fType_ != S,
           std::logic_error,
           "Pimpact::InterpolateV2S:: can only interpolate to Scalar!!!\n");


    int m = (int)x.fType_;

    x.exchange( m );

    OP_interpolateV2S(
        m+1,
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->dl(),
        space_->du(),
        x.sInd(),
        x.eInd(),
        c_[m],
        x.s_,
        y.s_ );


    y.changed();

  }

  void assignField( const RangeFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }


  void print( std::ostream& out=std::cout ) const {
    out << " --- stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ni: " << i << "\n( ";
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->du(i) - space_->dl(i) + 1);
      for( int j=0; j<nTemp; ++j )
        out << c_[i][j] <<"\t";
      out << ")\n";
    }
  }

};


/// \relates InterpolateV2S
template< class S, class O, int d=3 >
Teuchos::RCP< InterpolateV2S<S,O,d> > createInterpolateV2S(
    const Teuchos::RCP<const Space<S,O,d> >& space ) {
  return( Teuchos::rcp( new InterpolateV2S<S,O,d>( space ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
