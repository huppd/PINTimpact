#pragma once
#ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
#define PIMPACT_INTERPOLATEVTOSOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_Domain.hpp"
#include "Pimpact_GridCoordinatesLocal.hpp"


//#include "Pimpact_ScalarField.hpp"
//#include "Pimpact_VectorField.hpp"


namespace Pimpact{


extern "C" {

void OP_interpolateV2S(
    const int& m,
    const int* const N,
    const int* const bl,
    const int* const bu,
    const int&       dl,
    const int&       du,
    const int* const ss,
    const int* const nn,
    const double* const c,
    const double* const phi,
    double* const inter );

}

template<class S, class O, int d>
class ScalarField;

/// \brief Interpolation operator.
/// \ingroup BaseOperator
/// relates io+ nonlinear
/// \todo move to space
template<class Scalar,class Ordinal,int dimension=3>
class InterpolateV2S {

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  TO c_;

public:

  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  RangeFieldT;

  InterpolateV2S(
      const Teuchos::RCP<const ProcGrid<Ordinal,dimension> >&  procGrid,
      const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
      const Teuchos::RCP<const FieldSpace<Ordinal,dimension> >& fieldSpace,
      const Teuchos::RCP<const Domain<Scalar> >& domain,
      const Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> >& coordinatesLocal ) {


    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( gridSizeLocal->get(i) + 1 )*( fieldSpace->getDU(i) - fieldSpace->getDL(i) + 1);

      c_[i] = new Scalar[ nTemp ];
      if( i<domain->getDomainSize()->getDim() )
        FD_getDiffCoeff(
            procGrid->getRank(),
            gridSizeLocal->get(i),
            fieldSpace->getBL(i),
            fieldSpace->getBU(i),
            fieldSpace->getDL(i),
            fieldSpace->getDU(i),
            domain->getBCLocal()->getBCL(i),
            domain->getBCLocal()->getBCU(i),
            procGrid->getShift(i),
            3,
            i+1,
            0,
            0,
            true,
            fieldSpace->getDimNcbD(i),
            fieldSpace->getNcbD(i),
            coordinatesLocal->getX( i, i ),
            coordinatesLocal->getX( i, EField::S ),
            c_[i] );
    }
  };


  ~InterpolateV2S() {
    for( int i=0; i<3; ++i ) {
      delete[] c_[i];
    }
  }


  void apply( const DomainFieldT& x, RangeFieldT& y, Belos::ETrans trans=Belos::NOTRANS ) const {

    TEUCHOS_TEST_FOR_EXCEPTION(
        x.fType_ == S,
        std::logic_error,
        "Pimpact::InterpolateV2S:: can only interpolate from VectorField!!!\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
        y.fType_ != S,
        std::logic_error,
        "Pimpact::InterpolateV2S:: can only interpolate to Scalar!!!\n");

    auto space = x.getSpace();

    int m = (int)x.fType_;

    x.exchange( m );
//    x.exchange();

    OP_interpolateV2S(
        m+1,
        space->nLoc(),
        space->bl(),
        space->bu(),
        space->dl(m),
        space->du(m),
        y.sInd(),
        y.eInd(),
        c_[m],
        x.s_,
        y.s_ );


    y.changed();

  }

  void assignField( const RangeFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }


  void print( std::ostream& out=std::cout ) const {}

};



/// \relates InterpolateV2S
template< class S, class O, int d=3 >
Teuchos::RCP< InterpolateV2S<S,O,d> > createInterpolateV2S(
    const Teuchos::RCP<const ProcGrid<O,d> >&  procGrid,
    const Teuchos::RCP<const GridSizeLocal<O,d> >& gridSizeLocal,
    const Teuchos::RCP<const FieldSpace<O,d> >& fieldSpace,
    const Teuchos::RCP<const Domain<S> >& domain,
    const Teuchos::RCP<const GridCoordinatesLocal<S,O,d> >& coordinatesLocal ) {

  return(
      Teuchos::rcp(
          new InterpolateV2S<S,O,d>(
              procGrid,
              gridSizeLocal,
              fieldSpace,
              domain,
              coordinatesLocal ) ) );
}



template<class S,class O, int d>
class Space;



/// \relates InterpolateV2S
template< class S, class O, int d=3 >
Teuchos::RCP< InterpolateV2S<S,O,d> > createInterpolateV2S(
    const Teuchos::RCP<const Space<S,O,d> >& space ) {
  return( Teuchos::rcp( new InterpolateV2S<S,O,d>(
      space->getProcGrid(),
      space->getGridSizeLocal(),
      space->getFieldSpace(),
      space->getDomain(),
      space->getCoordinatesLocal()
  ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
