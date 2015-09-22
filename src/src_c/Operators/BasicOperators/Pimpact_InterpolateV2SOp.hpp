#pragma once
#ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
#define PIMPACT_INTERPOLATEVTOSOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"
#include "Pimpact_GridCoordinatesLocal.hpp"




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


template< class S,class O, int d, int dimNC>
class Space;

template<class SpaceT>
class ScalarField;


/// \brief Interpolation operator.
/// \ingroup BaseOperator
/// \ingroup Space
///
/// is used in the \c ScalarField::write method to interpolate the velocity to the pressure points, also used in \c ConvectionVOp
/// \todo check boundaries some fixing necessary
template< class Scalar, class Ordinal, int dimension, int dimNC >
class InterpolateV2S {

public:

  typedef Space<Scalar,Ordinal,dimension,dimNC> SpaceT;

  typedef ScalarField< SpaceT > DomainFieldT;
  typedef ScalarField< SpaceT > RangeFieldT;

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  TO c_;

public:

  InterpolateV2S(
      const Teuchos::RCP<const ProcGrid<Ordinal,dimension> >&  procGrid,
      const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
      const Teuchos::RCP<const StencilWidths<dimension,dimNC> >& fieldSpace,
			const Teuchos::RCP<const DomainSize<Scalar> >& domainSize,
			const Teuchos::RCP<const BoundaryConditionsLocal>& boundaryConditionsLocal,
      const Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> >& coordinatesLocal ) {


    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( gridSizeLocal->get(i) + 1 )*( fieldSpace->getDU(i) - fieldSpace->getDL(i) + 1);

      c_[i] = new Scalar[ nTemp ];
      if( i<domainSize->getDim() )
        FD_getDiffCoeff(
            procGrid->getRank(),
            gridSizeLocal->get(i),
            fieldSpace->getBL(i),
            fieldSpace->getBU(i),
            fieldSpace->getDL(i),
            fieldSpace->getDU(i),
            boundaryConditionsLocal->getBCL(i),
            boundaryConditionsLocal->getBCU(i),
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
        x.getType() == S,
        std::logic_error,
        "Pimpact::InterpolateV2S:: can only interpolate from VectorField!!!\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
        y.getType() != S,
        std::logic_error,
        "Pimpact::InterpolateV2S:: can only interpolate to Scalar!!!\n");

    auto space = x.space();

    int m = (int)x.getType();

    x.exchange( m );

    OP_interpolateV2S(
        m+1,
        space->nLoc(),
        space->bl(),
        space->bu(),
        space->dl(m),
        space->du(m),
        space->sInd(S),
        space->eInd(S),
        getC((ECoord)m),
        x.getConstRawPtr(),
        y.getRawPtr() );

    y.changed();

  }

  void assignField( const RangeFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
//    out << " --- InterpolateV2S stencil: ---";
//    for( int i=0; i<3; ++i ) {
//      out << "\ndir: " << i << "\n( ";
//      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->du(i) - space_->dl(i) + 1);
//      for( int j=0; j<nTemp; ++j )
//        out << c_[i][j] <<"\t";
//      out << ")\n";
    }

  const Scalar* getC( const ECoord& dir ) const  {
      return( c_[(int)dir] );
  }

	const std::string getLabel() const { return( "InterpolateV2S" ); };

};



/// \relates InterpolateV2S
template< class S, class O, int d, int dimNC >
Teuchos::RCP<const InterpolateV2S<S,O,d,dimNC> > createInterpolateV2S(
    const Teuchos::RCP<const ProcGrid<O,d> >&  procGrid,
    const Teuchos::RCP<const GridSizeLocal<O,d> >& gridSizeLocal,
    const Teuchos::RCP<const StencilWidths<d,dimNC> >& fieldSpace,
		const Teuchos::RCP<const DomainSize<S> >& domainSize,
		const Teuchos::RCP<const BoundaryConditionsLocal>& boundaryConditionsLocal,
    const Teuchos::RCP<const GridCoordinatesLocal<S,O,d> >& coordinatesLocal ) {

  return(
      Teuchos::rcp(
          new InterpolateV2S<S,O,d,dimNC>(
              procGrid,
              gridSizeLocal,
              fieldSpace,
              domainSize,
							boundaryConditionsLocal,
              coordinatesLocal ) ) );
}



/// \relates InterpolateV2S
/// \todo make specalization of create<Inter>( space)
template< class S, class O, int d, int dimNC >
Teuchos::RCP<const InterpolateV2S<S,O,d,dimNC> > createInterpolateV2S(
    const Teuchos::RCP<const Space<S,O,d,dimNC> >& space ) {
  return( space->getInterpolateV2S() );
}


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::InterpolateV2S<double,int,3,2>;
extern template class Pimpact::InterpolateV2S<double,int,3,4>;
extern template class Pimpact::InterpolateV2S<double,int,4,2>;
extern template class Pimpact::InterpolateV2S<double,int,4,4>;
#endif


#endif // end of #ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
