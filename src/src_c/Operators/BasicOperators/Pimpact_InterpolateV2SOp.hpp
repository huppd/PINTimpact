#pragma once
#ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
#define PIMPACT_INTERPOLATEVTOSOP_HPP


#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_CoordinatesLocal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{



extern "C"
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

	using SpaceT = Space<Scalar,Ordinal,dimension,dimNC>;

	using DomainFieldT = ScalarField< SpaceT >;
	using RangeFieldT = ScalarField< SpaceT >;

protected:

	using TO = const Teuchos::Tuple< Teuchos::ArrayRCP<Scalar>, 3 >;

	TO c_;
	TO cm_;

public:

	InterpolateV2S(
			const Teuchos::RCP<const IndexSpace<Ordinal,dimension> >&  indexSpace,
			const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
			const Teuchos::RCP<const StencilWidths<dimension,dimNC> >& stencilWidths,
			const Teuchos::RCP<const DomainSize<Scalar> >& domainSize,
			const Teuchos::RCP<const BoundaryConditionsLocal<dimension> >& boundaryConditionsLocal,
			const Teuchos::RCP<const CoordinatesLocal<Scalar,Ordinal,dimension,dimNC> >& coordinatesLocal ) {


		for( int i=0; i<3; ++i ) {

			Ordinal nTemp = ( gridSizeLocal->get(i) + 1 )*( stencilWidths->getDU(i) - stencilWidths->getDL(i) + 1);
			//c_[i] = new Scalar[ nTemp ];
			c_[i] = Teuchos::arcp<Scalar>( nTemp );
			cm_[i] = Teuchos::arcp<Scalar>( nTemp );

			if( i<domainSize->getDim() ) {
				FD_getDiffCoeff(
						gridSizeLocal->get(i),
						stencilWidths->getBL(i),
						stencilWidths->getBU(i),
						stencilWidths->getDL(i),
						stencilWidths->getDU(i),
						boundaryConditionsLocal->getBCL(i),
						boundaryConditionsLocal->getBCU(i),
						indexSpace->getShift(i),
						3,
						i+1,  // direction
						0,    // 0-derivative
						0,    // central
						true, // not working with stretching mapping
						stencilWidths->getDimNcbD(i),
						stencilWidths->getNcbD(i),
						coordinatesLocal->getX( i, i ),
						coordinatesLocal->getX( i, EField::S ),
						cm_[i].get() );
				FD_getDiffCoeff(
						gridSizeLocal->get(i),
						stencilWidths->getBL(i),
						stencilWidths->getBU(i),
						stencilWidths->getDL(i),
						stencilWidths->getDU(i),
						boundaryConditionsLocal->getBCL(i),
						boundaryConditionsLocal->getBCU(i),
						indexSpace->getShift(i),
						3,
						i+1,  // direction
						0,    // 0-derivative
						0,    // central
						false, // mapping, works with interpolateV2S
						stencilWidths->getDimNcbD(i),
						stencilWidths->getNcbD(i),
						coordinatesLocal->getX( i, i ),
						coordinatesLocal->getX( i, EField::S ),
						c_[i].get() );
			}
		}
	};


	//~InterpolateV2S() {
		//for( int i=0; i<3; ++i ) {
			//delete[] c_[i];
		//}
	//}


	void apply( const DomainFieldT& x, RangeFieldT& y, Belos::ETrans trans=Belos::NOTRANS ) const {

		TEUCHOS_TEST_FOR_EXCEPTION(
				x.getType() == S,
				std::logic_error,
				"Pimpact::InterpolateV2S:: can only interpolate from VectorField!!!\n");

		TEUCHOS_TEST_FOR_EXCEPTION(
				y.getType() != S,
				std::logic_error,
				"Pimpact::InterpolateV2S:: can only interpolate to Scalar!!!\n");

		Teuchos::RCP<const SpaceT> space = x.space();

		int m = static_cast<int>( x.getType() );

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
				getC( static_cast<ECoord>(m) ),
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
		return( c_[static_cast<int>(dir)].get() );
	}
	const Scalar* getCM( const ECoord& dir ) const  {
		return( cm_[static_cast<int>(dir)].get() );
	}

	const std::string getLabel() const { return( "InterpolateV2S" ); };

	};



	/// \relates InterpolateV2S
	template< class S, class O, int d, int dimNC >
		Teuchos::RCP<const InterpolateV2S<S,O,d,dimNC> > createInterpolateV2S(
				const Teuchos::RCP<const IndexSpace<O,d> >&  iS,
				const Teuchos::RCP<const GridSizeLocal<O,d> >& gridSizeLocal,
				const Teuchos::RCP<const StencilWidths<d,dimNC> >& stencilWidths,
				const Teuchos::RCP<const DomainSize<S> >& domainSize,
				const Teuchos::RCP<const BoundaryConditionsLocal<d> >& boundaryConditionsLocal,
				const Teuchos::RCP<const CoordinatesLocal<S,O,d,dimNC> >& coordinatesLocal ) {

			return(
					Teuchos::rcp(
						new InterpolateV2S<S,O,d,dimNC>(
							iS,
							gridSizeLocal,
							stencilWidths,
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
