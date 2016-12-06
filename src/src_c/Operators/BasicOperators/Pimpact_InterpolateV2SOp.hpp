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
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{




template< class S,class O, int sd, int d, int dimNC>
class Space;



template<class SpaceT>
class ScalarField;



/// \brief Interpolation operator.
/// \ingroup BaseOperator
/// \ingroup Space
///
/// is used in the \c ScalarField::write method to interpolate the velocity to the pressure points, also used in \c ConvectionVOp
template< class Scalar, class Ordinal, int sdim, int dimension, int dimNC >
class InterpolateV2S {

public:

	using SpaceT = Space<Scalar,Ordinal,sdim,dimension,dimNC>;

	using DomainFieldT = ScalarField< SpaceT >;
	using RangeFieldT = ScalarField< SpaceT >;

protected:

	using SW = StencilWidths<dimension,dimNC>;

	using Stenc = Stencil< Scalar, Ordinal, 0, SW::DL(0), SW::DU(0) >;
	using TS = const Teuchos::Tuple< Stenc, sdim >;

	TS c_;
	TS cm_;

public:

	InterpolateV2S(
			const Teuchos::RCP<const IndexSpace<Ordinal,dimension> >&  indexSpace,
			const Teuchos::RCP<const GridSizeLocal<Ordinal,sdim,dimension> >& gridSizeLocal,
			const Teuchos::RCP<const StencilWidths<dimension,dimNC> >& stencilWidths,
			const Teuchos::RCP<const DomainSize<Scalar,sdim> >& domainSize,
			const Teuchos::RCP<const BoundaryConditionsLocal<dimension> >& boundaryConditionsLocal,
			const Teuchos::RCP<const CoordinatesLocal<Scalar,Ordinal,dimension,dimNC> >& coordinatesLocal ) {


		for( int i=0; i<sdim; ++i ) {

			c_[i]  = Stenc( gridSizeLocal->get(i) );
			cm_[i] = Stenc( gridSizeLocal->get(i) );

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
					//true, // not working with stretching mapping
					false, // not working with stretching mapping
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
					//false, // mapping, works with interpolateV2S
					true, // mapping, works with interpolateV2S
					stencilWidths->getDimNcbD(i),
					stencilWidths->getNcbD(i),
					coordinatesLocal->getX( i, i ),
					coordinatesLocal->getX( i, EField::S ),
					c_[i].get() );
		}
	};


	void apply( const DomainFieldT& x, RangeFieldT& y, const Belos::ETrans&
			trans=Belos::NOTRANS ) const {

		assert( x.getType() != S );
		assert( y.getType() == S );

		Teuchos::RCP<const SpaceT> space = x.space();

		ECoord m = static_cast<ECoord>( x.getType() );
		
		x.exchange( m );

		if( X==m ) {
			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i ) {
						y.at(i,j,k) = 0.;
						for( Ordinal ii=c_[m].bl(); ii<=c_[m].bu(); ++ii )
							y.at(i,j,k) += getC( m, i, ii )*x.at(i+ii,j,k);
					}
		}

		if( Y==m ) {
			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i ) {
						y.at(i,j,k) = 0.;
						for( Ordinal jj=c_[m].bl(); jj<=c_[m].bu(); ++jj )
							y.at(i,j,k) += getC( m, j, jj )*x.at(i,j+jj,k);
					}
		}

		if( Z==m ) {
			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i ) {
						y.at(i,j,k) = 0.;
						for( Ordinal kk=c_[m].bl(); kk<=c_[m].bu(); ++kk )
							y.at(i,j,k) += getC( m, k, kk )*x.at(i,j,k+kk);
					}
		}

		y.changed();
	}


	void assignField( const RangeFieldT& mv ) {};

	bool hasApplyTranspose() const { return( false ); }

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	void print( std::ostream& out=std::cout ) const {
		out << "--- " << getLabel() << " ---\n";
		for( int dir=0; dir<sdim; ++dir ) {
			out << "\ndir: " << toString(static_cast<ECoord>(dir)) << "\n";
			c_[dir].print( out );
		}
	}

	constexpr const Scalar* getC( const int& dir ) const  {
		return( c_[dir].get() );
	}
	constexpr const Scalar* getCM( const int& dir ) const  {
		return( cm_[dir].get() );
	}
	constexpr const Scalar& getC( const int& dir, Ordinal i, Ordinal off ) const {
		return( c_[dir].at( i, off ) );
	}
	constexpr const Scalar& getCM( const int& dir, Ordinal i, Ordinal off ) const {
		return( cm_[dir].at( i, off ) );
	}

	const std::string getLabel() const { return( "InterpolateV2S" ); };

};



/// \relates InterpolateV2S
template< class S, class O, int sd,int d, int dimNC >
Teuchos::RCP<const InterpolateV2S<S,O,sd,d,dimNC> > createInterpolateV2S(
		const Teuchos::RCP<const IndexSpace<O,d> >&  iS,
		const Teuchos::RCP<const GridSizeLocal<O,sd,d> >& gridSizeLocal,
		const Teuchos::RCP<const StencilWidths<d,dimNC> >& stencilWidths,
		const Teuchos::RCP<const DomainSize<S,sd> >& domainSize,
		const Teuchos::RCP<const BoundaryConditionsLocal<d> >& boundaryConditionsLocal,
		const Teuchos::RCP<const CoordinatesLocal<S,O,d,dimNC> >& coordinatesLocal ) {

	return(
			Teuchos::rcp(
				new InterpolateV2S<S,O,sd,d,dimNC>(
					iS,
					gridSizeLocal,
					stencilWidths,
					domainSize,
					boundaryConditionsLocal,
					coordinatesLocal ) ) );
}



/// \relates InterpolateV2S
template< class S, class O, int sd, int d, int dimNC >
Teuchos::RCP<const InterpolateV2S<S,O,sd,d,dimNC> > createInterpolateV2S(
		const Teuchos::RCP<const Space<S,O,sd,d,dimNC> >& space ) {

	return( space->getInterpolateV2S() );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTERPOLATEVTOSOP_HPP
