#pragma once
#ifndef PIMPACT_COORDINATESLOCAL_HPP
#define PIMPACT_COORDINATESLOCAL_HPP


#include <cmath>
#include <ostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"
#include "Pimpact_CoordinatesGlobal.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{



extern "C" 
void PI_getLocalCoordinates(
    const double& L,
    const int& M,
    const int& N,
    const int& bL,
    const int& bU,
    const int& BC_L_global,
    const int& BC_U_global,
    const int& BC_L,
    const int& BC_U,
    const int& iB,
    const double* const ys,
    const double* const yv,
    double* const xs,
    double* const xv,
    double* const dxs,
    double* const dxv );



/// \brief local grid coordinates
///
/// Coordinate | index
/// -----------| --------------------------
/// xS         | bl..nLoc+bu+(ib-1)*(NB-1)
/// xV         | bl..nLoc+bu+(ib-1)*(NB-1)
///
/// \ingroup SpaceObject
template<class ScalarT, class Ordinal, int dim>
class CoordinatesLocal {

	template<class ST,class OT,int dT, int dNC>
	friend Teuchos::RCP<const CoordinatesLocal<ST,OT,dT> > createCoordinatesLocal(
			const Teuchos::RCP<const StencilWidths<dT,dNC> >& fieldSpace,
			const Teuchos::RCP<const DomainSize<ST> >& domainSize,
			const Teuchos::RCP<const GridSizeGlobal<OT> >& gridSizeGlobal,
			const Teuchos::RCP<const GridSizeLocal<OT,dT> >& gridSizeLocal,
			const Teuchos::RCP<const BoundaryConditionsGlobal<dT> >& bcGlobal,
			const Teuchos::RCP<const BoundaryConditionsLocal<dT> >& bcLocal,
			const Teuchos::RCP<const ProcGrid<OT,dT> >& procGrid,
			const Teuchos::RCP<const CoordinatesGlobal<ST,OT,dT> >& coordGlobal );

protected:

	using TO = const Teuchos::Tuple< Teuchos::ArrayRCP<ScalarT>, dim >;

	TO xS_;
	TO xV_;

	TO dxS_;
	TO dxV_;

	template<int dimNC>
	CoordinatesLocal(
			const Teuchos::RCP<const StencilWidths<dim,dimNC> >& stencilWidths,
			const Teuchos::RCP<const DomainSize<ScalarT> >& domainSize,
			const Teuchos::RCP<const GridSizeGlobal<Ordinal> >& gridSizeGlobal,
			const Teuchos::RCP<const GridSizeLocal<Ordinal,dim> >& gridSizeLocal,
			const Teuchos::RCP<const BoundaryConditionsGlobal<dim> >& bcGlobal,
			const Teuchos::RCP<const BoundaryConditionsLocal<dim> >& bcLocal,
			const Teuchos::RCP<const ProcGrid<Ordinal,dim> >& procGrid,
			const Teuchos::RCP<const CoordinatesGlobal<ScalarT,Ordinal,dim> >& coordGlobal ) {

		for( int i=0; i<dim; ++i ) {

			Ordinal nTemp = gridSizeLocal->get(i) + stencilWidths->getBU(i) - stencilWidths->getBL(i) + 1;

			xS_[i]  = Teuchos::arcp<ScalarT>( nTemp                   );
			xV_[i]  = Teuchos::arcp<ScalarT>( nTemp                   );
			dxS_[i] = Teuchos::arcp<ScalarT>( gridSizeLocal->get(i)   );
			dxV_[i] = Teuchos::arcp<ScalarT>( gridSizeLocal->get(i)+1 );

			if( i<3 )
				PI_getLocalCoordinates(
						domainSize->getSize(i),
						gridSizeGlobal->get(i),
						gridSizeLocal->get(i),
						stencilWidths->getBL(i),
						stencilWidths->getBU(i),
						bcGlobal->getBCL(i),
						bcGlobal->getBCU(i),
						bcLocal->getBCL(i),
						bcLocal->getBCU(i),
						procGrid->getIB(i),
						coordGlobal->get( i, EField::S),
						coordGlobal->get( i, i),
						xS_[i].getRawPtr(),
						xV_[i].getRawPtr(),
						dxS_[i].getRawPtr(),
						dxV_[i].getRawPtr() );
			else if( 3==i );
				PI_getLocalCoordinates(
						4.*std::atan(1.),
						gridSizeGlobal->get(i),
						gridSizeLocal->get(i),
						stencilWidths->getBL(i),
						stencilWidths->getBU(i),
						bcGlobal->getBCL(i),
						bcGlobal->getBCU(i),
						bcLocal->getBCL(i),
						bcLocal->getBCU(i),
						procGrid->getIB(i),
						coordGlobal->get( i, EField::S ),
						coordGlobal->get( i, i ),
						xS_[i].getRawPtr(),
						xV_[i].getRawPtr(),
						dxS_[i].getRawPtr(),
						dxV_[i].getRawPtr() );
		}
	}

public:

	/// \name getter
	/// @{ 

  const ScalarT* getX( ECoord dir, EField ftype ) const  {
    if( EField::S==ftype )
      return( xS_[dir].getRawPtr() );
    else if( (int)dir==(int)ftype )
      return( xV_[dir].getRawPtr() );
    else
      return( xS_[dir].getRawPtr() );
  }
  const ScalarT* getX( ECoord dir, int ftype ) const  {
    return( getX( dir, (EField) ftype ) );
  }
  const ScalarT* getX( int dir, EField ftype ) const  {
    return( getX( (ECoord) dir, ftype ) );
  }
  const ScalarT* getX( int dir, int ftype ) const  {
    return( getX( (ECoord) dir, (EField) ftype ) );
  }

	///  @} 

  void print( std::ostream& out=std::cout ) const {

    for( int i=0; i<dim; ++i ) {
      out << "Local coordinates of scalars in dir: " << i << "\n";
      out << "i\txS\n";
			Ordinal j = 0;
			for( typename Teuchos::ArrayRCP<ScalarT>::iterator jp=xS_[i].begin(); jp<xS_[i].end(); ++jp )
				out << ++j << "\t" << *jp << "\n";
		}
    for( int i=0; i<dim; ++i ) {
      out << "Local coordinates of velocities in dir: " << i << "\n";
			Ordinal j = 0;
			for( typename Teuchos::ArrayRCP<ScalarT>::iterator jp=xV_[i].begin(); jp<xV_[i].end(); ++jp )
				out << j++ << "\t" << *jp << "\n";
		}

	};

}; // end of class CoordinatesLocal



/// \brief create Grid coordinates Global
/// \relates CoordinatesLocal
template<class S, class O, int d, int dNC>
Teuchos::RCP<const CoordinatesLocal<S,O,d> >
createCoordinatesLocal(
		const Teuchos::RCP<const StencilWidths<d,dNC> >& stencilWidths,
		const Teuchos::RCP<const DomainSize<S> >& domainSize,
		const Teuchos::RCP<const GridSizeGlobal<O> >& gridSizeGlobal,
		const Teuchos::RCP<const GridSizeLocal<O,d> >& gridSizeLocal,
		const Teuchos::RCP<const BoundaryConditionsGlobal<d> >& bcGlobal,
		const Teuchos::RCP<const BoundaryConditionsLocal<d> >& bcLocal,
		const Teuchos::RCP<const ProcGrid<O,d> >& procGrid,
		const Teuchos::RCP<const CoordinatesGlobal<S,O,d> >& coordGlobal ) {

	return(
			Teuchos::rcp(
				new CoordinatesLocal<S,O,d>(
					stencilWidths,
					domainSize,
					gridSizeGlobal,
					gridSizeLocal,
					bcGlobal,
					bcLocal,
					procGrid,
					coordGlobal ) ) );

}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COORDINATESLOCAL_HPP
