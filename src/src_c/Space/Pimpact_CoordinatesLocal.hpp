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
template<class ScalarT, class OrdinalT, int dim, int dimNC>
class CoordinatesLocal {

	template<class ST,class OT,int dT, int dNC>
	friend Teuchos::RCP<const CoordinatesLocal<ST,OT,dT,dNC> > createCoordinatesLocal(
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

	Teuchos::RCP<const StencilWidths<dim,dimNC> > stencilWidths_;

	//template<int dimNC>
	CoordinatesLocal(
			const Teuchos::RCP<const StencilWidths<dim,dimNC> >& stencilWidths,
			const Teuchos::RCP<const DomainSize<ScalarT> >& domainSize,
			const Teuchos::RCP<const GridSizeGlobal<OrdinalT> >& gridSizeGlobal,
			const Teuchos::RCP<const GridSizeLocal<OrdinalT,dim> >& gridSizeLocal,
			const Teuchos::RCP<const BoundaryConditionsGlobal<dim> >& bcGlobal,
			const Teuchos::RCP<const BoundaryConditionsLocal<dim> >& bcLocal,
			const Teuchos::RCP<const ProcGrid<OrdinalT,dim> >& procGrid,
			const Teuchos::RCP<const CoordinatesGlobal<ScalarT,OrdinalT,dim> >& coordGlobal ):
	stencilWidths_(stencilWidths) {

		for( int i=0; i<dim; ++i ) {

			OrdinalT nTemp = gridSizeLocal->get(i) + stencilWidths->getBU(i) - stencilWidths->getBL(i) + 1;

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
						coordGlobal->getX( i, EField::S),
						coordGlobal->getX( i, i),
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
						coordGlobal->getX( i, EField::S ),
						coordGlobal->getX( i, i ),
						xS_[i].getRawPtr(),
						xV_[i].getRawPtr(),
						dxS_[i].getRawPtr(),
						dxV_[i].getRawPtr() );
		}
	}

public:

	/// \name getter
	/// @{ 

  inline constexpr const ScalarT* getX( const int& dir, const int& ftype ) const  {
		return(
				( EField::S==ftype || dir!=ftype ) ?
					xS_[dir].getRawPtr() :
					xV_[dir].getRawPtr()
				);
  }

  constexpr const ScalarT& getX( const int& ftype, const int& dir, const OrdinalT& i) const  {
		return(
				( EField::S==ftype || dir!=ftype )?
					xS_[dir][i-stencilWidths_->getBL(dir)]:
					xV_[dir][i-stencilWidths_->getBL(dir)]
				);
  }

	///  @} 

  void print( std::ostream& out=std::cout ) const {

    for( int i=0; i<dim; ++i ) {
      out << "Local coordinates of scalars in dir: " << i << "\n";
      out << "i\txS\n";
			OrdinalT j = 0;
			for( typename Teuchos::ArrayRCP<ScalarT>::iterator jp=xS_[i].begin(); jp<xS_[i].end(); ++jp )
				out << ++j << "\t" << *jp << "\n";
		}
    for( int i=0; i<dim; ++i ) {
      out << "Local coordinates of velocities in dir: " << i << "\n";
			OrdinalT j = 0;
			for( typename Teuchos::ArrayRCP<ScalarT>::iterator jp=xV_[i].begin(); jp<xV_[i].end(); ++jp )
				out << j++ << "\t" << *jp << "\n";
		}
	};

}; // end of class CoordinatesLocal



/// \brief create Grid coordinates Global
/// \relates CoordinatesLocal
template<class ST, class OT, int d, int dNC>
Teuchos::RCP<const CoordinatesLocal<ST,OT,d,dNC> >
createCoordinatesLocal(
		const Teuchos::RCP<const StencilWidths<d,dNC> >& stencilWidths,
		const Teuchos::RCP<const DomainSize<ST> >& domainSize,
		const Teuchos::RCP<const GridSizeGlobal<OT> >& gridSizeGlobal,
		const Teuchos::RCP<const GridSizeLocal<OT,d> >& gridSizeLocal,
		const Teuchos::RCP<const BoundaryConditionsGlobal<d> >& bcGlobal,
		const Teuchos::RCP<const BoundaryConditionsLocal<d> >& bcLocal,
		const Teuchos::RCP<const ProcGrid<OT,d> >& procGrid,
		const Teuchos::RCP<const CoordinatesGlobal<ST,OT,d> >& coordGlobal ) {

	return(
			Teuchos::rcp(
				new CoordinatesLocal<ST,OT,d,dNC>(
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
