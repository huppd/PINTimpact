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
#include "Pimpact_Utils.hpp"




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
/// \tparam ScalarT
/// \tparam OrdinalT
/// \tparam dim as soon as Time is own class ->sdim
/// \tparam dimNC
///
/// Coordinate | index
/// -----------| --------------------------
/// xS         | bl..nLoc+bu+(ib-1)*(NB-1)
/// xV         | bl..nLoc+bu+(ib-1)*(NB-1)
///
/// \todo make nice interface for getter
/// \ingroup SpaceObject
template<class ScalarT, class OrdinalT, int dim, int dimNC>
class CoordinatesLocal {

	template<class ST,class OT,int sdT,int dT, int dNC>
	friend Teuchos::RCP<const CoordinatesLocal<ST,OT,dT,dNC> > createCoordinatesLocal(
			const Teuchos::RCP<const StencilWidths<dT,dNC> >& fieldSpace,
			const Teuchos::RCP<const DomainSize<ST,sdT> >& domainSize,
			const Teuchos::RCP<const GridSizeGlobal<OT,sdT> >& gridSizeGlobal,
			const Teuchos::RCP<const GridSizeLocal<OT,sdT,dT> >& gridSizeLocal,
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
	template<int sdim>
	CoordinatesLocal(
			const Teuchos::RCP<const StencilWidths<dim,dimNC> >& stencilWidths,
			const Teuchos::RCP<const DomainSize<ScalarT,sdim> >& domainSize,
			const Teuchos::RCP<const GridSizeGlobal<OrdinalT,sdim> >& gridSizeGlobal,
			const Teuchos::RCP<const GridSizeLocal<OrdinalT,sdim,dim> >& gridSizeLocal,
			const Teuchos::RCP<const BoundaryConditionsGlobal<dim> >& bcGlobal,
			const Teuchos::RCP<const BoundaryConditionsLocal<dim> >& bcLocal,
			const Teuchos::RCP<const ProcGrid<OrdinalT,dim> >& procGrid,
			const Teuchos::RCP<const CoordinatesGlobal<ScalarT,OrdinalT,dim> >& coordGlobal ):
	stencilWidths_(stencilWidths) {

		for( int i=0; i<dim; ++i ) {

			OrdinalT nTemp = gridSizeLocal->get(i) + stencilWidths->getBU(i) - stencilWidths->getBL(i) + 1;

			xS_[i]  = Teuchos::arcp<ScalarT>( nTemp                   );
			xV_[i]  = Teuchos::arcp<ScalarT>( nTemp                   );
			dxS_[i] = Teuchos::arcp<ScalarT>( gridSizeLocal->get(i) );
			dxV_[i] = Teuchos::arcp<ScalarT>( gridSizeLocal->get(i)+1 );

			F fi = static_cast<F>( i );

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
						coordGlobal->getX( i, F::S),
						coordGlobal->getX( i, fi),
						xS_[i].getRawPtr(),
						xV_[i].getRawPtr(),
						dxS_[i].getRawPtr(),
						dxV_[i].getRawPtr() );
			else if( 3==i )
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
						coordGlobal->getX( i, F::S ),
						coordGlobal->getX( i, fi ),
						xS_[i].getRawPtr(),
						xV_[i].getRawPtr(),
						dxS_[i].getRawPtr(),
						dxV_[i].getRawPtr() );
		}
	}

public:

	/// \name getter
	/// @{ 

	constexpr const ScalarT* getX( const int& dir, const F& ftype ) const  {
		return(
				( F::S==ftype || dir!=ftype ) ?
					xS_[dir].getRawPtr() :
					xV_[dir].getRawPtr()
				);
  }

  constexpr const ScalarT& getX( const F& ftype, const int& dir, const OrdinalT& i) const  {
		return(
				( F::S==ftype || dir!=ftype )?
					xS_[dir][i-stencilWidths_->getBL(dir)]:
					xV_[dir][i-stencilWidths_->getBL(dir)]
				);
  }

	///  @} 

  void print( std::ostream& out=std::cout ) const {

    for( int dir=0; dir<dim; ++dir ) {
      out << "Local coordinates of scalars in dir: " << toString( static_cast<ECoord>(dir) ) << "\n";
      out << "i\txS\n";
			OrdinalT j = 0;
			for( typename Teuchos::ArrayRCP<ScalarT>::iterator jp=xS_[dir].begin(); jp<xS_[dir].end(); ++jp )
				out << j++ + stencilWidths_->getBL(dir) << "\t" << *jp << "\n";
		}
    for( int dir=0; dir<dim; ++dir ) {
      out << "Local coordinates of velocities in dir: " << toString( static_cast<ECoord>(dir) ) << "\n";
			OrdinalT j = 0;
			for( typename Teuchos::ArrayRCP<ScalarT>::iterator jp=xV_[dir].begin(); jp<xV_[dir].end(); ++jp )
				out << j++ + stencilWidths_->getBL(dir)<< "\t" << *jp << "\n";
		}
	};

}; // end of class CoordinatesLocal



/// \brief create Grid coordinates Global
/// \relates CoordinatesLocal
template<class ST, class OT, int sd, int d, int dNC>
Teuchos::RCP<const CoordinatesLocal<ST,OT,d,dNC> >
createCoordinatesLocal(
		const Teuchos::RCP<const StencilWidths<d,dNC> >& stencilWidths,
		const Teuchos::RCP<const DomainSize<ST,sd> >& domainSize,
		const Teuchos::RCP<const GridSizeGlobal<OT,sd> >& gridSizeGlobal,
		const Teuchos::RCP<const GridSizeLocal<OT,sd,d> >& gridSizeLocal,
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
