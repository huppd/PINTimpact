#pragma once
#ifndef PIMPACT_GRIDCOORDINATESLOCAL_HPP
#define PIMPACT_GRIDCOORDINATESLOCAL_HPP


#include<cmath>
#include<ostream>

#include"Teuchos_RCP.hpp"
#include"Teuchos_Tuple.hpp"

#include"Pimpact_Types.hpp"



namespace Pimpact{



extern "C" {

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

}

/// \brief local grid coordinates
///
/// Coordinate | index
/// -----------| --------------------------
/// xS         | bl..nLoc+bu+(ib-1)*(NB-1)
/// xV         | bl..nLoc+bu+(ib-1)*(NB-1)
///
/// \ingroup SpaceObject
template<class ScalarT, class OrdinalT, int dim>
class GridCoordinatesLocal {

  template<class ST,class OT,int dT, int dNC>
  friend Teuchos::RCP<const GridCoordinatesLocal<ST,OT,dT> > createGridCoordinatesLocal(
      const Teuchos::RCP<const StencilWidths<dT,dNC> >& fieldSpace,
      const Teuchos::RCP<const DomainSize<ST> >& domainSize,
      const Teuchos::RCP<const GridSizeGlobal<OT,dT> >& gridSizeGlobal,
      const Teuchos::RCP<const GridSizeLocal<OT,dT> >& gridSize,
      const Teuchos::RCP<const BoundaryConditionsGlobal<dT> >& bcGlobal,
      const Teuchos::RCP<const BoundaryConditionsLocal >& bcLocal,
      const Teuchos::RCP<const ProcGrid<OT,dT> >& procGrid,
      const Teuchos::RCP<const GridCoordinatesGlobal<ST,OT,dT> >& coordGlobal );

public:

  using TO = const Teuchos::Tuple<ScalarT*,dim>;

protected:

  Teuchos::RCP<const GridSizeLocal<OrdinalT,dim> > gridSize_;

  TO xS_;
  TO xV_;

  TO dxS_;
  TO dxV_;

	template<int dimNC>
		GridCoordinatesLocal(
				const Teuchos::RCP<const StencilWidths<dim,dimNC> >& stencilWidths,
				const Teuchos::RCP<const DomainSize<ScalarT> >& domainSize,
				const Teuchos::RCP<const GridSizeGlobal<OrdinalT,dim> >& gridSizeGlobal,
				const Teuchos::RCP<const GridSizeLocal<OrdinalT,dim> >& gridSize,
				const Teuchos::RCP<const BoundaryConditionsGlobal<dim> >& bcGlobal,
				const Teuchos::RCP<const BoundaryConditionsLocal >& bcLocal,
				const Teuchos::RCP<const ProcGrid<OrdinalT,dim> >& procGrid,
				const Teuchos::RCP<const GridCoordinatesGlobal<ScalarT,OrdinalT,dim> >& coordGlobal ):
			gridSize_( gridSize ) {

				for( int i=0; i<dim; ++i ) {

					if( i<3 ) {
						OrdinalT nTemp = gridSize_->get(i)+stencilWidths->getBU(i)-stencilWidths->getBL(i)+1;
						xS_[i]  = new ScalarT[ nTemp ];
						xV_[i]  = new ScalarT[ nTemp ];
						dxS_[i] = new ScalarT[ gridSize_->get(i) ];
						dxV_[i] = new ScalarT[ gridSize_->get(i)+1 ];

						PI_getLocalCoordinates(
								domainSize->getSize(i),
								gridSizeGlobal->get(i),
								gridSize_->get(i),
								stencilWidths->getBL(i),
								stencilWidths->getBU(i),
								bcGlobal->getBCL(i),
								bcGlobal->getBCU(i),
								bcLocal->getBCL(i),
								bcLocal->getBCU(i),
								procGrid->getIB(i),
								coordGlobal->get( i, EField::S ),
								coordGlobal->get( i, i),
								xS_[i],
								xV_[i],
								dxS_[i],
								dxV_[i] );
					}
					else {

						OrdinalT nTemp = gridSize_->get(i)+stencilWidths->getBU(i)-stencilWidths->getBL(i);

						xS_[i]  = new ScalarT[ nTemp ];
						xV_[i]  = new ScalarT[ nTemp ];
						dxS_[i] = new ScalarT[ gridSize_->get(i) ];
						dxV_[i] = new ScalarT[ gridSize_->get(i)+1 ];

						OrdinalT nt = gridSizeGlobal->get(i);
						ScalarT pi = 4.*std::atan(1.);
						ScalarT offset = procGrid->getShift(i) + stencilWidths->getBL(i);

						for( OrdinalT it=0; it<nTemp; ++it ) {
							xS_[i][it] = 2.*pi*( static_cast<ScalarT>(it) + offset )/nt;
							xV_[i][it] = 2.*pi*( static_cast<ScalarT>(it) + offset )/nt;
						}
					}
				}
			}

public:

  ~GridCoordinatesLocal() {
    for( int i=0; i<dim; ++i ) {
      delete[] xS_[i];
      delete[] xV_[i] ;
      delete[] dxS_[i];
      delete[] dxV_[i];
    }
  };


  constexpr const ScalarT* getX( const int& dir, const int& ftype ) const  {
		return(
				( EField::S==static_cast<EField>(ftype) || dir!=ftype )?
					xS_[dir]:
					xV_[dir]
				);
  }


  void print( std::ostream& out=std::cout ) const {
    for( int i=0; i<dim; ++i ) {
      out << "Local coordinates of scalars in dir: " << i << "\n";
      out << "i\txS\n";
      for( int j=0; j<gridSize_->get(i); ++j )
        out << j+1 << "\t" << xS_[i][j] << "\n";
    }
    for( int i=0; i<dim; ++i ) {
      out << "Local coordinates of velocities in dir: " << i << "\n";
      for( int j=0; j<gridSize_->get(i)+1; ++j )
        out << j<< "\t" << xV_[i][j] << "\n";
    }
  };

}; // end of class GridCoordinatesLocal



/// \brief create Grid coordinates Global
/// \relates GridCoordinatesLocal
template<class S, class O, int d, int dNC>
Teuchos::RCP<const GridCoordinatesLocal<S,O,d> > createGridCoordinatesLocal(
    const Teuchos::RCP<const StencilWidths<d,dNC> >& stencilWidths,
    const Teuchos::RCP<const DomainSize<S> >& domainSize,
    const Teuchos::RCP<const GridSizeGlobal<O,d> >& gridSizeGlobal,
    const Teuchos::RCP<const GridSizeLocal<O,d> >& gridSize,
    const Teuchos::RCP<const BoundaryConditionsGlobal<d> >& bcGlobal,
    const Teuchos::RCP<const BoundaryConditionsLocal >& bcLocal,
    const Teuchos::RCP<const ProcGrid<O,d> >& procGrid,
    const Teuchos::RCP<const GridCoordinatesGlobal<S,O,d> >& coordGlobal
) {

  return(
      Teuchos::rcp(
          new GridCoordinatesLocal<S,O,d>(
              stencilWidths,
              domainSize,
              gridSizeGlobal,
              gridSize,
              bcGlobal,
              bcLocal,
              procGrid,
              coordGlobal ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDCOORDINATESLOCAL_HPP
