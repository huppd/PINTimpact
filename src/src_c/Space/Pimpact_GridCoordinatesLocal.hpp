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
template<class Scalar, class Ordinal, int dim>
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

  using TO = const Teuchos::Tuple<Scalar*,dim>;

protected:

  Teuchos::RCP<const GridSizeLocal<Ordinal,dim> > gridSize_;

  TO xS_;
  TO xV_;

  TO dxS_;
  TO dxV_;

	template<int dimNC>
		GridCoordinatesLocal(
				const Teuchos::RCP<const StencilWidths<dim,dimNC> >& stencilWidths,
				const Teuchos::RCP<const DomainSize<Scalar> >& domainSize,
				const Teuchos::RCP<const GridSizeGlobal<Ordinal,dim> >& gridSizeGlobal,
				const Teuchos::RCP<const GridSizeLocal<Ordinal,dim> >& gridSize,
				const Teuchos::RCP<const BoundaryConditionsGlobal<dim> >& bcGlobal,
				const Teuchos::RCP<const BoundaryConditionsLocal >& bcLocal,
				const Teuchos::RCP<const ProcGrid<Ordinal,dim> >& procGrid,
				const Teuchos::RCP<const GridCoordinatesGlobal<Scalar,Ordinal,dim> >& coordGlobal ):
			gridSize_( gridSize ) {

				for( int i=0; i<dim; ++i ) {

					if( i<3 ) {
						Ordinal nTemp = gridSize_->get(i)+stencilWidths->getBU(i)-stencilWidths->getBL(i)+1;
						xS_[i]  = new Scalar[ nTemp ];
						xV_[i]  = new Scalar[ nTemp ];
						dxS_[i] = new Scalar[ gridSize_->get(i) ];
						dxV_[i] = new Scalar[ gridSize_->get(i)+1 ];

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

						Ordinal nTemp = gridSize_->get(i)+stencilWidths->getBU(i)-stencilWidths->getBL(i);

						xS_[i]  = new Scalar[ nTemp ];
						xV_[i]  = new Scalar[ nTemp ];
						dxS_[i] = new Scalar[ gridSize_->get(i) ];
						dxV_[i] = new Scalar[ gridSize_->get(i)+1 ];

						Ordinal nt = gridSizeGlobal->get(i);
						Scalar pi = 4.*std::atan(1.);
						Scalar offset = procGrid->getShift(i) + stencilWidths->getBL(i);

						for( Ordinal it=0; it<nTemp; ++it ) {
							xS_[i][it] = 2.*pi*( (Scalar)it + offset )/nt;
							xV_[i][it] = 2.*pi*( (Scalar)it + offset )/nt;
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


  const Scalar* getX( ECoord dir, EField ftype=EField::S ) const  {
    if( EField::S==ftype )
      return( xS_[dir] );
    else if( (int)dir==(int)ftype )
      return( xV_[dir] );
    else
      return( xS_[dir] );
  }
  const Scalar* getX( ECoord dir, int ftype ) const  {
    return( getX( dir, (EField) ftype ) );
  }
  const Scalar* getX( int dir, EField ftype ) const  {
    return( getX( (ECoord) dir, ftype ) );
  }
  const Scalar* getX( int dir, int ftype ) const  {
    return( getX( (ECoord) dir, (EField) ftype ) );
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
