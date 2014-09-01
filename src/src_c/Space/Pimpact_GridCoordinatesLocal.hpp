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

/// \brief global grid coordinates( independent of Field Type)
/// \ingroup Space
template<class Scalar=double, class Ordinal=int, int dim=3 >
class GridCoordinatesLocal {

  template<class ST,class OT,int dT>
  friend Teuchos::RCP<GridCoordinatesLocal<ST,OT,dT> > createGridCoordinatesLocal(
      const Teuchos::RCP< const FieldSpace<OT,dT> >& fieldSpace,
      const Teuchos::RCP< DomainSize<ST> >& domainSize,
      const Teuchos::RCP< GridSizeGlobal<OT,dT> >& gridSizeGlobal,
      const Teuchos::RCP< GridSizeLocal<OT,dT> >& gridSize,
      const Teuchos::RCP< BoundaryConditionsGlobal >& bcGlobal,
      const Teuchos::RCP< BoundaryConditionsLocal >& bcLocal,
      const Teuchos::RCP< ProcGrid<OT,dT> >& procGrid,
      const Teuchos::RCP< GridCoordinatesGlobal<ST,OT,dT> >& coordGlobal );

public:

  typedef const Teuchos::Tuple<Scalar*,dim> TO;

protected:

  Teuchos::RCP< GridSizeLocal<Ordinal,dim> > gridSize_;

  TO xS_;
  TO xV_;

  TO dxS_;
  TO dxV_;

  GridCoordinatesLocal(
      const Teuchos::RCP< const FieldSpace<Ordinal,dim> >& fieldSpace,
      const Teuchos::RCP< DomainSize<Scalar> >& domainSize,
      const Teuchos::RCP< GridSizeGlobal<Ordinal,dim> >& gridSizeGlobal,
      const Teuchos::RCP< GridSizeLocal<Ordinal,dim> >& gridSize,
      const Teuchos::RCP< BoundaryConditionsGlobal >& bcGlobal,
      const Teuchos::RCP< BoundaryConditionsLocal >& bcLocal,
      const Teuchos::RCP< ProcGrid<Ordinal,dim> >& procGrid,
      const Teuchos::RCP< GridCoordinatesGlobal<Scalar,Ordinal,dim> >& coordGlobal
      ):
    gridSize_( gridSize ) {

    for( int i=0; i<dim; ++i ) {
      Ordinal nTemp = gridSize_->get(i)+fieldSpace->getBU(i)-fieldSpace->getBL(i)+1;
      xS_[i]  = new Scalar[ nTemp ];
      xV_[i]  = new Scalar[ nTemp ];
      dxS_[i] = new Scalar[ gridSize_->get(i) ];
      dxV_[i] = new Scalar[ gridSize_->get(i)+1 ];

      if( i<3 )
        PI_getLocalCoordinates(
            domainSize->getSize(i),
            gridSizeGlobal->get(i),
            gridSize_->get(i),
            fieldSpace->getBL(i),
            fieldSpace->getBU(i),
            bcGlobal->getBCL(i),
            bcGlobal->getBCU(i),
            bcLocal->getBCL(i),
            bcLocal->getBCU(i),
            procGrid->getIB(i),
            coordGlobal->get( i, EField::S),
            coordGlobal->get( i, i),
            xS_[i],
            xV_[i],
            dxS_[i],
            dxV_[i] );
      else if( 4==i )
        PI_getLocalCoordinates(
            4.*std::atan(1.),
            gridSizeGlobal->get(i),
            gridSize_->get(i),
            fieldSpace->getBL(i),
            fieldSpace->getBU(i),
            bcGlobal->getBCL(i),
            bcGlobal->getBCU(i),
            bcLocal->getBCL(i),
            bcLocal->getBCU(i),
            procGrid->getIB(i),
            coordGlobal->get( i, EField::S),
            coordGlobal->get( i, i),
            xS_[i],
            xV_[i],
            dxS_[i],
            dxV_[i] );
    }
  };

public:

  ~GridCoordinatesLocal() {
    for( int i=0; i<dim; ++i ) {
      delete[] xS_[i];
      delete[] xV_[i] ;
      delete[] dxS_[i];
      delete[] dxV_[i];
    }
  };


  const Scalar* getX( ECoord dir, EField ftype ) const  {
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


  void print( std::ostream& out=std::cout ) {
    for( int i=0; i<dim; ++i ) {
      out << "ScalarField dir: " << i << ":\n(";
      for( int j=0; j<gridSize_->get(i); ++j )
        out << xS_[i][j] << "\t";
      out << ")\n";
    }
    for( int i=0; i<dim; ++i ) {
      out << "VectorField dir: " << i << ":\n(";
      for( int j=0; j<gridSize_->get(i)+1; ++j )
        out << xV_[i][j] << "\t";
      out << ")\n";
    }
  };

}; // end of class GridCoordinatesLocal



/// \brief create Grid coordinates Global
/// \relates GridCoordinatesLocal
template<class S=double, class O=int, int d=3 >
Teuchos::RCP<GridCoordinatesLocal<S,O,d> > createGridCoordinatesLocal(
    const Teuchos::RCP< const FieldSpace<O,d> >& fieldSpace,
    const Teuchos::RCP< DomainSize<S> >& domainSize,
    const Teuchos::RCP< GridSizeGlobal<O,d> >& gridSizeGlobal,
    const Teuchos::RCP< GridSizeLocal<O,d> >& gridSize,
    const Teuchos::RCP< BoundaryConditionsGlobal >& bcGlobal,
    const Teuchos::RCP< BoundaryConditionsLocal >& bcLocal,
    const Teuchos::RCP< ProcGrid<O,d> >& procGrid,
    const Teuchos::RCP< GridCoordinatesGlobal<S,O,d> >& coordGlobal
    ) {

  return(
      Teuchos::rcp(
          new GridCoordinatesLocal<S,O,d>(
              fieldSpace,
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
