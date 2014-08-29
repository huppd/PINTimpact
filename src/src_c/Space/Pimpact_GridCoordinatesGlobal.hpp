#pragma once
#ifndef PIMPACT_GRIDCOORDINATESGLOBAL_HPP
#define PIMPACT_GRIDCOORDINATESGLOBAL_HPP

#include<ostream>

#include"Teuchos_RCP.hpp"
#include"Teuchos_Tuple.hpp"

#include"Pimpact_Types.hpp"




namespace Pimpact{



extern "C" {

void PI_getGlobalCoordinates(
    const int& dimens,
    const double* const L,
    const int* const M,
    const double* const y_origin,
    double* const y1p,
    double* const y2p,
    double* const y3p,
    double* const y1u,
    double* const y2v,
    double* const y3w,
    double* const dy1p,
    double* const dy2p,
    double* const dy3p,
    double* const dy1u,
    double* const dy2v,
    double* const dy3w );


}

/// \brief global grid coordinates( independent of Field Type)
/// \ingroup Space
template<class Scalar=double, class Ordinal=int, int dim=3 >
class GridCoordinatesGlobal {

//  template<class OT,int dT>
//  friend Teuchos::RCP<GridCoordinatesGlobal<OT,dT> > createGridCoordinatesGlobal();

//  template<class OT,int dT>
//  friend Teuchos::RCP<GridCoordinatesGlobal<OT,dT> > createGridCoordinatesGlobal( OT n1, OT n2, OT n3, OT nt=1 );

  template<class ST,class OT,int dT>
  friend Teuchos::RCP<GridCoordinatesGlobal<ST,OT,dT> > createGridCoordinatesGlobal(
      const Teuchos::RCP< GridSizeGlobal<OT,dT> >& gridSize,
      const Teuchos::RCP< DomainSize<ST> >& domainSize );

public:

  typedef const Teuchos::Tuple<Scalar*,dim> TO;

protected:

  Teuchos::RCP< GridSizeGlobal<Ordinal,dim> > gridSize_;

  TO xS_;
  TO xV_;

  TO dxS_;
  TO dxV_;

  GridCoordinatesGlobal(
      const Teuchos::RCP< GridSizeGlobal<Ordinal,dim> >& gridSize,
      const Teuchos::RCP< DomainSize<Scalar> >& domainSize ):
    gridSize_( gridSize ) {

    for( int i=0; i<dim; ++i ) {
      xS_[i] = new Scalar[ gridSize_->getSizeP()[i] ];
      xV_[i] = new Scalar[ gridSize_->getSizeP()[i] + 1 ];
      dxS_[i] = new Scalar[ gridSize_->getSizeP()[i] ];
      dxV_[i] = new Scalar[ gridSize_->getSizeP()[i] + 1 ];
    }
    PI_getGlobalCoordinates(
        domainSize->getDim(),
        domainSize->getSizeP(),
        gridSize_->getSizeP(),
        domainSize->getOriginP(),
        xS_[X],
        xS_[Y],
        xS_[Z],
        xV_[X],
        xV_[Y],
        xV_[Z],
        dxS_[X],
        dxS_[Y],
        dxS_[Z],
        dxV_[X],
        dxV_[Y],
        dxV_[Z] );
  };

public:

  ~GridCoordinatesGlobal() {
    for( int i=0; i<dim; ++i ) {
      delete[] xS_[i];
      delete[] xV_[i] ;
      delete[] dxS_[i];
      delete[] dxV_[i];
    }
  };

  enum EField { U=0, V=1, W=2, S=4 };

  enum ECoord { X=0, Y=1, Z=2, T=4 };

  const Scalar* getX( ECoord dir, EField ftype ) const  {
    if( EField::S==ftype )
      return( xS_[dir] );
    else if( (int)dir==(int)ftype )
      return( xV_[dir] );
    else
      return( xS_[dir] );
  }


  void print( std::ostream& out=std::cout ) {
    for( int i=0; i<dim; ++i ) {
      out << "ScalarField dir: " << i << ":\n(";
      for( int j=0; j<gridSize_->getSizeP()[i]; ++j )
        out << xS_[i][j] << "\t";
      out << ")\n";
    }
    for( int i=0; i<dim; ++i ) {
      out << "VectorField dir: " << i << ":\n(";
      for( int j=0; j<gridSize_->getSizeP()[i]+1; ++j )
        out << xV_[i][j] << "\t";
      out << ")\n";
    }
  };

}; // end of class GridCoordinatesGlobal


///// \brief create GridSize Global from Impact
///// \relates GridCoordinatesGlobal
//template< class O=int, int d=3 >
//Teuchos::RCP<GridCoordinatesGlobal<O,d> > createGridCoordinatesGlobal() {
//
//  Teuchos::Tuple<O,d> bla;
//  SVS_get_nGlo(bla[0],bla[1],bla[2]);
//  if( 4==d ) bla[3] = 2;
//
//  return(
//      Teuchos::rcp(
//          new GridCoordinatesGlobal<O,d>( bla ) ) );
//}


///// \brief create GridSize Global
///// \relates GridCoordinatesGlobal
//template< class O=int, int d=3 >
//Teuchos::RCP<GridCoordinatesGlobal<O,d> > createGridCoordinatesGlobal( O n1, O n2, O n3, O nt=1 ) {
//  Teuchos::Tuple<O,d> temp;
//  if( 3==d ) {
//    temp[0] = n1;
//    temp[1] = n2;
//    temp[2] = n3;
//  }
//  if( 4==d ) {
//    temp[0] = n1;
//    temp[1] = n2;
//    temp[2] = n3;
//    temp[3] = nt;
//  }
//  return(
//      Teuchos::rcp(
//          new GridCoordinatesGlobal<O,d>( temp ) ) );
//}


/// \brief create Grid coordinates Global
/// \relates GridCoordinatesGlobal
template<class S=double, class O=int, int d=3 >
Teuchos::RCP<GridCoordinatesGlobal<S,O,d> > createGridCoordinatesGlobal(
    const Teuchos::RCP< GridSizeGlobal<O,d> >& gridSize,
    const Teuchos::RCP< DomainSize<S> >& domainSize ) {

  return(
      Teuchos::rcp(
          new GridCoordinatesGlobal<S,O,d>( gridSize, domainSize ) ) );

}





} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDCOORDINATESGLOBAL_HPP
