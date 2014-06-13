#pragma once
#ifndef PIMPACT_GRIDSIZELOCAL_HPP
#define PIMPACT_GRIDSIZELOCAL_HPP

#include<ostream>

#include"Teuchos_RCP.hpp"
#include"Teuchos_Tuple.hpp"

#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_Types.hpp"



extern "C" {

void fsetLS(
    const int& n1,
    const int& n2,
    const int& n3 );

void SVS_get_nLoc(int&,int&,int&);

}



namespace Pimpact{

/// \todo merge with Space should make nGlow in FieldSpace obsolete, maybe rename
template< class Ordinal=int, int dim=3 >
class GridSizeLocal {

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;

protected:

  TO gridSize_;

public:

  GridSizeLocal( TO gridSize ): gridSize_( gridSize ) {
    test();
  }

  GridSizeLocal(
      const Teuchos::RCP< GridSizeGlobal<Ordinal,dim> >& gridSizeGlobal,
      const Teuchos::RCP< ProcGridSize  <Ordinal,dim> >& procGridSize ):
        gridSize_() {

    for( int i=0; i<dim; ++i )
      gridSize_[i] = 1+(gridSizeGlobal->get(i)-1)/procGridSize->get(i);
    test();
    set_Impact();
  };


  Ordinal get( int i ) const {
    return( gridSize_[i] );
  }

  const Ordinal* const getPtr() const {
    return( gridSize_.getRawPtr() );
  }


  void set_Impact(){
    fsetLS(
        gridSize_[0],
        gridSize_[1],
        gridSize_[2] );
  };

  void print( std::ostream& out=std::cout ) {
    out << " \tnx=" << gridSize_[0] ;
    out << " \tny=" << gridSize_[1] ;
    out << " \tnz=" << gridSize_[2] ;
    if( 4==dim ) out << "\tnt=" << gridSize_[3];
    out << "\n";
  };

protected:

  void test() {
    for( int i=0; i<2; ++i ) {
      if( gridSize_[i] < 3) {
        std::cout << "!!!ERROR!GridSizeLocal: gridSize_["<<i<<"] < 3 !!!\n";
        return;
      }
      if( (gridSize_[i]-1)%2 != 0) {
        std::cout << "!!!ERROR! GridSizeLocal: (gridSize_["<<i<<"]-1)%2 != 0: cannot be used for multigrid!!!\n";
        return;
      }
    }
  }

}; // end of class GridSizeLocal


/// \relates GridSizeLocal
template< class O=int, int d=3 >
Teuchos::RCP<GridSizeLocal<O,d> > createGridSizeLocal(
    const Teuchos::RCP< GridSizeGlobal<O,d> >& gsg,
    const Teuchos::RCP< ProcGridSize<O,d> >& pgs ) {
  return(
      Teuchos::rcp(
          new GridSizeLocal<O,d>( gsg, pgs ) ) );
}


/// \relates GridSizeLocal
template< class O=int, int d=3 >
Teuchos::RCP<GridSizeLocal<O,d> > createGridSizeLocal() {

  Teuchos::Tuple<O,d> bla;
  SVS_get_nLoc(bla[0],bla[1],bla[2]);
  if( 4==d ) bla[3] = 1;
  return(
      Teuchos::rcp(
          new GridSizeLocal<O,d>(bla) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZELOCAL_HPP
