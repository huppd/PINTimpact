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




/// \brief local grid size
/// generated from \c GridSizeGlobal and \c ProcGridSize
/// \ingroup Space
template< class Ordinal=int, int dim=3 >
class GridSizeLocal {

  template< class OT, int dT >
  friend Teuchos::RCP<const GridSizeLocal<OT,dT> > createGridSizeLocal(
      const Teuchos::RCP<const GridSizeGlobal<OT,dT> >& gsg,
      const Teuchos::RCP<const ProcGridSize<OT,dT> >& pgs );

  template< class OT, int dT >
  friend Teuchos::RCP<const GridSizeLocal<OT,dT> > createGridSizeLocal();

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;

protected:

  TO gridSize_;


  GridSizeLocal( TO gridSize ): gridSize_( gridSize ) {
    test();
  }

  GridSizeLocal(
      const Teuchos::RCP<const GridSizeGlobal<Ordinal,dim> >& gridSizeGlobal,
      const Teuchos::RCP<const ProcGridSize  <Ordinal,dim> >& procGridSize ):
        gridSize_() {

    for( int i=0; i<dim; ++i )
      TEUCHOS_TEST_FOR_EXCEPTION(
          gridSizeGlobal->get(i) < procGridSize->get(i),
          std::logic_error,
          "!!!ERROR!GridSizeLocal!!!\n" );

    for( int i=0; i<3; ++i )
      gridSize_[i] = 1 + ( gridSizeGlobal->get(i)-1 )/procGridSize->get(i);

    if( 4==dim )
      gridSize_[3] = gridSizeGlobal->get(3)/procGridSize->get(3);

    test();
  };


public:

  const Ordinal& get( int i ) const {
    return( gridSize_[i] );
  }

  const Ordinal* get() const {
    return( gridSize_.getRawPtr() );
  }


  void set_Impact() const {
    fsetLS(
        gridSize_[0],
        gridSize_[1],
        gridSize_[2] );
  };

  void print( std::ostream& out=std::cout ) const {
    out << " \tnx=" << gridSize_[0] ;
    out << " \tny=" << gridSize_[1] ;
    out << " \tnz=" << gridSize_[2] ;
    if( 4==dim ) out << "\tnt=" << gridSize_[3];
    out << "\n";
  };

protected:

  void test() {
    for( int i=0; i<2; ++i )
      TEUCHOS_TEST_FOR_EXCEPTION(
          gridSize_[i] < 3,
          std::logic_error,
          "!!!ERROR!GridSizeLocal!!!\n" );
    for( int i=0; i<2; ++i )
      TEUCHOS_TEST_FOR_EXCEPTION(
          (gridSize_[i]-1)%2 != 0,
          std::logic_error,
          "!!!ERROR! GridSizeLocal: cannot be used for multigrid!!!\n" );
  }

}; // end of class GridSizeLocal



/// \brief creates GridSizeLocal, and sets Impact
/// \relates GridSizeLocal
template< class O=int, int d=3 >
Teuchos::RCP<const GridSizeLocal<O,d> > createGridSizeLocal(
    const Teuchos::RCP<const GridSizeGlobal<O,d> >& gsg,
    const Teuchos::RCP<const ProcGridSize<O,d> >& pgs ) {
  return(
      Teuchos::rcp(
          new GridSizeLocal<O,d>( gsg, pgs ) ) );
}


/// \brief creates GridSizeLocal from Impact
/// \relates GridSizeLocal
template< class O=int, int d=3 >
Teuchos::RCP<const GridSizeLocal<O,d> > createGridSizeLocal() {

  Teuchos::Tuple<O,d> bla;
  SVS_get_nLoc(bla[0],bla[1],bla[2]);
  if( 4==d ) bla[3] = 1;
  return(
      Teuchos::rcp(
          new GridSizeLocal<O,d>(bla) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZELOCAL_HPP
