#pragma once
#ifndef PIMPACT_GRIDSIZELOCAL_HPP
#define PIMPACT_GRIDSIZELOCAL_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief local grid size
///
/// \tparam OrdinalT
/// \tparam dim as soon as Time is own class ->sdim
/// generated from \c GridSizeGlobal and \c ProcGridSize
/// \f$ nLoc = (nGlo-1)/nProc + 1 \f$
/// \ingroup SpaceObject
template<class OrdinalT, int sd, int dim>
class GridSizeLocal : public Teuchos::Tuple<OrdinalT,dim> {

  template< class OT, int sdT, int dT, int dNC >
  friend Teuchos::RCP<const GridSizeLocal<OT,sdT,dT> > createGridSizeLocal(
    const Teuchos::RCP<const GridSizeGlobal<OT,sdT> >& gsg,
    const Teuchos::RCP<const ProcGrid<OT,dT> >& pg,
    const Teuchos::RCP<const StencilWidths<dT,dNC> >& sW );

protected:

  /// \brief constructor
  ///
  /// \tparam dNC stencil width
  /// \param gridSizeGlobal global grid size
  /// \param procGrid processor grid
  /// \param stencilWidths only necessary for size checking, and defining gs[3]
  template<int dNC>
  GridSizeLocal(
    const Teuchos::RCP<const GridSizeGlobal<OrdinalT,sd> >& gridSizeGlobal,
    const Teuchos::RCP<const ProcGrid      <OrdinalT,dim> >& procGrid,
    const Teuchos::RCP<const StencilWidths     <dim,dNC> >& stencilWidths ):
    Teuchos::Tuple<OrdinalT,dim>() {

    for( int i=0; i<3; ++i )
      TEUCHOS_TEST_FOR_EXCEPT( gridSizeGlobal->get(i) < procGrid->getNP(i) );
    if( 4==dim )
      TEUCHOS_TEST_FOR_EXCEPT( (gridSizeGlobal->get(3)+1) < procGrid->getNP(3) );

    for( int i=0; i<3; ++i )
      (*this)[i] = 1 + ( gridSizeGlobal->get(i)-1 )/procGrid->getNP(i);

    if( 4==dim ) {
      if( stencilWidths->spectralT() )
        (*this)[3] = gridSizeGlobal->get(3);
      else
        (*this)[3] = gridSizeGlobal->get(3)/procGrid->getNP(3);
    }

    // tests if local grid size is usable for multigrid and is big enough
    for( int i=0; i<sd; ++i ) {
      TEUCHOS_TEST_FOR_EXCEPT( (*this)[i] < stencilWidths->getBL(i) );
      TEUCHOS_TEST_FOR_EXCEPT( (*this)[i] < stencilWidths->getBU(i) );
    }
    // problem is used for MultiGridSpace
    //      for( int i=0; i<2; ++i )
    //        if( (gridSize_[i]-1)!=1 )
    //          assert( gridSize_[i]-1)%2 == 0 );
  }


public:

  constexpr const OrdinalT& get( const int i ) const {
    return (*this)[i];
  }

  void print( std::ostream& out=std::cout ) const {
    out << " \tlocal grid size= " << (*this) << "\n";
  };


}; // end of class GridSizeLocal



/// \brief creates GridSizeLocal
/// \relates GridSizeLocal
template< class O, int sd, int d, int dNC>
Teuchos::RCP<const GridSizeLocal<O,sd,d> > createGridSizeLocal(
    const Teuchos::RCP<const GridSizeGlobal<O,sd> >& gsg,
    const Teuchos::RCP<const ProcGrid<O,d> >& pg,
    const Teuchos::RCP<const StencilWidths<d,dNC> >& sW ) {

  return Teuchos::rcp( new GridSizeLocal<O,sd,d>( gsg, pg, sW ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZELOCAL_HPP
