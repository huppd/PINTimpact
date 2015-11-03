#pragma once
#ifndef PIMPACT_GRIDSIZELOCAL_HPP
#define PIMPACT_GRIDSIZELOCAL_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{



/// \brief local grid size
///
/// generated from \c GridSizeGlobal and \c ProcGridSize
/// \f$ nLoc = (nGlo-1)/nProc + 1 \f$
/// \ingroup SpaceObject
/// \todo think about inheriting from Tuple or generalize for use Global/Local
template<class Ordinal, int dim>
class GridSizeLocal : public Teuchos::Tuple<Ordinal,dim> {

  template< class OT, int dT, int dNC >
  friend Teuchos::RCP<const GridSizeLocal<OT,dT> > createGridSizeLocal(
      const Teuchos::RCP<const GridSizeGlobal<OT> >& gsg,
      const Teuchos::RCP<const ProcGrid<OT,dT> >& pg,
      const Teuchos::RCP<const StencilWidths<dT,dNC> >& sW );

protected:

	/// \todo change for pasp
	/// \param stencilWidths only necessary for size checking, and defining gs[3]
	template<int dNC>
	GridSizeLocal(
			const Teuchos::RCP<const GridSizeGlobal<Ordinal> >& gridSizeGlobal,
			const Teuchos::RCP<const ProcGrid      <Ordinal,dim> >& procGrid,
			const Teuchos::RCP<const StencilWidths     <dim,dNC> >& stencilWidths ):
		Teuchos::Tuple<Ordinal,dim>() {

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
			for( int i=0; i<dim; ++i ) {
				TEUCHOS_TEST_FOR_EXCEPT( (*this)[i] < stencilWidths->getBL(i) );
				TEUCHOS_TEST_FOR_EXCEPT( (*this)[i] < stencilWidths->getBU(i) );
			}
    // problem is used for MultiGridSpace
		//      for( int i=0; i<2; ++i )
		//        if( (gridSize_[i]-1)!=1 )
		//          TEUCHOS_TEST_FOR_EXCEPTION(
		//              (gridSize_[i]-1)%2 != 0,
		//              std::logic_error,
		//              "!!!ERROR! GridSizeLocal: "<< gridSize_[i] << " cannot be used for multigrid!!!\n" );
  }


public:

	const Ordinal& get( int i ) const { return( (*this)[i] ); }

  void print( std::ostream& out=std::cout ) const {
    out << " \tlocal grid size= " << (*this) << "\n";
  };


}; // end of class GridSizeLocal



/// \brief creates GridSizeLocal
/// \relates GridSizeLocal
template< class O, int d, int dNC>
Teuchos::RCP<const GridSizeLocal<O,d> > createGridSizeLocal(
    const Teuchos::RCP<const GridSizeGlobal<O> >& gsg,
    const Teuchos::RCP<const ProcGrid<O,d> >& pg,
    const Teuchos::RCP<const StencilWidths<d,dNC> >& sW ) {
  return(
      Teuchos::rcp(
          new GridSizeLocal<O,d>( gsg, pg, sW ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZELOCAL_HPP
