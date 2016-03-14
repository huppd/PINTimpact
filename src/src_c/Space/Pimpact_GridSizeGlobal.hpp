#pragma once
#ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
#define PIMPACT_GRIDSIZEGLOBAL_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"




namespace Pimpact{


/// \brief global grid size(independent of FieldType)
///
/// \todo remove tpara dimension and rm getters
/// \ingroup SpaceObject
template<class OrdinalT>
class GridSizeGlobal : public Teuchos::Tuple<OrdinalT,4> {

  template<class OT>
  friend Teuchos::RCP<const GridSizeGlobal<OT> > createGridSizeGlobal( OT n1, OT n2, OT n3, OT nt );

  template<class OT>
  friend Teuchos::RCP<const GridSizeGlobal<OT> > createGridSizeGlobal( const Teuchos::Tuple<OT,4>& tuple );


protected:

	/// \todo test also thrid dimension( not on coarser grids)
	GridSizeGlobal( const Teuchos::Tuple<OrdinalT,4>& gridSize ):
		Teuchos::Tuple<OrdinalT,4>( gridSize ) {

			for( int i=0; i<3; ++i )
				TEUCHOS_TEST_FOR_EXCEPT( ((*this)[i]-1)%2 != 0 );

		};

public:

	const OrdinalT& get( int i ) const  { return( (*this)[i] ); }

	void print( std::ostream& out=std::cout ) const {
		out << " --- GridSizeGlobal: " << *this << " ---\n";
  };


}; // end of class GridSizeGlobal



/// \brief create GridSize Global
/// \relates GridSizeGlobal
template< class O=int>
Teuchos::RCP<const GridSizeGlobal<O> > createGridSizeGlobal( O n1, O n2, O n3, O nt=1 ) {
  Teuchos::Tuple<O,4> temp;
    temp[0] = n1;
    temp[1] = n2;
    temp[2] = n3;
		temp[3] = nt;
  return(
      Teuchos::rcp(
          new GridSizeGlobal<O>( temp ) ) );
}



/// \brief create GridSize Global
/// \relates GridSizeGlobal
template<class O>
Teuchos::RCP<const GridSizeGlobal<O> > createGridSizeGlobal( const Teuchos::Tuple<O,4>& to  ) {

	return(
			Teuchos::rcp(
				new GridSizeGlobal<O>( to ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
