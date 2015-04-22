#pragma once
#ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
#define PIMPACT_GRIDSIZEGLOBAL_HPP

#include<ostream>

#include"Teuchos_RCP.hpp"
#include"Teuchos_Tuple.hpp"

#include"Pimpact_Types.hpp"




//extern "C" {
//
//void fsetGS(const int& n1, const int& n2, const int& n3 );
//
//}



namespace Pimpact{


/// \brief global grid size(independent of FieldType)
///
///
/// \todo include setter method, such that it gets update with enlarging Fourier modes!!!
/// one could think about inheriting from Tuple, or generalize for global and local use
/// \ingroup Space
template<class Ordinal, int dimension>
class GridSizeGlobal {

  template<class OT, int d>
  friend Teuchos::RCP<const GridSizeGlobal<OT,d> > createGridSizeGlobal( OT n1, OT n2, OT n3, OT nt );

  template<class OT, int d>
  friend Teuchos::RCP<const GridSizeGlobal<OT,d> > createGridSizeGlobal( const Teuchos::Tuple<OT,d>& tuple );

public:

  typedef Teuchos::Tuple<Ordinal,dimension> TO;

protected:

  TO gridSize_;

  GridSizeGlobal( const TO& gridSize ):
    gridSize_( gridSize ) {};

public:

	const Ordinal& get( int i ) const  { return( gridSize_[i] ); }

	const Ordinal* get() const { return( gridSize_.getRawPtr() ); }

	const TO& getTuple() const { return( gridSize_ ); }

  void print( std::ostream& out=std::cout ) const {
    out << " --- GridSizeGlobal: " << gridSize_ << " ---\n";

  };

}; // end of class GridSizeGlobal



/// \brief create GridSize Global
/// \relates GridSizeGlobal
template< class O=int, int d>
Teuchos::RCP<const GridSizeGlobal<O,d> > createGridSizeGlobal( O n1, O n2, O n3, O nt=1 ) {
  Teuchos::Tuple<O,d> temp;
    temp[0] = n1;
    temp[1] = n2;
    temp[2] = n3;
		if( 4==d )
			temp[3] = nt;
  return(
      Teuchos::rcp(
          new GridSizeGlobal<O,d>( temp ) ) );
}


/// \brief create GridSize Global
/// \relates GridSizeGlobal
template<class O, int d>
Teuchos::RCP<const GridSizeGlobal<O,d> > createGridSizeGlobal( const Teuchos::Tuple<O,d>& to  ) {

  return(
      Teuchos::rcp(
          new GridSizeGlobal<O,d>( to )
				)
			);

}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
