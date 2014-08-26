#pragma once
#ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
#define PIMPACT_GRIDSIZEGLOBAL_HPP

#include<ostream>

#include"Teuchos_RCP.hpp"
#include"Teuchos_Tuple.hpp"

#include"Pimpact_Types.hpp"

extern "C" {

void fsetGS(const int& n1, const int& n2, const int& n3 );
void SVS_get_nGlo(int&,int&,int&);

}

namespace Pimpact{


/// \brief global grid size( independent of Field Type)
/// \ingroup Space
template< class Ordinal=int, int dim=3 >
class GridSizeGlobal {

  template<class OT,int dT>
  friend Teuchos::RCP<GridSizeGlobal<OT,dT> > createGridSizeGlobal();

  template<class OT,int dT>
  friend Teuchos::RCP<GridSizeGlobal<OT,dT> > createGridSizeGlobal( OT n1, OT n2, OT n3, OT nt=1 );

  template<class OT,int dT>
  friend Teuchos::RCP<GridSizeGlobal<OT,dT> > createGridSizeGlobal( const Teuchos::Tuple<OT,dT>& tuple );

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;

protected:

  TO gridSize_;


  GridSizeGlobal( TO gridSize ):
    gridSize_( gridSize ) {};

public:

  void set_Impact(){
    fsetGS( gridSize_[0], gridSize_[1], gridSize_[2] );
  };

  Ordinal get( int i ) const  {
    return( gridSize_[i] );
  }

  const Ordinal* getPtr() const {
    return( gridSize_.getRawPtr() );
  }

  TO getTuple() const {
    return( gridSize_ );
  }

  void print( std::ostream& out=std::cout ) {
    out << " \tnx=" << gridSize_[0];
    out << " \tny=" << gridSize_[1];
    out << " \tnz=" << gridSize_[2];
    if( 4==dim )  out << "\tnt=" << gridSize_[3];
    out << "\n";
  };

}; // end of class GridSizeGlobal


/// \brief create GridSize Global from Impact
/// \relates GridSizeGlobal
template< class O=int, int d=3 >
Teuchos::RCP<GridSizeGlobal<O,d> > createGridSizeGlobal() {

  Teuchos::Tuple<O,d> bla;
  SVS_get_nGlo(bla[0],bla[1],bla[2]);
  if( 4==d ) bla[3] = 2;

  return(
      Teuchos::rcp(
          new GridSizeGlobal<O,d>( bla ) ) );
}


/// \brief create GridSize Global
/// \relates GridSizeGlobal
template< class O=int, int d=3 >
Teuchos::RCP<GridSizeGlobal<O,d> > createGridSizeGlobal( O n1, O n2, O n3, O nt=1 ) {
  Teuchos::Tuple<O,d> temp;
  if( 3==d ) {
    temp[0] = n1;
    temp[1] = n2;
    temp[2] = n3;
  }
  if( 4==d ) {
    temp[0] = n1;
    temp[1] = n2;
    temp[2] = n3;
    temp[3] = nt;
  }
  return(
      Teuchos::rcp(
          new GridSizeGlobal<O,d>( temp ) ) );
}


/// \brief create GridSize Global
/// \relates GridSizeGlobal
template< class O=int, int d=3 >
Teuchos::RCP<GridSizeGlobal<O,d> > createGridSizeGlobal( const Teuchos::Tuple<O,d>& to  ) {

  return(
      Teuchos::rcp(
          new GridSizeGlobal<O,d>( to ) ) );

}





} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
