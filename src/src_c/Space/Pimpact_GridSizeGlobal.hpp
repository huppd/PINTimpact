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


/// \todo merge with Space should make nGlow in FieldSpace obsolete, maybe rename
template< class Ordinal=int, int dim=3 >
class GridSizeGlobal {

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;

protected:

  TO gridSize_;

public:

  GridSizeGlobal( TO gridSize):
    gridSize_( gridSize ) {};

  GridSizeGlobal( Ordinal n1, Ordinal n2, Ordinal n3 ):
    gridSize_( Teuchos::tuple( n1, n2, n3) ) {
    set_Impact();
  };

  GridSizeGlobal( Ordinal n1, Ordinal n2, Ordinal n3, Ordinal nt ):
    gridSize_( Teuchos::tuple( n1, n2, n3, nt ) ) {
    set_Impact();
  };

  void set_Impact(){
    fsetGS( gridSize_[0], gridSize_[1], gridSize_[2] );
  };

  Ordinal get( int i ) const  {
    return( gridSize_[i] );
  }

  const Ordinal* const getPtr() const {
    return( gridSize_.getRawPtr() );
  }

  void print( std::ostream& out=std::cout ) {
    out << " \tnx=" << gridSize_[0];
    out << " \tny=" << gridSize_[1];
    out << " \tnz=" << gridSize_[2] << "\n";
  };

}; // end of class GridSizeGlobal


/// \relates GridSizeGlobal
/// from Impact
template< class O=int>
Teuchos::RCP<GridSizeGlobal<O,3> > createGridSizeGlobal() {

  Teuchos::Tuple<O,3> bla;
  SVS_get_nGlo(bla[0],bla[1],bla[2]);

  return(
      Teuchos::rcp(
          new GridSizeGlobal<O,3>( bla ) ) );
}


/// \relates GridSizeGlobal
/// sets impact
template< class O=int>
Teuchos::RCP<GridSizeGlobal<O,3> > createGridSizeGlobal( O n1, O n2, O n3 ) {
  return(
      Teuchos::rcp(
          new GridSizeGlobal<O,3>( n1, n2, n3 ) ) );
}

/// \relates GridSizeGlobal
/// sets impact
template< class O=int>
Teuchos::RCP<GridSizeGlobal<O,4> > createGridSizeGlobal( O n1, O n2, O n3, O nt ) {
  return(
      Teuchos::rcp(
          new GridSizeGlobal<O,4>( n1, n2, n3, nt ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
