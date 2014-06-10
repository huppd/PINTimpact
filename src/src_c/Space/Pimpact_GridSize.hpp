#pragma once
#ifndef PIMPACT_GRIDSIZE_HPP
#define PIMPACT_GRIDSIZE_HPP

#include<ostream>

#include"Teuchos_RCP.hpp"
#include"Teuchos_Tuple.hpp"

#include"Pimpact_Types.hpp"

extern "C" {
void fsetGS(const int& n1, const int& n2, const int& n3 );
}

namespace Pimpact{


/// \todo merge with Space should make nGlow in FieldSpace obsolete, maybe rename
template<class Ordinal>
class GridSize {

public:

  typedef const Teuchos::Tuple<Ordinal,3> TO3;

protected:

  TO3 gridSize_;

public:

  GridSize( Ordinal n1, Ordinal n2, Ordinal n3 ):
    gridSize_( Teuchos::tuple( n1, n2, n3) ) {
    set_Impact();
  };

  GridSize( TO3 domainSize ):
    gridSize_( domainSize ) {};

  void set_Impact(){
    fsetGS( gridSize_[0], gridSize_[1], gridSize_[2] );
  };

  void print( std::ostream& out ) {
    out  << " \tnx=" << gridSize_[0]
         << " \tny=" << gridSize_[1]
         << " \tnz=" << gridSize_[2] << "\n";
  };

}; // end of class GridSize


/// \relates GridSize
template<class Ordinal>
Teuchos::RCP<GridSize<Ordinal> > createGridSize( Ordinal n1=1, Ordinal n2=1, Ordinal n3=1 ) {
  return(
      Teuchos::rcp(
          new GridSize<Ordinal>( n1, n2, n3 ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZE_HPP
