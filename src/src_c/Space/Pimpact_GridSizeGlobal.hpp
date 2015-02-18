#pragma once
#ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
#define PIMPACT_GRIDSIZEGLOBAL_HPP

#include<ostream>

#include"Teuchos_RCP.hpp"
#include"Teuchos_Tuple.hpp"

#include"Pimpact_Types.hpp"




extern "C" {

void fsetGS(const int& n1, const int& n2, const int& n3 );

}



namespace Pimpact{


/// \brief global grid size(independent of FieldType)
///
///
/// \todo include setter method, such that it gets update with enlarging Fourier modes!!!
/// one could think about inheriting from Tuple, or generalize for global and local use
/// \ingroup Space
template<class Ordinal>
class GridSizeGlobal {

  template<class OT>
  friend Teuchos::RCP<const GridSizeGlobal<OT> > createGridSizeGlobal( OT n1, OT n2, OT n3, OT nt );

  template<class OT>
  friend Teuchos::RCP<const GridSizeGlobal<OT> > createGridSizeGlobal( const Teuchos::Tuple<OT,4>& tuple );

public:

  typedef const Teuchos::Tuple<Ordinal,4> TO;

protected:

  TO gridSize_;

  GridSizeGlobal( TO gridSize ):
    gridSize_( gridSize ) {};

public:

  void set_Impact() const {
    fsetGS( gridSize_[0], gridSize_[1], gridSize_[2] );
  };

  const Ordinal& get( int i ) const  {
    return( gridSize_[i] );
  }

  const Ordinal* get() const {
    return( gridSize_.getRawPtr() );
  }

  TO getTuple() const {
    return( gridSize_ );
  }

  void print( std::ostream& out=std::cout ) const {
    out << " --- GridSizeGlobal: " << "\n";
    out << " \tnx=" << gridSize_[0];
    out << " \tny=" << gridSize_[1];
    out << " \tnz=" << gridSize_[2];
    out << " \tnt=" << gridSize_[3];
    out << "\n";

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
