#pragma once
#ifndef PIMPACT_MGSPACES_HPP
#define PIMPACT_MGSPACES_HPP

#include <vector>

#include "Pimpact_Space.hpp"
#include "Pimpact_CoarsenStrategy.hpp"




namespace Pimpact {


/// \brief contains hierarchy of Spaces
/// \ingroup MG
/// \tparam FSpaceT \c Space type on finest level
/// \tparam CSpaceT \c Space type on coarser levels allows to have tighter boundaries, less memory used, less comunication
template<class FST, class CST>
class MGSpaces {

public:

  typedef FST FSpaceT;
  typedef CST CSpaceT;

protected:

  template<class FSpaceTT, class CSpaceTT, class CoarsenStrategyT>
  friend Teuchos::RCP<const MGSpaces<FSpaceTT,CSpaceTT> > createMGSpaces( const Teuchos::RCP<const FSpaceTT>& space, int nGridsMax=10 );

  Teuchos::RCP<const FSpaceT> space_;
  std::vector< Teuchos::RCP<const CSpaceT> > spaces_;

  MGSpaces(
      const Teuchos::RCP<const FSpaceT>& space,
      const std::vector<Teuchos::RCP<const CSpaceT> >& spaces ):
    space_(space),
    spaces_(spaces) {}

public:

  int getNGrids() const { return( spaces_.size() ); }
  Teuchos::RCP<const FSpaceT>  get()        const { return( space_ ); }

  /// \brief gets ith space, similar to python i=-1 is gets you the coarses space
  Teuchos::RCP<const CSpaceT>  get( int i ) const {
    if( i<0 )
      return( spaces_[ getNGrids()+i ] );
    else
      return( spaces_[i] );
  }

  void print(  std::ostream& out=std::cout ) const {

    for( int i = 0; i<spaces_.size(); ++i ) {
      out << "-------- space: "<< i << "--------\n";
      spaces_[i]->print(out);
    }
  }

}; // end of class MGSpaces



/// \relates MGSpaces
template<class FSpaceT, class CSpaceT, class CoarsenStrategy>
Teuchos::RCP<const MGSpaces<FSpaceT,CSpaceT> >
createMGSpaces(
    const Teuchos::RCP<const FSpaceT>& space,
    int maxGrids=10 ) {

  std::vector<Teuchos::RCP<const CSpaceT> > spaces =
      CoarsenStrategy::getMultiSpace( space, maxGrids );

  return(
      Teuchos::rcp(
          new MGSpaces<FSpaceT,CSpaceT>(
              space,
              spaces ) ) );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MGSPACES_HPP
