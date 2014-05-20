#pragma once
#ifndef PIMPACT_SPACE_HPP
#define PIMPACT_SPACE_HPP

//#include "mpi.h"

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"

#include <iostream>



namespace Pimpact {


/// \todo integrate this
template<class Ordinal>
class Space {

  typedef Teuchos::Tuple< Teuchos::RCP<const IndexSpace<Ordinal> >, 3 >  IndexSpaces;

public:

  Space(
      const Teuchos::RCP<const FieldSpace<Ordinal> > fieldSpace=Teuchos::null,
      const IndexSpace& innerIS=Teuchos::null,
      const IndexSpace& fullIS=Teuchos::null ):
        fieldSpace_(fieldSpace),
        innerIS_(innerIS),
        fullIS_(fullIS) {}

  Teuchos::RCP<const FieldSpace<Ordinal> > fieldSpace_;

  IndexSpaces innerIS_;
  IndexSpaces fullIS_;


  const MPI_Fint& commf() const { return( fieldSpace_->commf_ ); }
  const MPI_Comm& comm()  const { return( fieldSpace_->comm_  ); }
  const int&      dim()   const { return( fieldSpace_->dim_   ); }

  const Ordinal*  nGlo()  const { return( fieldSpace_->nGlo_.getRawPtr() ); }
  const Ordinal*  nLoc()  const { return( fieldSpace_->nLoc_.getRawPtr() ); }
  const Ordinal*  bl  ()  const { return( fieldSpace_->bl_.getRawPtr()   ); }
  const Ordinal*  bu  ()  const { return( fieldSpace_->bu_.getRawPtr()   ); }


  const Ordinal* sInd() const { return( fieldSpace_->sInd_.getRawPtr()  ); }
  const Ordinal* eInd() const { return( fieldSpace_->eInd_.getRawPtr()  ); }

  const Ordinal* sInd(  int fieldType ) const { return( innerIS_[fieldType]->sInd_.getRawPtr() ); }
  const Ordinal* eInd(  int fieldType ) const { return( innerIS_[fieldType]->eInd_.getRawPtr() ); }

  const Ordinal* sIndB( int fieldType ) const { return( fullIS_[fieldType]->sInd_.getRawPtr()  ); }
  const Ordinal* eIndB( int fieldType ) const { return( fullIS_[fieldType]->eInd_.getRawPtr()  ); }


}; // end of class Space


/// \relates Space
template<class Ordinal>
Teuchos::RCP< const Space<Ordinal> > createSpace(
    const Teuchos::RCP< const FieldSpace<Ordinal> > fieldSpace=Teuchos::null,
    const IndexSpace& innerIS=Teuchos::null,
    const IndexSpace& fullIS=Teuchos::null ) {

  return(
      Teuchos::rcp(
          new Space( fieldSpace, innerIS, fullIS ) ) );
}




} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SPACE_HPP
