#pragma once
#ifndef PIMPACT_SPACE_HPP
#define PIMPACT_SPACE_HPP


#include "Teuchos_Tuple.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"

#include <iostream>



namespace Pimpact {


/// \todo integrate this
template<class Ordinal>
class Space {

public:

  typedef Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > >  IndexSpaces;

  Space(
      const Teuchos::RCP<const FieldSpace<Ordinal> > fieldSpace=Teuchos::null,
      const IndexSpaces& innerIS=Teuchos::null,
      const IndexSpaces& fullIS=Teuchos::null ):
        fieldSpace_(fieldSpace),
        innerIS_(innerIS),
        fullIS_(fullIS) {}

protected:

  Teuchos::RCP<const FieldSpace<Ordinal> > fieldSpace_;

  IndexSpaces innerIS_;
  IndexSpaces fullIS_;

public:

  Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const { return( fieldSpace_ ); }
  IndexSpaces getInnerIndexSpace() const { return( innerIS_ ); }
  IndexSpaces getFullIndexSpace() const { return( fullIS_ ); }


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

  void print() const {
    std::cout << "\t---Space: ---\n";
        fieldSpace_->print();
    for( int i=0; i<dim(); ++i ) {
      innerIS_[i]->print();
      fullIS_[i]->print();
    }
  }

}; // end of class Space


/// \relates Space
template<class Ordinal>
Teuchos::RCP< const Space<Ordinal> > createSpace(
    const Teuchos::RCP< const FieldSpace<Ordinal> > fieldSpace=Teuchos::null,
    const typename Space<Ordinal>::IndexSpaces& innerIS=Teuchos::null,
    const typename Space<Ordinal>::IndexSpaces& fullIS=Teuchos::null ) {

  if( fieldSpace.is_null() &&  innerIS.is_null() && fullIS.is_null() )
    return(
        Teuchos::rcp(
            new Space<Ordinal>(
                createFieldSpace<Ordinal>(),
                createInnerFieldIndexSpaces<Ordinal>(),
                createFullFieldIndexSpaces<Ordinal>() ) ) );
  else
    return(
        Teuchos::rcp(
            new Space<Ordinal>( fieldSpace, innerIS, fullIS ) ) );
}




} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SPACE_HPP
