#pragma once
#ifndef PIMPACT_SPACE_HPP
#define PIMPACT_SPACE_HPP


#include "Teuchos_Tuple.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_IndexSpace.hpp"

#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_ProcGrid.hpp"

#include <iostream>



namespace Pimpact {


/// \todo integrate this
template<class Ordinal=int, int dimension=3 >
class Space {

public:

  typedef Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > >  IndexSpaces;

  Space(
      const Teuchos::RCP<const FieldSpace<Ordinal,dimension> >& fieldSpace,
      const Teuchos::RCP<const IndexSpace<Ordinal> >& scalarIS,
      const IndexSpaces& innerIS,
      const IndexSpaces& fullIS,
      const Teuchos::RCP< GridSizeGlobal<Ordinal,dimension> >& gridSizeGlobal,
      const Teuchos::RCP< GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
      const Teuchos::RCP< ProcGridSize<Ordinal,dimension> >& procGridSize,
      const Teuchos::RCP< ProcGrid<Ordinal,dimension> >& procGrid ):
        fieldSpace_(fieldSpace),
        scalarIS_(scalarIS),
        innerIS_(innerIS),
        fullIS_(fullIS),
        gridSizeGlobal_(gridSizeGlobal),
        gridSizeLocal_(gridSizeLocal),
        procGridSize_(procGridSize),
        procGrid_(procGrid)
  {}

protected:

  Teuchos::RCP<const FieldSpace<Ordinal,dimension> > fieldSpace_;

  Teuchos::RCP<const IndexSpace<Ordinal> > scalarIS_;

  IndexSpaces innerIS_;
  IndexSpaces fullIS_;

  Teuchos::RCP< GridSizeGlobal<Ordinal,dimension> > gridSizeGlobal_;
  Teuchos::RCP< GridSizeLocal<Ordinal,dimension> > gridSizeLocal_;

  Teuchos::RCP< ProcGridSize<Ordinal,dimension> > procGridSize_;

  Teuchos::RCP< ProcGrid<Ordinal,dimension> > procGrid_;

public:

  //  Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const { return( fieldSpace_ ); }
  //  IndexSpaces getInnerIndexSpace() const { return( innerIS_ ); }
  //  IndexSpaces getFullIndexSpace() const { return( fullIS_ ); }

  const MPI_Fint& commf() const { return( fieldSpace_->commf_ ); }
  const MPI_Comm& comm()  const { return( fieldSpace_->comm_  ); }

//  const MPI_Fint& commf() const { return( fieldSpace_->commf_ ); }
  const MPI_Comm& commST()  const { return( procGrid_->commSpaceTime_  ); }

  int rankST() const { return( procGrid_->rankST_ ); }

  const int&      dim()   const { return( fieldSpace_->dim_   ); }

  const Ordinal* nGlo()  const { return( gridSizeGlobal_->getPtr() ); }
  const Ordinal* nLoc()  const { return( gridSizeLocal_->getPtr() ); }
  const Ordinal* bl  ()  const { return( fieldSpace_->bl_.getRawPtr()   ); }
  const Ordinal* bu  ()  const { return( fieldSpace_->bu_.getRawPtr()   ); }


  const Ordinal* sInd() const { return( scalarIS_->sInd_.getRawPtr()  ); }
  const Ordinal* eInd() const { return( scalarIS_->eInd_.getRawPtr()  ); }

  const Ordinal* sInd(  int fieldType ) const { return( innerIS_[fieldType]->sInd_.getRawPtr() ); }
  const Ordinal* eInd(  int fieldType ) const { return( innerIS_[fieldType]->eInd_.getRawPtr() ); }

  const Ordinal* sIndB( int fieldType ) const { return( fullIS_[fieldType]->sInd_.getRawPtr()  ); }
  const Ordinal* eIndB( int fieldType ) const { return( fullIS_[fieldType]->eInd_.getRawPtr()  ); }

  const Ordinal* rankU() const { return( procGrid_->rankU_.getRawPtr()  ); }
  const Ordinal* rankL() const { return( procGrid_->rankL_.getRawPtr()  ); }

  const Ordinal* procCoordinate() const { return( procGrid_->iB_.getRawPtr()  ); }

  const Ordinal* shift() const { return( procGrid_->shift_.getRawPtr()  ); }

  const Ordinal& getNProc(int i) const { return( procGridSize_->get(i) ); }


  void print(  std::ostream& out=std::cout ) const {

    out << "\t---Space: ---\n";

    if( !fieldSpace_.is_null() )
      fieldSpace_->print( out );
    else
      out << "fieldSpace_ is null\n";
    MPI_Barrier( commST() );

    if( !gridSizeGlobal_.is_null() ) {
      out <<"---GridSizeGlobal: ---\n";
      gridSizeGlobal_->print( out );
    }
    else
      out << "gridSizeGlobal_ is null\n";
    MPI_Barrier( commST() );

    if( !gridSizeLocal_.is_null() ) {
      out <<"---GridSizeLocal: ---\n";
      gridSizeLocal_->print( out );
    }
    else
      out << "gridSizeLocal_ is null\n";
    MPI_Barrier( commST() );

    if( !scalarIS_.is_null() )
      scalarIS_->print(out);
    else
      out << "scalarIS_ is null\n";
    MPI_Barrier( commST() );

    for( int i=0; i<2; ++i ) {
      innerIS_[i]->print( out );
      fullIS_[i]->print( out );
    }
    MPI_Barrier( comm() );

    if( !procGridSize_.is_null() )
      procGridSize_->print( out );
    else
      out << "procGridSize_ is null\n";
    MPI_Barrier( commST() );

    if( !procGrid_.is_null() )
      procGrid_->print( out );
    else
      out << "procGrid_ is null\n";
    MPI_Barrier( commST() );

  }

}; // end of class Space


/// \relates Space
/// \todo wünschenswert initialization from parameterlist
template< class O=int, int d=3>
Teuchos::RCP< const Space<O,d> > createSpace(
    const Teuchos::RCP<const FieldSpace<O,d> >& fieldSpace,
    const Teuchos::RCP<const IndexSpace<O> >& scalarIS,
    const Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > >& innerIS,
    const Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > >& fullIS,
    const Teuchos::RCP< GridSizeGlobal<O,d> >& gridSizeGlobal,
    const Teuchos::RCP< GridSizeLocal<O,d> >& gridSizeLocal,
    const Teuchos::RCP< ProcGridSize<O,d> >& procGridSize,
    const Teuchos::RCP< ProcGrid<O,d> >& procGrid ) {

  return(
      Teuchos::rcp(
          new Space<O,d>(
              fieldSpace,
              scalarIS,
              innerIS,
              fullIS,
              gridSizeGlobal,
              gridSizeLocal,
              procGridSize,
              procGrid ) ) );
}


/// \relates Space
template< class O=int, int d=3>
Teuchos::RCP< const Space<O,d> > createSpace() {

  return(
      Teuchos::rcp(
          new Space<O,d>(
              createFieldSpace<O>(),
              createScalarIndexSpace<O>(),
              createInnerFieldIndexSpaces<O>(),
              createFullFieldIndexSpaces<O>(),
              createGridSizeGlobal<O>(),
              createGridSizeLocal<O,3>(),
              Teuchos::null,
              Teuchos::null ) ) );
}
} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SPACE_HPP
