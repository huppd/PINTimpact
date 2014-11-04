#pragma once
#ifndef PIMPACT_INDEXSPACE_HPP
#define PIMPACT_INDEXSPACE_HPP

#include "mpi.h"

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include <iostream>

#include "Pimpact_Types.hpp"

#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"



namespace Pimpact {


extern "C" {

void SVS_set_sInd(const int* const);
void SVS_set_eInd(const int* const);

void VS_set_sIndU(const int* const);
void VS_set_eIndU(const int* const);

void VS_set_sIndUB(const int* const);
void VS_set_eIndUB(const int* const);

void VS_set_sIndV(const int* const);
void VS_set_eIndV(const int* const);

void VS_set_sIndVB(const int* const);
void VS_set_eIndVB(const int* const);

void VS_set_sIndW(const int* const);
void VS_set_eIndW(const int* const);

void VS_set_sIndWB(const int* const);
void VS_set_eIndWB(const int* const);

}


/// \brief class that stores neccessary lower and upper indexes for \c ScalarField for different FieldTypes
///
/// Index         | BoundaryConditionslocal>0 | BoundaryConditionslocal<=0
/// ------------- | --------------------------| -------------
/// sInd          | 1                         | 1
/// eInd          | nLoc                      | nLoc-1
///
/// Index           | BoundaryConditionslocal>0 | BoundaryConditionslocal<=0 | symmetryBoundaryConditions
/// --------------- | --------------------------| ---------------------------| ---------------------------
/// sInd[dir==field]| 1                         | 2                          | 1
/// sInd[dir!=field]| 2                         | 2                          | 1
/// eInd[dir==field]| nLoc-1                    | nLoc-1                     | nLoc-1
/// eInd[dir!=field]| nLoc-1                    | nLoc-1                     | nLoc
///
/// Index            | BoundaryConditionslocal>0 | BoundaryConditionslocal<=0 | symmetryBoundaryConditions
/// ---------------- | --------------------------| ---------------------------| ---------------------------
/// sIndB[dir==field]| 0                         | 1                          | 1
/// sIndB[dir!=field]| 1                         | 1                          | 1
/// eIndB[dir==field]| nLoc                      | nLoc-1                     | nLoc-1
/// eIndB[dir!=field]| nLoc                      | nLoc-1                     | nLoc
///
/// and including Boundaries or excluding them.
/// \todo remove fieldType_
/// \todo move handling from space here
/// \todo make constructor private
/// \ingroup Space
template< class Ordinal = int >
class IndexSpace {

public:

  typedef Teuchos::Tuple<Ordinal,3> TO3;

  /// \brief constructor
  /// \param fieldType says which kind of field is taken
  /// \param sInd start index of gridpoints
  /// \param eInd last index of gridpoints
  IndexSpace( EField fieldType,TO3 sInd, TO3 eInd ):
    fieldType_(fieldType),sInd_(sInd),eInd_(eInd) {};


  void print( std::ostream& out=std::cout ) const {
    out << "\t---IndexSpace: ---\n";
    out << "fieldType: " << fieldType_ << "\n";
    out << "sInd: " << sInd_ << "\n";
    out << "eInd: " << eInd_ << "\n";
  }

  EField fieldType_;

  TO3 sInd_;
  TO3 eInd_;

}; // end of class IndexSpace



/// \brief function that creates IndexSpace for \c EFieldType::S
///
/// \relates IndexSpace
/// \deprecated  \param setImpact
template<class O=int, int d=3>
Teuchos::RCP<const IndexSpace<O> >
createScalarIndexSpace(
    const Teuchos::RCP<const StencilWidths<d> >& sW,
    const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
    const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
    bool setImpact=true ) {

  typedef typename IndexSpace<O>::TO3 TO3;

  TO3 sInd;
  TO3 eInd;

  for( int i=0; i<3; ++i ) {
    sInd[i] = 2            + sW->getLS(i);
    eInd[i] = nLoc->get(i) + sW->getLS(i);
    if( bc->BCL_int_[i] > 0 )
      sInd[i] = 1;
    if( bc->BCU_int_[i] > 0 )
      eInd[i] = nLoc->get(i);
  }

  if( setImpact ) {
    SVS_set_sInd( sInd.getRawPtr() );
    SVS_set_eInd( eInd.getRawPtr() );
  }

  return( Teuchos::rcp(
      new IndexSpace<O>( EField::S, sInd, eInd ) ) );

}


/// \brief function that creates IndexSpace for inner Velocity Field
///
/// \todo make construction clean without copying
/// \relates IndexSpace
template<class O=int, int d=3>
Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<O> > >
createInnerFieldIndexSpaces(
    const Teuchos::RCP<const StencilWidths<d> >& sW,
    const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
    const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
    bool setImpact=true ){

  typedef typename IndexSpace<O>::TO3 TO3;

  Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > > fIS(3);

  TO3 sIndU;
  TO3 eIndU;

  TO3 sIndV;
  TO3 eIndV;

  TO3 sIndW;
  TO3 eIndW;

  for( int i=0; i<3; ++i ) {
    sIndU[i] = 2            + sW->getLS(i);
    eIndU[i] = nLoc->get(i) + sW->getLS(i);

    sIndV[i] = 2            + sW->getLS(i);
    eIndV[i] = nLoc->get(i) + sW->getLS(i);

    sIndW[i] = 2            + sW->getLS(i);
    eIndW[i] = nLoc->get(i) + sW->getLS(i);
  }

  // Lower index in x-direction
  int i = 0;
  if( bc->BCL_int_[i] > 0 ) {
    sIndU[i] = 1;
    sIndV[i] = 2;
    sIndW[i] = 2;
  }
  // lower index in y-direction
  i = 1;
  if( bc->BCL_int_[i] > 0 ) {
    sIndU[i] = 2;
    sIndV[i] = 1;
    sIndW[i] = 2;
  }
  // lower index in z-direction
  i = 2;
  if( bc->BCL_int_[i] > 0 ) {
    sIndU[i] = 2;
    sIndV[i] = 2;
    sIndW[i] = 1;
  }

  // lower index for symmetriBC
  for( int i=0; i<3; ++i )
    if( bc->BCL_local_[i]==SymmetryBC ) {
      sIndU[i] = 1;
      sIndV[i] = 1;
      sIndW[i] = 1;
    }

  // upper index
  for( int i=0; i<3; ++i )
    if( bc->BCU_int_[i] > 0 ) {
      eIndU[i] = nLoc->get(i)-1;
      eIndV[i] = nLoc->get(i)-1;
      eIndW[i] = nLoc->get(i)-1;
    }

  // upper index in x-direction for symmetricBC
  i=0;
  if( bc->BCU_local_[i]==SymmetryBC ) {
    eIndU[i] = nLoc->get(i)-1;
    eIndV[i] = nLoc->get(i);
    eIndW[i] = nLoc->get(i);
  }

  // upper index in y-direction for symmetricBC
  i=1;
  if( bc->BCU_local_[i]==SymmetryBC ) {
    eIndU[i] = nLoc->get(i);
    eIndV[i] = nLoc->get(i)-1;
    eIndW[i] = nLoc->get(i);
  }

  // upper index in z-direction for symmetricBC
  i=2;
  if( bc->BCU_local_[i]==SymmetryBC ) {
    eIndU[i] = nLoc->get(i);
    eIndV[i] = nLoc->get(i);
    eIndW[i] = nLoc->get(i)-1;
  }


  if( setImpact ) {
    VS_set_sIndU( sIndU.getRawPtr() );
    VS_set_eIndU( eIndU.getRawPtr() );
  }
  fIS[0] =  Teuchos::rcp( new IndexSpace<O>( EField::U, sIndU, eIndU ) );

  if( setImpact ) {
    VS_set_sIndV( sIndV.getRawPtr() );
    VS_set_eIndV( eIndV.getRawPtr() );
  }
  fIS[1] =  Teuchos::rcp( new IndexSpace<O>( EField::V, sIndV, eIndV ) );

  if( setImpact ) {
    VS_set_sIndW( sIndW.getRawPtr() );
    VS_set_eIndW( eIndW.getRawPtr() );
  }
  fIS[2] =  Teuchos::rcp( new IndexSpace<O>( EField::W, sIndW, eIndW ) );

  return( fIS );
}



/// \brief function that creates ScaparIndexSpace
///
/// by getting values from impact
/// \todo make construction clean without copying
/// \relates IndexSpace
template<class O=int, int d=3>
Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<O> > >
createFullFieldIndexSpaces(
    const Teuchos::RCP<const StencilWidths<d> >& sW,
    const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
    const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
    bool setImpact=true ){


  typedef typename IndexSpace<O>::TO3 TO3;

  Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > > fIS(3);

  TO3 sIndU;
  TO3 eIndU;

  TO3 sIndV;
  TO3 eIndV;

  TO3 sIndW;
  TO3 eIndW;

  for( int i=0; i<3; ++i ) {
    sIndU[i] = 2            + sW->getLS(i);
    eIndU[i] = nLoc->get(i) + sW->getLS(i);

    sIndV[i] = 2            + sW->getLS(i);
    eIndV[i] = nLoc->get(i) + sW->getLS(i);

    sIndW[i] = 2            + sW->getLS(i);
    eIndW[i] = nLoc->get(i) + sW->getLS(i);
  }

  // Lower index in x-direction
  if( bc->BCL_int_[0] > 0 ) {
    sIndU[0] = 0;
    sIndV[0] = 1;
    sIndW[0] = 1;
  }
  // lower index in y-direction
  if( bc->BCL_int_[1] > 0 ) {
    sIndU[1] = 1;
    sIndV[1] = 0;
    sIndW[1] = 1;
  }
  // lower index in z-direction
  if( bc->BCL_int_[2] > 0 ) {
    sIndU[2] = 1;
    sIndV[2] = 1;
    sIndW[2] = 0;
  }

  // lower index for symmetriBC
  for( int i=0; i<3; ++i )
    if( bc->BCL_local_[i]==SymmetryBC ) {
      sIndU[i] = 1;
      sIndV[i] = 1;
      sIndW[i] = 1;
    }

  // upper index
  for( int i=0; i<3; ++i )
    if( bc->BCU_int_[i] > 0 ) {
      eIndU[i] = nLoc->get(i);
      eIndV[i] = nLoc->get(i);
      eIndW[i] = nLoc->get(i);
    }

  // upper index in x-direction for symmetricBC
  if( bc->BCU_local_[0]==SymmetryBC ) {
    eIndU[0] = nLoc->get(0)-1;
    eIndV[0] = nLoc->get(0);
    eIndW[0] = nLoc->get(0);
  }

  // upper index in y-direction for symmetricBC
  if( bc->BCU_local_[1]==SymmetryBC ) {
    eIndU[1] = nLoc->get(1);
    eIndV[1] = nLoc->get(1)-1;
    eIndW[1] = nLoc->get(1);
  }

  // upper index in z-direction for symmetricBC
  if( bc->BCU_local_[2]==SymmetryBC ) {
    eIndU[2] = nLoc->get(2);
    eIndV[2] = nLoc->get(2);
    eIndW[2] = nLoc->get(2)-1;
  }


  if( setImpact ) {
    VS_set_sIndUB( sIndU.getRawPtr() );
    VS_set_eIndUB( eIndU.getRawPtr() );
  }
  fIS[0] =  Teuchos::rcp( new IndexSpace<O>( EField::U, sIndU, eIndU ) );

  if( setImpact ) {
    VS_set_sIndVB( sIndV.getRawPtr() );
    VS_set_eIndVB( eIndV.getRawPtr() );
  }
  fIS[1] =  Teuchos::rcp( new IndexSpace<O>( EField::V, sIndV, eIndV ) );

  if( setImpact ) {
    VS_set_sIndWB( sIndW.getRawPtr() );
    VS_set_eIndWB( eIndW.getRawPtr() );
  }
  fIS[2] =  Teuchos::rcp( new IndexSpace<O>( EField::W, sIndW, eIndW ) );

  return( fIS );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INDEXSPACE_HPP
