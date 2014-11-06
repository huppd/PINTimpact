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
/// and including Boundaries or excluding them.
///
/// Index         | BoundaryConditionslocal>0 | BoundaryConditionslocal<=0 | symmetryBoundaryConditions
/// ------------- | --------------------------| ---------------------------| ---------------------------
/// sInd          | 1                         | 1                          | 1
/// eInd          | nLoc                      | nLoc-1                     | nLoc
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
/// \ingroup Space
template< class Ordinal = int >
class IndexSpace {

  template<class O, int d>
  friend Teuchos::RCP<const IndexSpace<O> >
  createIndexSpace(
      const Teuchos::RCP<const StencilWidths<d> >& sW,
      const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
      const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
      bool setImpact=true );

public:

  typedef Teuchos::Tuple<Ordinal,3> TO3;

  typedef Teuchos::Tuple< Teuchos::Tuple<Ordinal, 3>, 3 > TTO3;

protected:
  /// \brief constructor
  /// \param fieldType says which kind of field is taken
  /// \param sInd start index of gridpoints
  /// \param eInd last index of gridpoints
  template<int d=3>
  IndexSpace(
      const Teuchos::RCP<const StencilWidths<d> >& sW,
      const Teuchos::RCP<const GridSizeLocal<Ordinal,d> >& nLoc,
      const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
      bool setImpact=true ) {

    // ------------------init IndS_--------------------
    for( int i=0; i<3; ++i ) {
      sIndS_[i] = 2            + sW->getLS(i);
      eIndS_[i] = nLoc->get(i) + sW->getLS(i);
      if( bc->BCL_int_[i] > 0 )
        sIndS_[i] = 1;
      if( bc->BCU_int_[i] > 0 )
        eIndS_[i] = nLoc->get(i);
      if( bc->BCL_local_[i]==SymmetryBC )
        sIndS_[i] = 1;
      if( bc->BCU_local_[i]==SymmetryBC )
        eIndS_[i] = nLoc->get(i);
    }

    if( setImpact ) {
      SVS_set_sInd( sIndS_.getRawPtr() );
      SVS_set_eInd( eIndS_.getRawPtr() );
    }
    // --- init IndU_ -------------------------------

    for( int field=0; field<3; ++field )
      for( int dir=0; dir<3; ++dir ) {
        sIndU_[field][dir] = 2              + sW->getLS(dir);
        eIndU_[field][dir] = nLoc->get(dir) + sW->getLS(dir);
    }

    // Lower index in x-direction
    if( bc->BCL_int_[0] > 0 ) {
      sIndU_[U][0] = 1;
      sIndU_[V][0] = 2;
      sIndU_[W][0] = 2;
    }
    // lower index in y-direction
    if( bc->BCL_int_[1] > 0 ) {
      sIndU_[U][1] = 2;
      sIndU_[V][1] = 1;
      sIndU_[W][1] = 2;
    }
    // lower index in z-direction
    if( bc->BCL_int_[2] > 0 ) {
      sIndU_[U][2] = 2;
      sIndU_[V][2] = 2;
      sIndU_[W][2] = 1;
    }

    // lower index for symmetriBC
    for( int i=0; i<3; ++i )
      if( bc->BCL_local_[i]==SymmetryBC ) {
        sIndU_[U][i] = 1;
        sIndU_[V][i] = 1;
        sIndU_[W][i] = 1;
      }

    // upper index
    for( int i=0; i<3; ++i )
      if( bc->BCU_int_[i] > 0 ) {
        eIndU_[U][i] = nLoc->get(i)-1;
        eIndU_[V][i] = nLoc->get(i)-1;
        eIndU_[W][i] = nLoc->get(i)-1;
      }

    // upper index in x-direction for symmetricBC
    if( bc->BCU_local_[0]==SymmetryBC ) {
      eIndU_[U][0] = nLoc->get(0)-1;
      eIndU_[V][0] = nLoc->get(0);
      eIndU_[W][0] = nLoc->get(0);
    }

    // upper index in y-direction for symmetricBC
    if( bc->BCU_local_[1]==SymmetryBC ) {
      eIndU_[U][1] = nLoc->get(1);
      eIndU_[V][1] = nLoc->get(1)-1;
      eIndU_[W][1] = nLoc->get(1);
    }

    // upper index in z-direction for symmetricBC
    if( bc->BCU_local_[2]==SymmetryBC ) {
      eIndU_[U][2] = nLoc->get(2);
      eIndU_[V][2] = nLoc->get(2);
      eIndU_[W][2] = nLoc->get(2)-1;
    }

    if( setImpact ) {
      VS_set_sIndU( sIndU_[U].getRawPtr() );
      VS_set_eIndU( eIndU_[U].getRawPtr() );
      VS_set_sIndV( sIndU_[V].getRawPtr() );
      VS_set_eIndV( eIndU_[V].getRawPtr() );
      VS_set_sIndW( sIndU_[W].getRawPtr() );
      VS_set_eIndW( eIndU_[W].getRawPtr() );
    }

    // --- init IndUB_ -------------------------------
    for( int i=0; i<3; ++i ) {
      sIndUB_[U][i] = 2            + sW->getLS(i);
      eIndUB_[U][i] = nLoc->get(i) + sW->getLS(i);

      sIndUB_[V][i] = 2            + sW->getLS(i);
      eIndUB_[V][i] = nLoc->get(i) + sW->getLS(i);

      sIndUB_[W][i] = 2            + sW->getLS(i);
      eIndUB_[W][i] = nLoc->get(i) + sW->getLS(i);
    }

    // Lower index in x-direction
    if( bc->BCL_int_[0] > 0 ) {
      sIndUB_[U][0] = 0;
      sIndUB_[V][0] = 1;
      sIndUB_[W][0] = 1;
    }
    // lower index in y-direction
    if( bc->BCL_int_[1] > 0 ) {
      sIndUB_[U][1] = 1;
      sIndUB_[V][1] = 0;
      sIndUB_[W][1] = 1;
    }
    // lower index in z-direction
    if( bc->BCL_int_[2] > 0 ) {
      sIndUB_[U][2] = 1;
      sIndUB_[V][2] = 1;
      sIndUB_[W][2] = 0;
    }

    // lower index for symmetriBC
    for( int i=0; i<3; ++i )
      if( bc->BCL_local_[i]==SymmetryBC ) {
        sIndUB_[U][i] = 1;
        sIndUB_[V][i] = 1;
        sIndUB_[W][i] = 1;
      }

    // upper index
    for( int i=0; i<3; ++i )
      if( bc->BCU_int_[i] > 0 ) {
        eIndUB_[U][i] = nLoc->get(i);
        eIndUB_[V][i] = nLoc->get(i);
        eIndUB_[W][i] = nLoc->get(i);
      }

    // upper index in x-direction for symmetricBC
    if( bc->BCU_local_[0]==SymmetryBC ) {
      eIndUB_[U][0] = nLoc->get(0)-1;
      eIndUB_[V][0] = nLoc->get(0);
      eIndUB_[W][0] = nLoc->get(0);
    }

    // upper index in y-direction for symmetricBC
    if( bc->BCU_local_[1]==SymmetryBC ) {
      eIndUB_[U][1] = nLoc->get(1);
      eIndUB_[V][1] = nLoc->get(1)-1;
      eIndUB_[W][1] = nLoc->get(1);
    }

    // upper index in z-direction for symmetricBC
    if( bc->BCU_local_[2]==SymmetryBC ) {
      eIndUB_[U][2] = nLoc->get(2);
      eIndUB_[V][2] = nLoc->get(2);
      eIndUB_[W][2] = nLoc->get(2)-1;
    }


    if( setImpact ) {
      VS_set_sIndUB( sIndUB_[U].getRawPtr() );
      VS_set_eIndUB( eIndUB_[U].getRawPtr() );
      VS_set_sIndVB( sIndUB_[V].getRawPtr() );
      VS_set_eIndVB( eIndUB_[V].getRawPtr() );
      VS_set_sIndWB( sIndUB_[W].getRawPtr() );
      VS_set_eIndWB( eIndUB_[W].getRawPtr() );
    }

  }

public:

  void print( std::ostream& out=std::cout ) const {
    out << "\t---IndexSpace: ---\n";
    out << "fieldType: S\n";
    out << "sInd: " << sIndS_ << "\n";
    out << "eInd: " << eIndS_ << "\n";
    for( int field=0; field<3; ++field ) {
      out << "inner field: " << field << "\n";
      out << "sInd: " << sIndU_[field] << "\n";
      out << "eInd: " << eIndU_[field] << "\n";
    }
    for( int field=0; field<3; ++field ) {
      out << "full field: " << field << "\n";
      out << "sInd: " << sIndUB_[field] << "\n";
      out << "eInd: " << eIndUB_[field] << "\n";
    }
  }

  const Ordinal* sInd( int fieldType ) const {
    if( EField::S == (EField)fieldType )
      return( sIndS_.getRawPtr()  );
    else
      return( sIndU_[fieldType].getRawPtr() );
  }
  const Ordinal* eInd(  int fieldType ) const {
    if( EField::S == (EField)fieldType )
      return( eIndS_.getRawPtr()  );
    else
      return( eIndU_[fieldType].getRawPtr() );
  }

  const Ordinal* sIndB( int fieldType ) const {
    if( EField::S == (EField)fieldType )
      return( sIndS_.getRawPtr()  );
    else
      return( sIndUB_[fieldType].getRawPtr()  );
  }
  const Ordinal* eIndB( int fieldType ) const {
    if( EField::S == (EField)fieldType )
      return( eIndS_.getRawPtr()  );
    else
      return( eIndUB_[fieldType].getRawPtr()  );
  }

  const Ordinal& sInd( int fieldType, int dir ) const {
    if( EField::S == (EField)fieldType )
      return( sIndS_[dir]  );
    else
      return( sIndU_[fieldType][dir] );
  }
  const Ordinal& eInd(  int fieldType, int dir ) const {
    if( EField::S == (EField)fieldType )
      return( eIndS_[dir]  );
    else
      return( eIndU_[fieldType][dir] );
  }

  const Ordinal& sIndB( int fieldType, int dir ) const {
    if( EField::S == (EField)fieldType )
      return( sIndS_[dir]  );
    else
      return( sIndUB_[fieldType].getRawPtr()  );
  }
  const Ordinal& eIndB( int fieldType, int dir ) const {
    if( EField::S == (EField)fieldType )
      return( eIndS_[dir]  );
    else
      return( eIndUB_[fieldType][dir]  );
  }

protected:

  TO3 sIndS_;
  TO3 eIndS_;

  TTO3 sIndU_;
  TTO3 eIndU_;

  TTO3 sIndUB_;
  TTO3 eIndUB_;

}; // end of class IndexSpace

template<class O=int, int d=3>
Teuchos::RCP<const IndexSpace<O> >
createIndexSpace(
    const Teuchos::RCP<const StencilWidths<d> >& sW,
    const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
    const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
    bool setImpact=true ) {
  return(
      Teuchos::rcp(
          new IndexSpace<O>( sW,nLoc,bc )
      )
  );

}

///// \brief function that creates IndexSpace for \c EFieldType::S
/////
///// \relates IndexSpace
///// \deprecated  \param setImpact
//template<class O=int, int d=3>
//Teuchos::RCP<const IndexSpace<O> >
//createScalarIndexSpace(
//    const Teuchos::RCP<const StencilWidths<d> >& sW,
//    const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
//    const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
//    bool setImpact=true ) {
//
//  typedef typename IndexSpace<O>::TO3 TO3;
//
//  TO3 sInd;
//  TO3 eInd;
//
//  for( int i=0; i<3; ++i ) {
//    sInd[i] = 2            + sW->getLS(i);
//    eInd[i] = nLoc->get(i) + sW->getLS(i);
//    if( bc->BCL_int_[i] > 0 )
//      sInd[i] = 1;
//    if( bc->BCU_int_[i] > 0 )
//      eInd[i] = nLoc->get(i);
//    if( bc->BCL_local_[i]==SymmetryBC )
//      sInd[i] = 1;
//    if( bc->BCU_local_[i]==SymmetryBC )
//      eInd[i] = nLoc->get(i);
//  }
//
//  if( setImpact ) {
//    SVS_set_sInd( sInd.getRawPtr() );
//    SVS_set_eInd( eInd.getRawPtr() );
//  }
//
//  return( Teuchos::rcp(
//      new IndexSpace<O>( EField::S, sInd, eInd ) ) );
//
//}
//
//
///// \brief function that creates IndexSpace for inner Velocity Field
/////
///// \todo make construction clean without copying
///// \relates IndexSpace
//template<class O=int, int d=3>
//Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<O> > >
//createInnerFieldIndexSpaces(
//    const Teuchos::RCP<const StencilWidths<d> >& sW,
//    const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
//    const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
//    bool setImpact=true ){
//
//  typedef typename IndexSpace<O>::TO3 TO3;
//
//  Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > > fIS(3);
//
//  TO3 sIndU;
//  TO3 eIndU;
//
//  TO3 sIndV;
//  TO3 eIndV;
//
//  TO3 sIndW;
//  TO3 eIndW;
//
//  for( int i=0; i<3; ++i ) {
//    sIndU[i] = 2            + sW->getLS(i);
//    eIndU[i] = nLoc->get(i) + sW->getLS(i);
//
//    sIndV[i] = 2            + sW->getLS(i);
//    eIndV[i] = nLoc->get(i) + sW->getLS(i);
//
//    sIndW[i] = 2            + sW->getLS(i);
//    eIndW[i] = nLoc->get(i) + sW->getLS(i);
//  }
//
//  // Lower index in x-direction
//  int i = 0;
//  if( bc->BCL_int_[i] > 0 ) {
//    sIndU[i] = 1;
//    sIndV[i] = 2;
//    sIndW[i] = 2;
//  }
//  // lower index in y-direction
//  i = 1;
//  if( bc->BCL_int_[i] > 0 ) {
//    sIndU[i] = 2;
//    sIndV[i] = 1;
//    sIndW[i] = 2;
//  }
//  // lower index in z-direction
//  i = 2;
//  if( bc->BCL_int_[i] > 0 ) {
//    sIndU[i] = 2;
//    sIndV[i] = 2;
//    sIndW[i] = 1;
//  }
//
//  // lower index for symmetriBC
//  for( int i=0; i<3; ++i )
//    if( bc->BCL_local_[i]==SymmetryBC ) {
//      sIndU[i] = 1;
//      sIndV[i] = 1;
//      sIndW[i] = 1;
//    }
//
//  // upper index
//  for( int i=0; i<3; ++i )
//    if( bc->BCU_int_[i] > 0 ) {
//      eIndU[i] = nLoc->get(i)-1;
//      eIndV[i] = nLoc->get(i)-1;
//      eIndW[i] = nLoc->get(i)-1;
//    }
//
//  // upper index in x-direction for symmetricBC
//  i=0;
//  if( bc->BCU_local_[i]==SymmetryBC ) {
//    eIndU[i] = nLoc->get(i)-1;
//    eIndV[i] = nLoc->get(i);
//    eIndW[i] = nLoc->get(i);
//  }
//
//  // upper index in y-direction for symmetricBC
//  i=1;
//  if( bc->BCU_local_[i]==SymmetryBC ) {
//    eIndU[i] = nLoc->get(i);
//    eIndV[i] = nLoc->get(i)-1;
//    eIndW[i] = nLoc->get(i);
//  }
//
//  // upper index in z-direction for symmetricBC
//  i=2;
//  if( bc->BCU_local_[i]==SymmetryBC ) {
//    eIndU[i] = nLoc->get(i);
//    eIndV[i] = nLoc->get(i);
//    eIndW[i] = nLoc->get(i)-1;
//  }
//
//
//  if( setImpact ) {
//    VS_set_sIndU( sIndU.getRawPtr() );
//    VS_set_eIndU( eIndU.getRawPtr() );
//  }
//  fIS[0] =  Teuchos::rcp( new IndexSpace<O>( EField::U, sIndU, eIndU ) );
//
//  if( setImpact ) {
//    VS_set_sIndV( sIndV.getRawPtr() );
//    VS_set_eIndV( eIndV.getRawPtr() );
//  }
//  fIS[1] =  Teuchos::rcp( new IndexSpace<O>( EField::V, sIndV, eIndV ) );
//
//  if( setImpact ) {
//    VS_set_sIndW( sIndW.getRawPtr() );
//    VS_set_eIndW( eIndW.getRawPtr() );
//  }
//  fIS[2] =  Teuchos::rcp( new IndexSpace<O>( EField::W, sIndW, eIndW ) );
//
//  return( fIS );
//}
//
//
//
///// \brief function that creates ScaparIndexSpace
/////
///// by getting values from impact
///// \todo make construction clean without copying
///// \relates IndexSpace
//template<class O=int, int d=3>
//Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<O> > >
//createFullFieldIndexSpaces(
//    const Teuchos::RCP<const StencilWidths<d> >& sW,
//    const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
//    const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
//    bool setImpact=true ){
//
//
//  typedef typename IndexSpace<O>::TO3 TO3;
//
//  Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > > fIS(3);
//
//  TO3 sIndU;
//  TO3 eIndU;
//
//  TO3 sIndV;
//  TO3 eIndV;
//
//  TO3 sIndW;
//  TO3 eIndW;
//
//  for( int i=0; i<3; ++i ) {
//    sIndU[i] = 2            + sW->getLS(i);
//    eIndU[i] = nLoc->get(i) + sW->getLS(i);
//
//    sIndV[i] = 2            + sW->getLS(i);
//    eIndV[i] = nLoc->get(i) + sW->getLS(i);
//
//    sIndW[i] = 2            + sW->getLS(i);
//    eIndW[i] = nLoc->get(i) + sW->getLS(i);
//  }
//
//  // Lower index in x-direction
//  if( bc->BCL_int_[0] > 0 ) {
//    sIndU[0] = 0;
//    sIndV[0] = 1;
//    sIndW[0] = 1;
//  }
//  // lower index in y-direction
//  if( bc->BCL_int_[1] > 0 ) {
//    sIndU[1] = 1;
//    sIndV[1] = 0;
//    sIndW[1] = 1;
//  }
//  // lower index in z-direction
//  if( bc->BCL_int_[2] > 0 ) {
//    sIndU[2] = 1;
//    sIndV[2] = 1;
//    sIndW[2] = 0;
//  }
//
//  // lower index for symmetriBC
//  for( int i=0; i<3; ++i )
//    if( bc->BCL_local_[i]==SymmetryBC ) {
//      sIndU[i] = 1;
//      sIndV[i] = 1;
//      sIndW[i] = 1;
//    }
//
//  // upper index
//  for( int i=0; i<3; ++i )
//    if( bc->BCU_int_[i] > 0 ) {
//      eIndU[i] = nLoc->get(i);
//      eIndV[i] = nLoc->get(i);
//      eIndW[i] = nLoc->get(i);
//    }
//
//  // upper index in x-direction for symmetricBC
//  if( bc->BCU_local_[0]==SymmetryBC ) {
//    eIndU[0] = nLoc->get(0)-1;
//    eIndV[0] = nLoc->get(0);
//    eIndW[0] = nLoc->get(0);
//  }
//
//  // upper index in y-direction for symmetricBC
//  if( bc->BCU_local_[1]==SymmetryBC ) {
//    eIndU[1] = nLoc->get(1);
//    eIndV[1] = nLoc->get(1)-1;
//    eIndW[1] = nLoc->get(1);
//  }
//
//  // upper index in z-direction for symmetricBC
//  if( bc->BCU_local_[2]==SymmetryBC ) {
//    eIndU[2] = nLoc->get(2);
//    eIndV[2] = nLoc->get(2);
//    eIndW[2] = nLoc->get(2)-1;
//  }
//
//
//  if( setImpact ) {
//    VS_set_sIndUB( sIndU.getRawPtr() );
//    VS_set_eIndUB( eIndU.getRawPtr() );
//  }
//  fIS[0] =  Teuchos::rcp( new IndexSpace<O>( EField::U, sIndU, eIndU ) );
//
//  if( setImpact ) {
//    VS_set_sIndVB( sIndV.getRawPtr() );
//    VS_set_eIndVB( eIndV.getRawPtr() );
//  }
//  fIS[1] =  Teuchos::rcp( new IndexSpace<O>( EField::V, sIndV, eIndV ) );
//
//  if( setImpact ) {
//    VS_set_sIndWB( sIndW.getRawPtr() );
//    VS_set_eIndWB( eIndW.getRawPtr() );
//  }
//  fIS[2] =  Teuchos::rcp( new IndexSpace<O>( EField::W, sIndW, eIndW ) );
//
//  return( fIS );
//}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INDEXSPACE_HPP
