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
/// \note if one would remove ls_ one could decouple with \c StencilWidths
/// \ingroup Space
template<class Ordinal, int dimension>
class IndexSpace {

  template<class O, int d, int dimNC>
  friend Teuchos::RCP<const IndexSpace<O,d> >
  createIndexSpace(
      const Teuchos::RCP<const StencilWidths<d,dimNC> >& sW,
      const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
      const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
      bool setImpact=true );

public:

  typedef Teuchos::Tuple<Ordinal,dimension> TO;

  typedef Teuchos::Tuple< Teuchos::Tuple<Ordinal, dimension>, 3 > TTO;

protected:

  TO sIndS_;
  TO eIndS_;

  TTO sIndU_;
  TTO eIndU_;

  TTO sIndUB_;
  TTO eIndUB_;

  /// \brief constructor
  /// \param sW with of differenct stencils
  /// \param nLoc amount of grid points stored localy on this node
  /// \param bc local boundary conditions
  template<int dimNC>
  IndexSpace(
      const Teuchos::RCP<const StencilWidths<dimension,dimNC> >& sW,
      const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >& nLoc,
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

    if( 4==dimension ) {
      sIndS_[3] = 0 - sW->getBL(3);
      eIndS_[3] = nLoc->get(3) + sW->getBU(3) - sW->getBL(3);
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

    if( 4==dimension ) {
      for( int i=0; i<3; ++i ) {
        sIndU_[i][3] = 0 - sW->getBL(3);
        eIndU_[i][3] = nLoc->get(3) + sW->getBU(3) - sW->getBL(3);
      }
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

    // time direction automatically periodic BC
    if( 4==dimension ) {
      for( int i=0; i<3; ++i ) {
        sIndUB_[i][3] = 0 - sW->getBL(3);
        eIndUB_[i][3] = nLoc->get(3) + sW->getBU(3) - sW->getBL(3);
      }
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
      return( sIndS_.getRawPtr() );
    else
      return( sIndU_[fieldType].getRawPtr() );
  }
  const Ordinal* eInd( int fieldType ) const {
    if( EField::S == (EField)fieldType )
      return( eIndS_.getRawPtr() );
    else
      return( eIndU_[fieldType].getRawPtr() );
  }

  const Ordinal* sIndB( int fieldType ) const {
    if( EField::S == (EField)fieldType )
      return( sIndS_.getRawPtr() );
    else
      return( sIndUB_[fieldType].getRawPtr()  );
  }
  const Ordinal* eIndB( int fieldType ) const {
    if( EField::S == (EField)fieldType )
      return( eIndS_.getRawPtr() );
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


}; // end of class IndexSpace



template<class O, int d, int dimNC>
Teuchos::RCP<const IndexSpace<O,d> >
createIndexSpace(
    const Teuchos::RCP<const StencilWidths<d,dimNC> >& sW,
    const Teuchos::RCP<const GridSizeLocal<O,d> >& nLoc,
    const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
    bool setImpact=true ) {
  return(
      Teuchos::rcp(
          new IndexSpace<O,d>( sW,nLoc,bc )
      )
  );

}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INDEXSPACE_HPP
