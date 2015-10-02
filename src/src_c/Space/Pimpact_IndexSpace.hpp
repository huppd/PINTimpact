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
#include "Pimpact_ProcGrid.hpp"



namespace Pimpact {



/// \brief class that stores neccessary lower and upper indexes for \c ScalarField for different FieldTypes
/// and including Boundaries or excluding them.
///
/// Index         | BoundaryConditionslocal>0 | BoundaryConditionslocal<=0 | symmetryBoundaryConditions
/// ------------- | --------------------------| ---------------------------| ---------------------------
/// Ind          | 1                         | 1                          | 1
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
/// \ingroup SpaceObject
template<class Ordinal, int dimension>
class IndexSpace {

  template<class O, int d, int dimNC>
  friend Teuchos::RCP<const IndexSpace<O,d> >
  createIndexSpace(
      const Teuchos::RCP<const StencilWidths<d,dimNC> >& sW,
      const Teuchos::RCP<const GridSizeLocal<O,d> >& gL,
      const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
		 	const Teuchos::RCP<const ProcGrid<O,d> >& pG );

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

  Teuchos::Tuple<Ordinal,dimension> shift_;

  /// \brief constructor
  /// \param sW with of differenct stencils
  /// \param gridSizeLocal amount of grid points stored localy on this node
  /// \param bc local boundary conditions
  template<int dimNC>
  IndexSpace(
      const Teuchos::RCP<const StencilWidths<dimension,dimNC> >& sW,
      const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
      const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
		 	const Teuchos::RCP<const ProcGrid<Ordinal,dimension> >& procGrid ) {

    // ------------------init IndS_--------------------
    for( int i=0; i<3; ++i ) {
      sIndS_[i] = 2            + sW->getLS(i);
      eIndS_[i] = gridSizeLocal->get(i) + sW->getLS(i);
      if( bc->BCL_int_[i] > 0 )
        sIndS_[i] = 1;
      if( bc->BCU_int_[i] > 0 )
        eIndS_[i] = gridSizeLocal->get(i);
      if( bc->BCL_local_[i]==SymmetryBC )
        sIndS_[i] = 1;
      if( bc->BCU_local_[i]==SymmetryBC )
        eIndS_[i] = gridSizeLocal->get(i);
    }

    // time direction automatically periodic BC
    if( 4==dimension ) {
			if( sW->spectralT() ) {
				Ordinal nl  = (gridSizeLocal->get(3)+1)/procGrid->getNP(3);
				Ordinal rem = (gridSizeLocal->get(3)+1)%procGrid->getNP(3);
				int rank = procGrid->getIB(3)-1;
				Ordinal sI = -1 +  rank   *nl + (( rank   <rem)? rank   :rem);
				Ordinal eI = -1 + (rank+1)*nl + (((rank+1)<rem)?(rank+1):rem);
				std::cout << "\tnl: " << nl << "\trem: " << rem << "\trank: " << rank << "\tsI: " << sI << "\teI: " << eI << "\n";
				sIndS_[3] = sI;
				eIndS_[3] = eI;
				for( int i=0; i<3; ++i ) {
					sIndU_[i][3] = sI;
					eIndU_[i][3] = eI;
				}
				for( int i=0; i<3; ++i ) {
					sIndUB_[i][3] = sI;
					eIndUB_[i][3] = eI;
				}
			}
			else{
				Ordinal sI = 0 - sW->getBL(3);
				Ordinal eI = gridSizeLocal->get(3) + sW->getBU(3) - sW->getBL(3);
				sIndS_[3] = sI;
				eIndS_[3] = eI;
				for( int i=0; i<3; ++i ) {
					sIndU_[i][3] = sI;
					eIndU_[i][3] = eI;
				}
				for( int i=0; i<3; ++i ) {
					sIndUB_[i][3] = sI;
					eIndUB_[i][3] = eI;
				}
			}
    }

		//
    // --- init IndU_ -------------------------------
		//
    for( int field=0; field<3; ++field )
      for( int dir=0; dir<3; ++dir ) {
        sIndU_[field][dir] = 2              + sW->getLS(dir);
        eIndU_[field][dir] = gridSizeLocal->get(dir) + sW->getLS(dir);
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
        eIndU_[U][i] = gridSizeLocal->get(i)-1;
        eIndU_[V][i] = gridSizeLocal->get(i)-1;
        eIndU_[W][i] = gridSizeLocal->get(i)-1;
      }

    // upper index in x-direction for symmetricBC
    if( bc->BCU_local_[0]==SymmetryBC ) {
      eIndU_[U][0] = gridSizeLocal->get(0)-1;
      eIndU_[V][0] = gridSizeLocal->get(0);
      eIndU_[W][0] = gridSizeLocal->get(0);
    }

    // upper index in y-direction for symmetricBC
    if( bc->BCU_local_[1]==SymmetryBC ) {
      eIndU_[U][1] = gridSizeLocal->get(1);
      eIndU_[V][1] = gridSizeLocal->get(1)-1;
      eIndU_[W][1] = gridSizeLocal->get(1);
    }

    // upper index in z-direction for symmetricBC
    if( bc->BCU_local_[2]==SymmetryBC ) {
      eIndU_[U][2] = gridSizeLocal->get(2);
      eIndU_[V][2] = gridSizeLocal->get(2);
      eIndU_[W][2] = gridSizeLocal->get(2)-1;
    }

		//
    // --- init IndUB_ -------------------------------
		// 
    for( int i=0; i<3; ++i ) {
      sIndUB_[U][i] = 2            + sW->getLS(i);
      eIndUB_[U][i] = gridSizeLocal->get(i) + sW->getLS(i);

      sIndUB_[V][i] = 2            + sW->getLS(i);
      eIndUB_[V][i] = gridSizeLocal->get(i) + sW->getLS(i);

      sIndUB_[W][i] = 2            + sW->getLS(i);
      eIndUB_[W][i] = gridSizeLocal->get(i) + sW->getLS(i);
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
        eIndUB_[U][i] = gridSizeLocal->get(i);
        eIndUB_[V][i] = gridSizeLocal->get(i);
        eIndUB_[W][i] = gridSizeLocal->get(i);
      }

    // upper index in x-direction for symmetricBC
    if( bc->BCU_local_[0]==SymmetryBC ) {
      eIndUB_[U][0] = gridSizeLocal->get(0)-1;
      eIndUB_[V][0] = gridSizeLocal->get(0);
      eIndUB_[W][0] = gridSizeLocal->get(0);
    }

    // upper index in y-direction for symmetricBC
    if( bc->BCU_local_[1]==SymmetryBC ) {
      eIndUB_[U][1] = gridSizeLocal->get(1);
      eIndUB_[V][1] = gridSizeLocal->get(1)-1;
      eIndUB_[W][1] = gridSizeLocal->get(1);
    }

    // upper index in z-direction for symmetricBC
    if( bc->BCU_local_[2]==SymmetryBC ) {
      eIndUB_[U][2] = gridSizeLocal->get(2);
      eIndUB_[V][2] = gridSizeLocal->get(2);
      eIndUB_[W][2] = gridSizeLocal->get(2)-1;
    }


		//
    // --- computes index offset -----------------------------------
		//
    for( int i=0; i<3; ++i )
      shift_[i] = (procGrid->getIB(i)-1)*( gridSizeLocal->get(i)-1 );
    if( 4==dimension )
      shift_[3] = (procGrid->getIB(3)-1)*( gridSizeLocal->get(3) );

  }

public:

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
      return( sIndUB_[fieldType].getRawPtr()[dir]  );
  }
  const Ordinal& eIndB( int fieldType, int dir ) const {
    if( EField::S == (EField)fieldType )
      return( eIndS_[dir]  );
    else
      return( eIndUB_[fieldType][dir]  );
  }

	const int& getShift( int i ) const { return( shift_[i] ); }
	const int* getShift() const { return( shift_.getRawPtr() ); }

  void print( std::ostream& out=std::cout ) const {
    out << "\t---IndexSpace: ---\n";
    out << "\tfieldType: S\n";
    out << "\tsInd: " << sIndS_ << "\n";
    out << "\teInd: " << eIndS_ << "\n";
    for( int field=0; field<3; ++field ) {
      out << "\tinner field: " << field << "\n";
      out << "\tsInd: " << sIndU_[field] << "\n";
      out << "\teInd: " << eIndU_[field] << "\n";
    }
    for( int field=0; field<3; ++field ) {
      out << "\tfull field: " << field << "\n";
      out << "\tsInd: " << sIndUB_[field] << "\n";
      out << "\teInd: " << eIndUB_[field] << "\n";
    }
		out << "\toffset: " << shift_ << "\n";
  }

}; // end of class IndexSpace



template<class O, int d, int dimNC>
Teuchos::RCP<const IndexSpace<O,d> >
createIndexSpace(
    const Teuchos::RCP<const StencilWidths<d,dimNC> >& sW,
    const Teuchos::RCP<const GridSizeLocal<O,d> >& gSL,
    const Teuchos::RCP<const BoundaryConditionsLocal>& bc,
		const Teuchos::RCP<const ProcGrid<O,d> >& pG ) {
  return(
      Teuchos::rcp(
          new IndexSpace<O,d>( sW, gSL, bc,pG )
      )
  );

}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INDEXSPACE_HPP
