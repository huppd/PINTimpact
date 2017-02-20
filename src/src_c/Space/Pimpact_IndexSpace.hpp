#pragma once
#ifndef PIMPACT_INDEXSPACE_HPP
#define PIMPACT_INDEXSPACE_HPP


#include <iostream>

#include "mpi.h"

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_BoundaryConditionsLocal.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief class that stores neccessary lower and upper indexes for \c ScalarField for different FieldTypes
/// and including Boundaries or excluding them.
///
/// Index         | BoundaryConditionslocal>0 | BoundaryConditionslocal<=0 | symmetryBoundaryConditions
/// ------------- | --------------------------| ---------------------------| ---------------------------
/// Ind           | 1                         | 1                          | 1
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
/// \tparam OrdinalT
/// \tparam dimension as soon as Time is own class ->spatial dimension
/// \ingroup SpaceObject
template<class OrdinalT, int dimension>
class IndexSpace {

  template<class O, int sd, int d, int dimNC>
  friend Teuchos::RCP<const IndexSpace<O,d> >
  createIndexSpace(
      const Teuchos::RCP<const StencilWidths<d,dimNC> >& sW,
      const Teuchos::RCP<const GridSizeLocal<O,sd,d> >& gL,
      const Teuchos::RCP<const BoundaryConditionsLocal<d> >& bc,
		 	const Teuchos::RCP<const ProcGrid<O,d> >& pG );

public:

  using TO = const Teuchos::Tuple<OrdinalT,dimension>;

  using TTO = const Teuchos::Tuple< Teuchos::Tuple<OrdinalT, dimension>, 3 >;

protected:

  TO sIndS_;
  TO eIndS_;

  TTO sIndU_;
  TTO eIndU_;

  TTO sIndUB_;
  TTO eIndUB_;

  Teuchos::Tuple<OrdinalT,dimension> shift_;

  /// \brief constructor
	///
	/// \tparam dimNC stencil widths
  /// \param sW with of differenct stencils
  /// \param gridSizeLocal amount of grid points stored localy on this node
  /// \param bc local boundary conditions
	/// \param procGrid processor grid
  template<int sd, int dimNC>
  IndexSpace(
      const Teuchos::RCP<const StencilWidths<dimension,dimNC> >& sW,
      const Teuchos::RCP<const GridSizeLocal<OrdinalT,sd,dimension> >& gridSizeLocal,
      const Teuchos::RCP<const BoundaryConditionsLocal<dimension> >& bc,
		 	const Teuchos::RCP<const ProcGrid<OrdinalT,dimension> >& procGrid ) {

    // ------------------init IndS_--------------------
    for( int i=0; i<3; ++i ) {
      sIndS_[i] = 2            + sW->getLS(i);
      eIndS_[i] = gridSizeLocal->get(i) + sW->getLS(i);
      if( bc->getBCL(i) > 0 )
        sIndS_[i] = 1;
      if( bc->getBCU(i) > 0 )
        eIndS_[i] = gridSizeLocal->get(i);
      if( bc->getBCL(i)==BC::Symmetry )
        sIndS_[i] = 1;
      if( bc->getBCU(i)==BC::Symmetry )
        eIndS_[i] = gridSizeLocal->get(i);
    }

    // time direction automatically periodic BC
    if( 4==dimension ) {
			if( sW->spectralT() ) {

				OrdinalT nl  = (gridSizeLocal->get(3)+1)/procGrid->getNP(3);
				OrdinalT rem = (gridSizeLocal->get(3)+1)%procGrid->getNP(3);
				int rank = procGrid->getIB(3)-1;
				OrdinalT sI =       rank   *nl + (( rank   <rem)? rank   :rem);
				OrdinalT eI = -1 + (rank+1)*nl + (((rank+1)<rem)?(rank+1):rem);

//				std::cout << "\tnl: " << nl << "\trem: " << rem << "\trank: " << rank << "\tsI: " << sI << "\teI: " << eI << "\n";
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
				OrdinalT sI = 0 - sW->getBL(3);
				OrdinalT eI = gridSizeLocal->get(3) - sW->getBL(3)-1;
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
    if( bc->getBCL(0) > 0 ) {
      sIndU_[0][0] = 1;
      sIndU_[1][0] = 2;
      sIndU_[2][0] = 2;
    }
    // lower index in y-direction
    if( bc->getBCL(1) > 0 ) {
      sIndU_[0][1] = 2;
      sIndU_[1][1] = 1;
      sIndU_[2][1] = 2;
    }
    // lower index in z-direction
    if( bc->getBCL(2) > 0 ) {
      sIndU_[0][2] = 2;
      sIndU_[1][2] = 2;
      sIndU_[2][2] = 1;
    }

    // lower index for symmetriBC
    for( int i=0; i<3; ++i )
      if( bc->getBCL(i)==BC::Symmetry ) {
        sIndU_[0][i] = 1;
        sIndU_[1][i] = 1;
        sIndU_[2][i] = 1;
      }

    // upper index
    for( int i=0; i<3; ++i )
      if( bc->getBCU(i) > 0 ) {
        eIndU_[0][i] = gridSizeLocal->get(i)-1;
        eIndU_[1][i] = gridSizeLocal->get(i)-1;
        eIndU_[2][i] = gridSizeLocal->get(i)-1;
      }

    // upper index in x-direction for symmetricBC
    if( bc->getBCU(0)==BC::Symmetry ) {
      eIndU_[0][0] = gridSizeLocal->get(0)-1;
      eIndU_[1][0] = gridSizeLocal->get(0);
      eIndU_[2][0] = gridSizeLocal->get(0);
    }

    // upper index in y-direction for symmetricBC
    if( bc->getBCU(1)==BC::Symmetry ) {
      eIndU_[0][1] = gridSizeLocal->get(1);
      eIndU_[1][1] = gridSizeLocal->get(1)-1;
      eIndU_[2][1] = gridSizeLocal->get(1);
    }

    // upper index in z-direction for symmetricBC
    if( bc->getBCU(2)==BC::Symmetry ) {
      eIndU_[0][2] = gridSizeLocal->get(2);
      eIndU_[1][2] = gridSizeLocal->get(2);
      eIndU_[2][2] = gridSizeLocal->get(2)-1;
    }

		//
    // --- init IndUB_ -------------------------------
		// 
    for( int i=0; i<3; ++i ) {
      sIndUB_[0][i] = 2            + sW->getLS(i);
      eIndUB_[0][i] = gridSizeLocal->get(i) + sW->getLS(i);

      sIndUB_[1][i] = 2            + sW->getLS(i);
      eIndUB_[1][i] = gridSizeLocal->get(i) + sW->getLS(i);

      sIndUB_[2][i] = 2            + sW->getLS(i);
      eIndUB_[2][i] = gridSizeLocal->get(i) + sW->getLS(i);
    }

    // Lower index in x-direction
    if( bc->getBCL(0) > 0 ) {
      sIndUB_[0][0] = 0;
      sIndUB_[1][0] = 1;
      sIndUB_[2][0] = 1;
    }
    // lower index in y-direction
    if( bc->getBCL(1) > 0 ) {
      sIndUB_[0][1] = 1;
      sIndUB_[1][1] = 0;
      sIndUB_[2][1] = 1;
    }
    // lower index in z-direction
    if( bc->getBCL(2) > 0 ) {
      sIndUB_[0][2] = 1;
      sIndUB_[1][2] = 1;
      sIndUB_[2][2] = 0;
    }

    // lower index for symmetriBC
    for( int i=0; i<3; ++i )
      if( bc->getBCL(i)==BC::Symmetry ) {
        sIndUB_[0][i] = 1;
        sIndUB_[1][i] = 1;
        sIndUB_[2][i] = 1;
      }

    // upper index
    for( int i=0; i<3; ++i )
      if( bc->getBCU(i) > 0 ) {
        eIndUB_[0][i] = gridSizeLocal->get(i);
        eIndUB_[1][i] = gridSizeLocal->get(i);
        eIndUB_[2][i] = gridSizeLocal->get(i);
      }

    // upper index in x-direction for symmetricBC
    if( bc->getBCU(0)==BC::Symmetry ) {
      eIndUB_[0][0] = gridSizeLocal->get(0)-1;
      eIndUB_[1][0] = gridSizeLocal->get(0);
      eIndUB_[2][0] = gridSizeLocal->get(0);
    }

    // upper index in y-direction for symmetricBC
    if( bc->getBCU(1)==BC::Symmetry ) {
      eIndUB_[0][1] = gridSizeLocal->get(1);
      eIndUB_[1][1] = gridSizeLocal->get(1)-1;
      eIndUB_[2][1] = gridSizeLocal->get(1);
    }

    // upper index in z-direction for symmetricBC
    if( bc->getBCU(2)==BC::Symmetry ) {
      eIndUB_[0][2] = gridSizeLocal->get(2);
      eIndUB_[1][2] = gridSizeLocal->get(2);
      eIndUB_[2][2] = gridSizeLocal->get(2)-1;
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

  constexpr const OrdinalT* sInd( const F& ft ) const {
		return(
				( F::S==ft )?
				sIndS_.getRawPtr():
				sIndU_[static_cast<int>(ft)].getRawPtr() );
	}
  constexpr const OrdinalT* eInd( const F& ft ) const {
		return(
				( F::S==ft )?
				eIndS_.getRawPtr():
				eIndU_[static_cast<int>(ft)].getRawPtr() );
  }

  constexpr const OrdinalT* sIndB( const F& ft ) const {
    return(
				( F::S==ft )?
				sIndS_.getRawPtr():
				sIndUB_[static_cast<int>(ft)].getRawPtr() );
  }
  constexpr const OrdinalT* eIndB( const F& ft ) const {
    return(
				( F::S== ft )?
				eIndS_.getRawPtr():
				eIndUB_[static_cast<int>(ft)].getRawPtr() );
  }

  constexpr const OrdinalT& sInd( const F& ft, const int& dir ) const {
    return(
				( F::S==ft )?
				sIndS_[dir]:
				sIndU_[static_cast<int>(ft)][dir] );
  }
  constexpr const OrdinalT& eInd(  const F& ft, const int& dir ) const {
		return(
				( F::S==ft )?
				eIndS_[dir]:
				eIndU_[static_cast<int>(ft)][dir] );
  }

  constexpr const OrdinalT& sIndB( const F& ft, const int& dir ) const {
		return(
				( F::S==ft )?
				sIndS_[dir]:
				sIndUB_[static_cast<int>(ft)][dir] );
  }
  constexpr const OrdinalT& eIndB( const F& ft, const int& dir ) const {
    return(
				( F::S==ft )?
				eIndS_[dir]:
				eIndUB_[static_cast<int>(ft)][dir] );
  }

	constexpr const int& getShift( const int& i ) const { return( shift_[i] ); }
	constexpr const int* getShift() const { return( shift_.getRawPtr() ); }

  void print( std::ostream& out=std::cout ) const {
    out << "\t---IndexSpace: ---\n";
    out << "\tfieldType: S\n";
    out << "\tsInd: " << sIndS_ << "\n";
    out << "\teInd: " << eIndS_ << "\n";
    for( int field=0; field<3; ++field ) {
      out << "\tinner field: " << static_cast<F>(field) << "\n";
      out << "\tsInd: " << sIndU_[field] << "\n";
      out << "\teInd: " << eIndU_[field] << "\n";
    }
    for( int field=0; field<3; ++field ) {
      out << "\tfull field: " << static_cast<F>(field) << "\n";
      out << "\tsInd: " << sIndUB_[field] << "\n";
      out << "\teInd: " << eIndUB_[field] << "\n";
    }
		out << "\toffset: " << shift_ << "\n";
  }

}; // end of class IndexSpace



template<class OT, int sd, int d, int dimNC>
Teuchos::RCP<const IndexSpace<OT,d> >
createIndexSpace(
    const Teuchos::RCP<const StencilWidths<d,dimNC> >& sW,
    const Teuchos::RCP<const GridSizeLocal<OT,sd,d> >& gSL,
    const Teuchos::RCP<const BoundaryConditionsLocal<d> >& bc,
		const Teuchos::RCP<const ProcGrid<OT,d> >& pG ) {
  return(
      Teuchos::rcp(
          new IndexSpace<OT,d>( sW, gSL, bc,pG )
      )
  );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INDEXSPACE_HPP
