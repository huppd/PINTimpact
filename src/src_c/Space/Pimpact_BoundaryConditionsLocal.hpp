#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP
#define PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{


/// \ingroup SpaceObject
//template< int dimension=3>
class BoundaryConditionsLocal {

protected:

	friend Teuchos::RCP<const BoundaryConditionsLocal> createBoudaryConditionsLocal();

	template< class OT, int dT >
	friend Teuchos::RCP<const BoundaryConditionsLocal>  createBoudaryConditionsLocal(
			const Teuchos::RCP<const BoundaryConditionsGlobal<dT> >& bcg,
			const Teuchos::RCP<const ProcGrid<OT,dT> >&  pg );

	typedef const Teuchos::Tuple<EBCType,3> TBC3;
	typedef const Teuchos::Tuple<int,3> Ti3;

	Ti3 BCL_int_;
	Ti3 BCU_int_;

	BoundaryConditionsLocal(
			EBCType BC1L=DirichletBC,
			EBCType BC1U=DirichletBC,
			EBCType BC2L=DirichletBC,
			EBCType BC2U=DirichletBC,
			EBCType BC3L=DirichletBC,
			EBCType BC3U=DirichletBC ) {

		TBC3 BCL_local_ = Teuchos::tuple(BC1L, BC2L, BC3L);
		TBC3 BCU_local_ = Teuchos::tuple(BC1U, BC2U, BC3U);

		for( int i=0; i<3; ++i ) {
			BCL_int_[i] = static_cast<int>( BCL_local_[i] );
			BCU_int_[i] = static_cast<int>( BCU_local_[i] );
		}

	};

	BoundaryConditionsLocal( TBC3 BCL_local, TBC3 BCU_local ) {

		for( int i=0; i<3; ++i ) {
			BCL_int_[i] = static_cast<int>( BCL_local[i] );
			BCU_int_[i] = static_cast<int>( BCU_local[i] );
		}

	}

public:

	/// \name getter
	/// @{ 
	
  EBCType getBCL( const int& dir ) const { return( static_cast<EBCType>(BCL_int_[dir]) ); }
	EBCType getBCU( const int& dir ) const { return( static_cast<EBCType>(BCU_int_[dir]) ); }

  const int* getBCL() const { return( BCL_int_.getRawPtr() ); }
  const int* getBCU() const { return( BCU_int_.getRawPtr() ); }

	///  @} 

  void print( std::ostream& out=std::cout ) const {
    out << "\t--- local BoundaryConditions: ---\n";
    out << "\tlower: " << BCL_int_ << "\n";
    out << "\tupper: " << BCU_int_ << "\n";
  }

}; // end of class BoundaryConditionsLocal



template< class O, int d >
Teuchos::RCP<const BoundaryConditionsLocal>  createBoudaryConditionsLocal(
      const Teuchos::RCP<const BoundaryConditionsGlobal<d> >& bcg,
      const Teuchos::RCP<const ProcGrid<O,d> >&  pg ) {

  typedef const Teuchos::Tuple<EBCType,3> TBC3;

  // Default: Neighbor block
  TBC3 BCL = Teuchos::tuple( NeighborBC, NeighborBC, NeighborBC );
  TBC3 BCU = Teuchos::tuple( NeighborBC, NeighborBC, NeighborBC );

  // special case periodic BC with only one block:
  for( int i=0; i<3; ++i ) {
    if( PeriodicBC==bcg->getBCL(i) && 1==pg->getNP(i) ){
      BCL[i] = PeriodicBC;
      BCU[i] = PeriodicBC;
    }
  }

  // boundary condition on procgrid boundar
  for( int i=0; i<3; ++i ) {
    if( pg->getRankL(i)<0 )
      BCL[i] = bcg->getBCL(i);
    if( pg->getRankU(i)<0 )
      BCU[i] = bcg->getBCU(i);
  }

  return(
      Teuchos::rcp(
          new BoundaryConditionsLocal( BCL, BCU ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP
