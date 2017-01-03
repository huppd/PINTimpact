#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP
#define PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact{


/// \brief local boundary conditions depending on processor grid
///
/// \tparam dimension  as soon as Time is own class dimension ->sd
/// \ingroup SpaceObject
template<int dimension>
class BoundaryConditionsLocal {

protected:

	template< class OT, int dT >
	friend Teuchos::RCP< const BoundaryConditionsLocal<dT> >  createBoudaryConditionsLocal(
			const Teuchos::RCP<const BoundaryConditionsGlobal<dT> >& bcg,
			const Teuchos::RCP<const ProcGrid<OT,dT> >&  pg );

	using TBC3 = const Teuchos::Tuple<EBCType,3>;
	using TO = const Teuchos::Tuple<int,dimension>;

	TO BCL_int_;
	TO BCU_int_;

	BoundaryConditionsLocal(
			EBCType BC1L=DirichletBC,
			EBCType BC1U=DirichletBC,
			EBCType BC2L=DirichletBC,
			EBCType BC2U=DirichletBC,
			EBCType BC3L=DirichletBC,
			EBCType BC3U=DirichletBC ) {

		TBC3 BCL_local_ = Teuchos::tuple( BC1L, BC2L, BC3L );
		TBC3 BCU_local_ = Teuchos::tuple( BC1U, BC2U, BC3U );

		for( int i=0; i<3; ++i ) {
			BCL_int_[i] = static_cast<int>( BCL_local_[i] );
			BCU_int_[i] = static_cast<int>( BCU_local_[i] );
		}
		if( 4==dimension ) {
			BCL_int_[3] = static_cast<int>( NeighborBC );
			BCU_int_[3] = static_cast<int>( NeighborBC );
		}
	};


	BoundaryConditionsLocal( TBC3 BCL_local, TBC3 BCU_local ) {

		for( int i=0; i<3; ++i ) {
			BCL_int_[i] = static_cast<int>( BCL_local[i] );
			BCU_int_[i] = static_cast<int>( BCU_local[i] );
		}
		if( 4==dimension ) {
			BCL_int_[3] = static_cast<int>( NeighborBC );
			BCU_int_[3] = static_cast<int>( NeighborBC );
		}
	}

public:

	/// \name getter
	/// @{ 
	
	constexpr const int& getBCL( const int& dir ) const { return( BCL_int_[dir] ); }
	constexpr const int& getBCU( const int& dir ) const { return( BCU_int_[dir] ); }

	constexpr const int* getBCL() const { return( BCL_int_.getRawPtr() ); }
	constexpr const int* getBCU() const { return( BCU_int_.getRawPtr() ); }

	///  @} 

	/// \brief prints local BC
	///
	/// \param out output stream
	void print( std::ostream& out=std::cout ) const {
		out << "\t--- local BoundaryConditions: ---\n";
		out << "\tlower: " << BCL_int_ << "\n";
		out << "\tupper: " << BCU_int_ << "\n";
	}

}; // end of class BoundaryConditionsLocal




/// \brief creates local boundary conditions
///
/// \tparam O Ordinal
/// \tparam d dimension
/// \param bcg global boundary conditions
/// \param pg processor grid
///
/// \return 
template< class O, int d >
Teuchos::RCP< const BoundaryConditionsLocal<d> >  createBoudaryConditionsLocal(
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
			BCL[i] = static_cast<EBCType>( bcg->getBCL(i) );
		if( pg->getRankU(i)<0 )
			BCU[i] = static_cast<EBCType>( bcg->getBCU(i) );
	}

	return(
			Teuchos::rcp(
				new BoundaryConditionsLocal<d>( BCL, BCU ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP
