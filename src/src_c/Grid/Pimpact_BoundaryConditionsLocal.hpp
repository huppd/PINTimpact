#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP
#define PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_ProcGrid.hpp"



extern "C" {
void fgetBCLoc(
    int& BC_1L_,
    int& BC_1U_,
    int& BC_2L_,
    int& BC_2U_,
    int& BC_3L_,
    int& BC_3U_ );
void fsetBCLoc(
    const int& BC_1L_,
    const int& BC_1U_,
    const int& BC_2L_,
    const int& BC_2U_,
    const int& BC_3L_,
    const int& BC_3U_ );
}



namespace Pimpact{



/// \ingroup domain
//template< int dimension=3>
class BoundaryConditionsLocal {

public:

  friend Teuchos::RCP<BoundaryConditionsLocal> createBoudaryConditionsLocal();

  template< class OT, int dT >
  friend Teuchos::RCP<BoundaryConditionsLocal>  createBoudaryConditionsLocal(
        const Teuchos::RCP< BoundaryConditionsGlobal >& bcg,
        const Teuchos::RCP< ProcGridSize<OT,dT> >&  pgs,
        const Teuchos::RCP< ProcGrid<OT,dT> >&  pg );


  typedef const Teuchos::Tuple<EBCType,3> TBC3;
  typedef const Teuchos::Tuple<int,3> Ti3;

  TBC3 BCL_local_;
  TBC3 BCU_local_;

  Ti3 BCL_int_;
  Ti3 BCU_int_;

protected:

  BoundaryConditionsLocal(
      EBCType BC1L=DirichletBC,
      EBCType BC1U=DirichletBC,
      EBCType BC2L=DirichletBC,
      EBCType BC2U=DirichletBC,
      EBCType BC3L=DirichletBC,
      EBCType BC3U=DirichletBC ):
        BCL_local_( Teuchos::tuple(BC1L, BC2L, BC3L) ),
        BCU_local_( Teuchos::tuple(BC1U, BC2U, BC3U) ) {
    for( int i=0; i<3; ++i ) {
      BCL_int_[i] = (int)( BCL_local_[i] );
      BCU_int_[i] = (int)( BCU_local_[i] );
    }

  };

  BoundaryConditionsLocal( TBC3 BCL_local, TBC3 BCU_local ):
    BCL_local_( BCL_local ),
    BCU_local_( BCU_local ) {
    for( int i=0; i<3; ++i ) {
      BCL_int_[i] = (int)( BCL_local_[i] );
      BCU_int_[i] = (int)( BCU_local_[i] );
    }
  }

public:

  void set_Impact(){
    fsetBCLoc(
        BCL_local_[0],
        BCU_local_[0],
        BCL_local_[1],
        BCU_local_[1],
        BCL_local_[2],
        BCU_local_[2] );
  }

  const int* getBCL() const { return( BCL_int_.getRawPtr() ); }
  const int& getBCL( int i ) const { return( BCL_int_[i] ); }

  const int* getBCU() const { return( BCU_int_.getRawPtr() ); }
  const int& getBCU( int i ) const { return( BCU_int_[i] ); }

  void print( std::ostream& out=std::cout ) const {
    out << "---BoundaryConditionsGlobal: ---\n";
    out << " BCL_local: " << BCL_local_ << "\n";
    out << " BCU_local: " << BCU_local_ << "\n";
  }

}; // end of class BoundaryConditionsLocal





/// \relates BoundaryConditionsLocal
Teuchos::RCP<BoundaryConditionsLocal> createBoudaryConditionsLocal() {
  typedef const Teuchos::Tuple<int,3> TBC3;

  TBC3 BCL;
  TBC3 BCU;

  fgetBCLoc(
    BCL[0],
    BCU[0],
    BCL[1],
    BCU[1],
    BCL[2],
    BCU[2] );

  return(
      Teuchos::rcp(
          new BoundaryConditionsLocal(
              EBCType( BCL[0] ),
              EBCType( BCU[0] ),
              EBCType( BCL[1] ),
              EBCType( BCU[1] ),
              EBCType( BCL[2] ),
              EBCType( BCU[2] )   ) ) );
}


template< class O, int d >
Teuchos::RCP<BoundaryConditionsLocal>  createBoudaryConditionsLocal(
      const Teuchos::RCP< BoundaryConditionsGlobal >& bcg,
      const Teuchos::RCP< ProcGridSize<O,d> >&  pgs,
      const Teuchos::RCP< ProcGrid<O,d> >&  pg ) {

  typedef const Teuchos::Tuple<EBCType,3> TBC3;

  // Default: Neighbor block
  TBC3 BCL = Teuchos::tuple( NeighborBC, NeighborBC, NeighborBC );
  TBC3 BCU = Teuchos::tuple( NeighborBC, NeighborBC, NeighborBC );

  // special case periodic BC with only one block:
  for( int i=0; i<3; ++i ) {
    if( PeriodicBC==bcg->getBCL(i) && 1==pgs->get(i) ){
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
