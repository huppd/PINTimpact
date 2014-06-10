#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP
#define PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



extern "C" {
void fgetBCloc(
    int& BC_1L_,
    int& BC_1U_,
    int& BC_2L_,
    int& BC_2U_,
    int& BC_3L_,
    int& BC_3U_ );
}



namespace Pimpact{

class BoundaryConditionsLocal {

public:

  typedef const Teuchos::Tuple<EBCType,3> TBC3;

//protected:

  TBC3 BCL_local_;
  TBC3 BCU_local_;

//public:

  BoundaryConditionsLocal(
      EBCType BC1L=DirichletBC,
      EBCType BC1U=DirichletBC,
      EBCType BC2L=DirichletBC,
      EBCType BC2U=DirichletBC,
      EBCType BC3L=DirichletBC,
      EBCType BC3U=DirichletBC ):
        BCL_local_( Teuchos::tuple(BC1L, BC2L, BC3L) ),
        BCU_local_( Teuchos::tuple(BC1U, BC2U, BC3U) ) {};

  BoundaryConditionsLocal( TBC3 BCL_local, TBC3 BCU_local ):
    BCL_local_( BCL_local ),
    BCU_local_( BCU_local ) {}

}; // end of class BoundaryConditionsLocal





/// \relates BoundaryConditions
Teuchos::RCP<BoundaryConditionsLocal> createBoudaryConditionsLocal(
    EDomainType dtype = Dirichelt2DChannel ) {

  return( Teuchos::null );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_BOUNDARYCONDITIONSLOCAL_HPP
