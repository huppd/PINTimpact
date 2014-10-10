#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP
#define PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



extern "C" {
void fgetBC(
    int& BC_1L_,
    int& BC_1U_,
    int& BC_2L_,
    int& BC_2U_,
    int& BC_3L_,
    int& BC_3U_ );
void fsetBC(
    const int& BC_1L_,
    const int& BC_1U_,
    const int& BC_2L_,
    const int& BC_2U_,
    const int& BC_3L_,
    const int& BC_3U_ );
}



namespace Pimpact{



/// \ingroup domain
class BoundaryConditionsGlobal {

public:

  typedef const Teuchos::Tuple<EBCType,3> TBC3;
  typedef const Teuchos::Tuple<int,3> Ti3;


  friend Teuchos::RCP<BoundaryConditionsGlobal> createBoudaryConditionsGlobal();

  friend Teuchos::RCP<BoundaryConditionsGlobal> createBoudaryConditionsGlobal( EDomainType dtype );

protected:

  TBC3 BCL_global_;
  TBC3 BCU_global_;

  Ti3 BCL_int_;
  Ti3 BCU_int_;

  BoundaryConditionsGlobal(
      EBCType BC1L=DirichletBC,
      EBCType BC1U=DirichletBC,
      EBCType BC2L=DirichletBC,
      EBCType BC2U=DirichletBC,
      EBCType BC3L=DirichletBC,
      EBCType BC3U=DirichletBC ):
        BCL_global_( Teuchos::tuple(BC1L, BC2L, BC3L) ),
        BCU_global_( Teuchos::tuple(BC1U, BC2U, BC3U) ) {
    for( int i=0; i<3; ++i ) {
         BCL_int_[i] = (int)( BCL_global_[i] );
         BCU_int_[i] = (int)( BCU_global_[i] );
    };
  }

  BoundaryConditionsGlobal( TBC3 BCL_global, TBC3 BCU_global ):
    BCL_global_( BCL_global ), BCU_global_( BCU_global ) {
    for( int i=0; i<3; ++i ) {
         BCL_int_[i] = (int)( BCL_global_[i] );
         BCU_int_[i] = (int)( BCU_global_[i] );
    };
  }

public:

  EBCType getBCL( int dim ) const { return( BCL_global_[dim] ); }

  EBCType getBCU( int dim ) const {
    return( BCU_global_[dim] );
  }


  const int* getBCL() const { return( BCL_int_.getRawPtr() ); }
//  const int& getBCL( int i ) const { return( BCL_int_[i] ); }

  const int* getBCU() const { return( BCU_int_.getRawPtr() ); }
//  const int& getBCU( int i ) const { return( BCU_int_[i] ); }

  void set_Impact(){
    fsetBC(
        BCL_global_[0],
        BCU_global_[0],
        BCL_global_[1],
        BCU_global_[1],
        BCL_global_[2],
        BCU_global_[2] );
  };

  void print( std::ostream& out=std::cout ) const {
    out << "---BoundaryConditionsGlobal: ---\n";
    out << " BCL_global: " << BCL_global_ << "\n";
    out << " BCU_global: " << BCU_global_ << "\n";

  }

}; // end of class BoundaryConditionsGlobal



/// \relates BoundaryConditionsGlobal
Teuchos::RCP<BoundaryConditionsGlobal>
createBoudaryConditionsGlobal() {

  typedef const Teuchos::Tuple<int,3> TBC3;

  TBC3 BCL;
  TBC3 BCU;

  fgetBC(
    BCL[0],
    BCU[0],
    BCL[1],
    BCU[1],
    BCL[2],
    BCU[2] );

  return(
      Teuchos::rcp(
          new BoundaryConditionsGlobal(
              EBCType( BCL[0] ),
              EBCType( BCU[0] ),
              EBCType( BCL[1] ),
              EBCType( BCU[1] ),
              EBCType( BCL[2] ),
              EBCType( BCU[2] )   ) ) );
}


/// \relates BoundaryConditionsGlobal
Teuchos::RCP<BoundaryConditionsGlobal> createBoudaryConditionsGlobal(
    EDomainType dtype ) {

  switch( dtype ) {
  case AllDirichlet:
    return(
        Teuchos::rcp(
            new BoundaryConditionsGlobal() ) );
    break;
  case Dirichelt2DChannel:
    return(
        Teuchos::rcp(
            new BoundaryConditionsGlobal(
                DirichletBC,
                DirichletBC,
                DirichletBC,
                DirichletBC,
                PeriodicBC,
                PeriodicBC ) ) );
    break;
  case Periodic2DChannel:
    return(
        Teuchos::rcp(
            new BoundaryConditionsGlobal(
                PeriodicBC,
                PeriodicBC,
                DirichletBC,
                DirichletBC,
                PeriodicBC,
                PeriodicBC ) ) );
    break;
  case AllNeumann2D:
    return(
        Teuchos::rcp(
            new BoundaryConditionsGlobal(
                NeumannBC,
                NeumannBC,
                NeumannBC,
                NeumannBC,
                PeriodicBC,
                PeriodicBC ) ) );
    break;
  case AllPeriodic:
    return(
        Teuchos::rcp(
            new BoundaryConditionsGlobal(
                PeriodicBC,
                PeriodicBC,
                PeriodicBC,
                PeriodicBC,
                PeriodicBC,
                PeriodicBC ) ) );
    break;
  case Neumann1Periodic2:
    return(
        Teuchos::rcp(
            new BoundaryConditionsGlobal(
                NeumannBC,
                NeumannBC,
                PeriodicBC,
                PeriodicBC,
                PeriodicBC,
                PeriodicBC ) ) );
    break;
  default:
    std::cout << "!!!Warning: unkown EDomainType:\t" <<dtype<<"\t!!!\n";
    return( Teuchos::null );
  }
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP
