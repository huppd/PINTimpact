#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP
#define PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"




namespace Pimpact{



/// \ingroup domain
template<int dim>
class BoundaryConditionsGlobal {

public:

  typedef const Teuchos::Tuple<EBCType,dim> TBC3;
  typedef const Teuchos::Tuple<int,dim> Ti3;


	template<int d>
  friend Teuchos::RCP<const BoundaryConditionsGlobal<d> > createBoudaryConditionsGlobal();

	template<int d>
  friend Teuchos::RCP<const BoundaryConditionsGlobal<d> > createBoudaryConditionsGlobal( EDomainType dtype );

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
      EBCType BC3U=DirichletBC ) {

		BCL_global_[0] = BC1L;
		BCU_global_[0] = BC1U;
		BCL_global_[1] = BC2L;
		BCU_global_[1] = BC2U;
		BCL_global_[2] = BC3L;
		BCU_global_[2] = BC3U;

    for( int i=0; i<3; ++i ) {
      BCL_int_[i] = (int)( BCL_global_[i] );
      BCU_int_[i] = (int)( BCU_global_[i] );
    };
		if( 4==dim ) {
      BCL_global_[3] = PeriodicBC;
      BCU_global_[3] = PeriodicBC;
      BCL_int_[3] = (int)( PeriodicBC );
      BCU_int_[3] = (int)( PeriodicBC );
		}
  }

  BoundaryConditionsGlobal( TBC3 BCL_global, TBC3 BCU_global ):
    BCL_global_( BCL_global ), BCU_global_( BCU_global ) {
    for( int i=0; i<3; ++i ) {
      BCL_int_[i] = (int)( BCL_global_[i] );
      BCU_int_[i] = (int)( BCU_global_[i] );
    };
		if( 4==dim ) {
      BCL_global_[3] = PeriodicBC;
      BCU_global_[3] = PeriodicBC;
      BCL_int_[3] = (int)( PeriodicBC );
      BCU_int_[3] = (int)( PeriodicBC );
		}
  }

public:

  EBCType getBCL( int dir ) const { return( BCL_global_[dir] ); }

	EBCType getBCU( int dir ) const { return( BCU_global_[dir] ); }


  const int* getBCL() const { return( BCL_int_.getRawPtr() ); }
  //  const int& getBCL( int i ) const { return( BCL_int_[i] ); }

  const int* getBCU() const { return( BCU_int_.getRawPtr() ); }
  //  const int& getBCU( int i ) const { return( BCU_int_[i] ); }


	const Teuchos::Tuple<int,dim> periodic() const {
		Teuchos::Tuple<int,dim> periodic;
    for( int i=0; i<3; ++i ) {
      if( getBCL(i)==PeriodicBC )
        periodic[i] = 1;
      else
        periodic[i] = 0;
    }
    if( 4==dim ) periodic[3] = 1;
		return( periodic );

	}


  void print( std::ostream& out=std::cout ) const {
    out << "---BoundaryConditionsGlobal: ---\n";
    out << " BCL_global: " << BCL_global_ << "\n";
    out << " BCU_global: " << BCU_global_ << "\n";

  }

}; // end of class BoundaryConditionsGlobal




/// \relates BoundaryConditionsGlobal
template<int dim=3>
Teuchos::RCP<const BoundaryConditionsGlobal<dim> >
createBoudaryConditionsGlobal(
    EDomainType dtype ) {

  switch( dtype ) {
  case AllDirichlet:
    return(
        Teuchos::rcp(
            new BoundaryConditionsGlobal<dim>() ) );
    break;
  case Dirichelt2DChannel:
    return(
        Teuchos::rcp(
            new BoundaryConditionsGlobal<dim>(
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
            new BoundaryConditionsGlobal<dim>(
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
            new BoundaryConditionsGlobal<dim>(
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
            new BoundaryConditionsGlobal<dim>(
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
            new BoundaryConditionsGlobal<dim>(
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
