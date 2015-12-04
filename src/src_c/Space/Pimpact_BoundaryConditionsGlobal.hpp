#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP
#define PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"




namespace Pimpact{



/// \brief global boundary conditions
///
/// \tparam dim
/// \ingroup SpaceObject
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

  Ti3 BCL_int_;
  Ti3 BCU_int_;

	BoundaryConditionsGlobal(
			EBCType BC1L=DirichletBC,
			EBCType BC1U=DirichletBC,
			EBCType BC2L=DirichletBC,
			EBCType BC2U=DirichletBC,
			EBCType BC3L=DirichletBC,
			EBCType BC3U=DirichletBC ) {

		BCL_int_[0] = static_cast<int>( BC1L );
		BCU_int_[0] = static_cast<int>( BC1U );
		BCL_int_[1] = static_cast<int>( BC2L );
		BCU_int_[1] = static_cast<int>( BC2U );
		BCL_int_[2] = static_cast<int>( BC3L );
		BCU_int_[2] = static_cast<int>( BC3U );

		if( 4==dim ) {
			BCL_int_[3] = static_cast<int>( PeriodicBC );
			BCU_int_[3] = static_cast<int>( PeriodicBC );
		}

	}

	BoundaryConditionsGlobal( TBC3 BCL_global, TBC3 BCU_global ) {

		for( int i=0; i<3; ++i ) {
			BCL_int_[i] = static_cast<int>( BCL_global[i] );
			BCU_int_[i] = static_cast<int>( BCU_global[i] );
		};
		if( 4==dim ) {
			BCL_int_[3] = static_cast<int>( PeriodicBC );
			BCU_int_[3] = static_cast<int>( PeriodicBC );
		}

	}

public:

	/// \name getter
	/// @{ 

  EBCType getBCL( const int& dir ) const { return( static_cast<EBCType>(BCL_int_[dir]) ); }
	EBCType getBCU( const int& dir ) const { return( static_cast<EBCType>(BCU_int_[dir]) ); }

  const int* getBCL() const { return( BCL_int_.getRawPtr() ); }
  const int* getBCU() const { return( BCU_int_.getRawPtr() ); }

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

	///  @} 
	
	/// \brief prints BC tuples
	///
	/// \param out output stream
  void print( std::ostream& out=std::cout ) const {
		out << "---BoundaryConditionsGlobal: ---\n";
		out << " BCL_global: " << BCL_int_ << "\n";
		out << " BCU_global: " << BCU_int_ << "\n";
  }

}; // end of class BoundaryConditionsGlobal




/// \brief creates global boundary conditions accordint to domaint Tyep
///
/// \tparam dim computational dimension
/// \param dtype domain type
///
/// \return 
/// \relates BoundaryConditionsGlobal
/// \relates EDomainType
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
