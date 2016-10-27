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

	using TBC3 = const Teuchos::Tuple<EBCType,dim>;
	using Ti3 = const Teuchos::Tuple<int,dim>;

	template<int d>
  friend Teuchos::RCP<const BoundaryConditionsGlobal<d> > createBoudaryConditionsGlobal();

	template<int d>
  friend Teuchos::RCP<const BoundaryConditionsGlobal<d> > createBoudaryConditionsGlobal( EDomainType dtype );

template<int d>
friend Teuchos::RCP<const BoundaryConditionsGlobal<d> >
createBoudaryConditionsGlobal( const Teuchos::RCP<Teuchos::ParameterList>& pl );


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

  constexpr const int& getBCL( const int& dir ) const { return( BCL_int_[dir] ); }
	constexpr const int& getBCU( const int& dir ) const { return( BCU_int_[dir] ); }

  constexpr const int* getBCL() const { return( BCL_int_.getRawPtr() ); }
  constexpr const int* getBCU() const { return( BCU_int_.getRawPtr() ); }

	constexpr const Teuchos::Tuple<int,dim> periodic() const {

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




template<int dim=3>
Teuchos::RCP<const BoundaryConditionsGlobal<dim> >
createBoudaryConditionsGlobal( const Teuchos::RCP<Teuchos::ParameterList>& pl ) {

	return( Teuchos::rcp( new BoundaryConditionsGlobal<dim>(
					pl->get<EBCType>( "lower X", DirichletBC ),
					pl->get<EBCType>( "upper X", DirichletBC ),
					pl->get<EBCType>( "lower Y", DirichletBC ),
					pl->get<EBCType>( "upper Y", DirichletBC ),
					pl->get<EBCType>( "lower Z", DirichletBC ),
					pl->get<EBCType>( "upper Z", DirichletBC ) ) ) );
}




} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP
