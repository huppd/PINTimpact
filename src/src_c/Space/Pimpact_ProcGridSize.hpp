#pragma once
#ifndef PIMPACT_PROCGRIDSIZE_HPP
#define PIMPACT_PROCGRIDSIZE_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



//extern "C" {
//
//void fsetPGS(const int& np1, const int& np2, const int& np3 );
//
//}




namespace Pimpact{


/// \brief size of processor grid
/// \ingroup SpaceObject
/// \todo maybe inherit from Tuple
template<class Ordinal, int dim>
class ProcGridSize {

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;

  template< class OT, int dT>
  friend Teuchos::RCP<const ProcGridSize<OT,dT> > createProcGridSize( OT np1, OT np2, OT np3, OT npt );

  template< class OT, int dT>
  friend Teuchos::RCP<const ProcGridSize<OT,dT> > createProcGridSize( const Teuchos::Tuple<OT,dT>  );

  template< class OT, int dT>
  friend Teuchos::RCP<const ProcGridSize<OT,dT> > createProcGridSize( const Teuchos::Tuple<OT,dT>, const MPI_Comm&  );

  template< class OT, int dT>
  friend Teuchos::RCP<const ProcGridSize<OT,dT> > createProcGridSize( const Teuchos::Tuple<OT,dT>, const MPI_Comm&, bool test_yes  );
protected:

  TO procGridSize_;

  void test( const MPI_Comm& comm ) const {

    int commSize;
    MPI_Comm_size( comm, &commSize );

    int procSize = 1;

    for( int i=0; i<dim; ++i )
      procSize *= procGridSize_[i];

		TEUCHOS_TEST_FOR_EXCEPT( procSize != commSize );

    for( int i=0; i<dim; ++i )
			TEUCHOS_TEST_FOR_EXCEPT( procGridSize_[i]<1 );
  }


  ProcGridSize( TO procGridSize, const MPI_Comm& comm, bool test_yes ):
		procGridSize_( procGridSize ) {

			if( test_yes ) test( comm );

  };

public:

  const Ordinal& get( int i ) const {
    return( procGridSize_[i] );
  }
  Ordinal* get() const {
    return( procGridSize_.getRawPtr() );
  }

  TO getTuple() const {
    return( procGridSize_ );
  }

  void print( std::ostream& out=std::cout ) const {
		out << "---ProcGridSize: "<<procGridSize_ << " ---\n";
  };

};



/// \relates ProcGridSize
template<class O=int, int d=3>
Teuchos::RCP<const ProcGridSize<O,d> > createProcGridSize( O np1, O np2, O np3, O npt=0 ) {

  Teuchos::Tuple<O,d> temp;

	temp[0] = np1;
	temp[1] = np2;
	temp[2] = np3;
	if( 4==d ) temp[3] = npt;

  return(
      Teuchos::rcp(
          new ProcGridSize<O,d>( temp, MPI_COMM_WORLD, true ) ) );
}



/// \relates ProcGridSize
template<class O=int, int d=3>
Teuchos::RCP<const ProcGridSize<O,d> >
createProcGridSize( const Teuchos::Tuple<O,d> temp, const MPI_Comm& comm=MPI_COMM_WORLD, bool test_yes=true ){ 

  return(
      Teuchos::rcp(
          new ProcGridSize<O,d>( temp, comm, test_yes ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_PROCGRIDSIZE_HPP
