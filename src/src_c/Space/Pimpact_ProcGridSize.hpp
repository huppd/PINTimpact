#pragma once
#ifndef PIMPACT_PROCGRIDSIZE_HPP
#define PIMPACT_PROCGRIDSIZE_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



extern "C" {
void fgetPGS(      int& np1,       int& np2,       int& np3 );
void fsetPGS(const int& np1, const int& np2, const int& np3 );
}



namespace Pimpact{



/// \brief size of processor grid
/// \ingroup Space
template<class Ordinal=int,int dim=3>
class ProcGridSize {

public:


  typedef const Teuchos::Tuple<Ordinal,dim> TO;


  template<class OT>
  friend Teuchos::RCP<ProcGridSize<OT,3> > createProcGridSize();

  template< class OT, int dT >
  friend Teuchos::RCP<ProcGridSize<OT,dT> > createProcGridSize( OT np1, OT np2, OT np3, OT npt=0 );


protected:

  TO procGridSize_;

  void test() {

    int commSize;
    MPI_Comm_size( MPI_COMM_WORLD, &commSize );

    int procSize = 1;

    for( int i=0; i<dim; ++i )
      procSize *= procGridSize_[i];

    TEUCHOS_TEST_FOR_EXCEPTION(
        procSize != commSize,
        std::logic_error,
        "!!!ERROR! ProcGridSize:  differs from number of allocated processors !!!" );
    for( int i=0; i<dim; ++i )
      TEUCHOS_TEST_FOR_EXCEPTION(
          procGridSize_[i]<1,
        std::logic_error,
        "!!!ERROR! ProcGridSize: has to be greater than one!!!" );
  }


  ProcGridSize( TO procGridSize ):
    procGridSize_( procGridSize ) {

    test();

  };

public:

  const Ordinal& get( int i ) {
    return( procGridSize_[i] );
  }

  void set_Impact(){
    fsetPGS( procGridSize_[0], procGridSize_[1], procGridSize_[2] );
  };

  Ordinal* getRawPtr() {
    return( procGridSize_.getRawPtr() );
  }

  void print( std::ostream& out=std::cout ) {
    for( int i=0; i<dim; ++i)
      out << "\t#proc"<<i<<": "<<procGridSize_[i];
    out<< "\n";
  };
};



/// \relates ProcGridSize
template<class Ordinal=int>
Teuchos::RCP<ProcGridSize<Ordinal,3> > createProcGridSize() {
  typedef const Teuchos::Tuple<Ordinal,3> TO;

  TO procGridSize;

  fgetPGS( procGridSize[0], procGridSize[1], procGridSize[2] );

  return(
      Teuchos::rcp(
          new ProcGridSize<Ordinal,3>( procGridSize ) ) );
}



/// \relates ProcGridSize
template< class O=int, int d=3 >
Teuchos::RCP<ProcGridSize<O,d> > createProcGridSize( O np1, O np2, O np3, O npt=0 ) {

  Teuchos::Tuple<O,d> temp;
  if( 3==d ) {
    temp[0] = np1;
    temp[1] = np2;
    temp[2] = np3;
  }
  else if( 4==d ) {
    temp[0] = np1;
    temp[1] = np2;
    temp[2] = np3;
    temp[3] = npt;
  }
  return(
      Teuchos::rcp(
          new ProcGridSize<O,d>( temp ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_PROCGRIDSIZE_HPP
