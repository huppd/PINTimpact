#pragma once
#ifndef PIMPACT_PROCGRIDSIZE_HPP
#define PIMPACT_PROCGRIDSIZE_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



extern "C" {

void fsetPGS(const int& np1, const int& np2, const int& np3 );

}



namespace Pimpact{


template<class Ordinal=int,int dim=3>
class ProcGridSize {

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;

protected:

  TO procGridSize_;

  void test() {

    int commSize;
    MPI_Comm_size( MPI_COMM_WORLD, &commSize );

    int procSize = 1;

    for( int i=0; i<dim; ++i )
      procSize *= procGridSize_[i];

    if( procSize != commSize)
      std::cout << "!!!ERROR! ProcGridSize: ( procSize("<<procSize<<") != commSize("<<commSize<<")  differs from number of allocated processors !!!";
      for( int i=0; i<dim; ++i )
        if( procGridSize_[i]<1 )
          std::cout << "!!!ERROR! ProcGridSize: procGridSize_["<<i<<"]="<<procGridSize_[i]<<"<1 !!!\n";
  }

public:

  ProcGridSize( Ordinal np1, Ordinal np2, Ordinal np3 ):
    procGridSize_( Teuchos::tuple(np1, np2, np3) ) {

    test();
    set_Impact();
  };

  ProcGridSize( Ordinal np1, Ordinal np2, Ordinal np3, Ordinal npt ):
    procGridSize_( Teuchos::tuple(np1, np2, np3, npt) ) {

    test();
    set_Impact();
  };


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
Teuchos::RCP<ProcGridSize<Ordinal,3> > createProcGridSize( Ordinal np1, Ordinal np2, Ordinal np3 ) {
  return(
      Teuchos::rcp(
          new ProcGridSize<Ordinal,3>( np1, np2, np3 ) ) );
}

/// \relates ProcGridSize
template<class Ordinal=int>
Teuchos::RCP<ProcGridSize<Ordinal,4> > createProcGridSize( Ordinal np1, Ordinal np2, Ordinal np3, Ordinal npt ) {
  return(
      Teuchos::rcp(
          new ProcGridSize<Ordinal,4>( np1, np2, np3, npt ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_PROCGRIDSIZE_HPP
