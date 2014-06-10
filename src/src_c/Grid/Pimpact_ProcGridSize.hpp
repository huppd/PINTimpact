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


template<class Ordinal,int dim=3>
class ProcGridSize {

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;

protected:

  TO procGridSize_;

public:

  ProcGridSize( Ordinal np1, Ordinal np2, Ordinal np3 ):
    procGridSize_( Teuchos::tuple(np1, np2, np3) ) {
    set_Impact();
  };


  ProcGridSize( TO procGridSize ):
    procGridSize_( procGridSize ) {
    set_Impact();
  };

  void set_Impact(){
    fsetPGS( procGridSize_[0], procGridSize_[1], procGridSize_[2] );
  };

  Ordinal* getRawPtr() {
    return( procGridSize_.getRawPtr() );
  }

  void print( std::ostream& out ) {
    for( int i=0; i<dim; ++i)
      out << "\t#proc"<<i<<": "<<procGridSize_[i];
    out<< "\n";
  };
};


/// \relates ProcGridSize
template<class Ordinal>
Teuchos::RCP<ProcGridSize<Ordinal,3> > createProcGridSize( Ordinal np1=2, Ordinal np2=2, Ordinal np3=1 ) {
  return(
      Teuchos::rcp(
          new ProcGridSize<Ordinal,3>( np1, np2, np3 ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_PROCGRIDSIZE_HPP
