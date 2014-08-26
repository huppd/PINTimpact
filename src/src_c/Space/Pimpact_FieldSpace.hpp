#pragma once
#ifndef PIMPACT_FIELDSPACE_HPP
#define PIMPACT_FIELDSPACE_HPP

#include "mpi.h"

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"

#include <iostream>



namespace Pimpact {


/// public class, that stores neccessary information for Vectors(Fields), which are common for \c ScalarField and \c VectorField
/// \todo check constant/ make variables protected SV friend class
/// \todo remove indexes(add indexSpace to ScalarField) sInd/eInd
template< class Ordinal=int, int dim=3>
class FieldSpace {

  template< class OT, int dT >
  friend const Teuchos::RCP<const FieldSpace<OT,dT> > createFieldSpace();

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;

  Ordinal dim_;

  TO bl_;
  TO bu_;

  Teuchos::Tuple<Ordinal,3> ls_; /// < \brief Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup)

protected:

  /// \brief constructor
  ///
  /// \param bl lower bound of storage
  /// \param bu upper bound of storage
  FieldSpace(
      Ordinal dimension,
      TO bl,
      TO bu ):
        dim_(dimension),
        bl_(bl),
        bu_(bu),
        ls_(Teuchos::tuple(-1,-1,-1))
  {};

  FieldSpace(
      const FieldSpace& fs ):
        dim_(fs.dim_),
        bl_(fs.bl_),
        bu_(fs.bu_),
        ls_(Teuchos::tuple(-1,-1,-1))
  {};

public:

  /// prints to \c std::cout, only for debuging purpose
  void print( std::ostream& out=std::cout ) const {
    out << "\t---FieldSpace: ---\n";
    out << "comput dim: " << dim << "\n";
    out << "pseudo dim: " << dim_ << "\n";
    out << "bl: " << bl_ << "\n";
    out << "bu: " << bu_ << "\n";
    out << "ls: " << ls_ << "\n";
  }


}; // end of class FieldSpace


extern "C" {
void FS_get_dim(int&);
void SVS_get_bl(int&,int&,int&);
void SVS_get_bu(int&,int&,int&);
}


/// \brief function that creates \c Pimpact:FieldSpace
/// by getting values from \c IMPACT
/// should be changed by getting vallues from \c ProcGridSize and \c GridSize
/// \relates FieldSpace
template< class O=int, int d=3 >
const Teuchos::RCP<const FieldSpace<O,d> > createFieldSpace(){

  typedef typename FieldSpace<O,d>::TO TO;

  int dim;
  FS_get_dim( dim );

  TO bl;
  SVS_get_bl( bl[0], bl[1], bl[2] );
  if( d==4 ) bl[3] = -1;

  TO bu;
  SVS_get_bu( bu[0], bu[1], bu[2] );
  if( d==4 ) bu[3] = 0;

  return(
      Teuchos::RCP<const FieldSpace<O,d> > (
          new FieldSpace<O,d>(
              dim,
              bl,
              bu ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_FIELDSPACE_HPP
