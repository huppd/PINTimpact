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
template< class Ordinal=int, int dim=3 >
class FieldSpace {

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;
  typedef const Teuchos::Tuple<Ordinal,3> TO3;
  //	using Teuchos::RCP;

  /// \brief constructor
  ///
  /// \param commf Fortran MPI communicator(an integer)
  /// \param comm MPI communicator
  /// \param nGlo amount of global gridpoints
  /// \param nLoc amount of local gridpoints
  /// \param bl lower bound of storage
  /// \param bu upper bound of storage
  FieldSpace(
      MPI_Fint commf,
      MPI_Comm comm,
      Ordinal dimension,
      TO nGlo,
      TO nLoc,
      TO3 bl,
      TO3 bu ):
        commf_(commf),
        comm_(comm),
        dim_(dimension),
        nGlo_(nGlo),
        nLoc_(nLoc),
        bl_(bl),
        bu_(bu)
  {};

  FieldSpace(const FieldSpace& fs):
    commf_(fs.commf_),
    comm_(fs.comm_),
    dim_(fs.dim_),
    nGlo_(fs.nGlo_),
    nLoc_(fs.nLoc_),
    bl_(fs.bl_),
    bu_(fs.bu_)
  {};

  /// prints to \c std::cout, only for debuging purpose
  void print( std::ostream& out=std::cout ) const {
    int rank;
    MPI_Comm_rank(comm_,&rank);
    if(0==rank) out << "\t---FieldSpace: ---\n";
    out << "rank: " << rank << " :commf: " << commf_ << "\n";
    out << "rank: " << rank << " :comput dim: " << dim << "\n";
    out << "rank: " << rank << " :pseudo dim: " << dim_ << "\n";
    out << "rank: " << rank << " :nGlo: " << nGlo_ << "\n";
    out << "rank: " << rank << " :nLoc: " << nLoc_ << "\n";
    out << "rank: " << rank << " :bl: " << bl_ << "\n";
    out << "rank: " << rank << " :bu: " << bu_ << "\n";
    MPI_Barrier(comm_);
  }


  const MPI_Fint commf_;
  MPI_Comm comm_;

  Ordinal dim_;

  TO nGlo_;
  TO nLoc_;
  TO3 bl_;
  TO3 bu_;


}; // end of class FieldSpace


extern "C" {
void SVS_get_comm(MPI_Fint&);
void FS_get_dim(int&);
void SVS_get_nGlo(int&,int&,int&);
void SVS_get_nLoc(int&,int&,int&);
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
  typedef typename FieldSpace<O,3>::TO3 TO3;

  MPI_Fint comm;
  SVS_get_comm( comm );

  int dim;
  FS_get_dim( dim );

  TO nGlo;
  SVS_get_nGlo( nGlo[0], nGlo[1], nGlo[2] );

  TO nLoc;
  SVS_get_nLoc( nLoc[0], nLoc[1], nLoc[2] );

  TO3 bl;
  SVS_get_bl( bl[0], bl[1], bl[2] );

  TO3 bu;
  SVS_get_bu( bu[0], bu[1], bu[2] );

  return(
      Teuchos::RCP<const FieldSpace<O,d> > (
          new FieldSpace<O,d>(
              comm,
              MPI_Comm_f2c(comm),
              dim,
              nGlo,
              nLoc,
              bl,
              bu ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_FIELDSPACE_HPP
