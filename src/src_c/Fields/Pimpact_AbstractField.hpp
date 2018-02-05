#pragma once
#ifndef PIMPACT_ABSTRACTFIELD_HPP
#define PIMPACT_ABSTRACTFIELD_HPP


#include "mpi.h"

#include "Teuchos_RCP.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Grid.hpp"




namespace Pimpact {


template<class GridT>
class AbstractField {

public:


  static const int dimension = GridT::dimension;

  AbstractField(const Teuchos::RCP<const GridT>& grid):grid_(grid) {};

protected:

  using ST = typename GridT::Scalar;

  Teuchos::RCP<const GridT> grid_;

  /// \note openMPI is picky her with const references :(
  /// \note think about making pubplic
  constexpr ST reduce(const MPI_Comm& comm, ST normLocal, const MPI_Op& op=MPI_SUM)  {

    ST normGlob;
    MPI_Allreduce(&normLocal, &normGlob, 1, MPI_REAL8, op, comm);
    return normGlob;
  }

}; // end of class AbstractField


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_ABSTRACTFIELD_HPP
