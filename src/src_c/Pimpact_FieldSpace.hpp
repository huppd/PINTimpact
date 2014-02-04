#pragma once
#ifndef PIMPACT_FIELDSPACE_HPP
#define PIMPACT_FIELDSPACE_HPP

#include "mpi.h"

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"


#include <iostream>

namespace Pimpact {

/**
 * public class, that stores neccessary information for Vectors(Fields), which are common for \c ScalarField and \c VectorField
 * \todo check constant/ make variables protected SV friend class
 * \todo remove indexes(add indexSpace to Vector)
 * \todo add procssor blocks
 * \todo add neibouring ranks
 */
template<class O>
class FieldSpace {


public:

	typedef O Ordinal ;
	typedef const Teuchos::Tuple<Ordinal,3> TO3;
//	using Teuchos::RCP;

	/** \brief constructor
	 *
	 * @param commf Fortran MPI communicator(an integer)
	 * @param comm MPI communicator
	 * @param nGlo amount of global gridpoints
	 * @param nLoc amount of local gridpoints
	 * @param sInd start index of gridpoints
	 * @param eInd last index of gridpoints
	 * @param bl lower bound of storage
	 * @param bu upper bound of storage
	 */
	FieldSpace(MPI_Fint commf, MPI_Comm comm, Ordinal dim, TO3 nGlo, TO3 nLoc, TO3 sInd, TO3 eInd, TO3 bl, TO3 bu):
		commf_(commf), comm_(comm), dim_(dim), nGlo_(nGlo), nLoc_(nLoc), sInd_(sInd), eInd_(eInd), bl_(bl), bu_(bu) {};

	FieldSpace(const FieldSpace& fs):
		commf_(fs.commf_), comm_(fs.comm_), dim_(fs.dim_), nGlo_(fs.nGlo_), nLoc_(fs.nLoc_), sInd_(fs.sInd_), eInd_(fs.eInd_), bl_(fs.bl_), bu_(fs.bu_) {};
	/** prints to \c std::cout, only for debuging purpose
	 *
	 */
	void print() const {
		int rank;
		MPI_Comm_rank(comm_,&rank);
		std::cout << "rank: " << rank << " :commf: " << commf_ << "\n";
		std::cout << "rank: " << rank << " :dim: " << dim_ << "\n";
		std::cout << "rank: " << rank << " :nGlo: " << nGlo_ << "\n";
		std::cout << "rank: " << rank << " :nLoc: " << nLoc_ << "\n";
		std::cout << "rank: " << rank << " :sInd: " << sInd_ << "\n";
		std::cout << "rank: " << rank << " :eInd: " << eInd_ << "\n";
		std::cout << "rank: " << rank << " :bl: " << bl_ << "\n";
		std::cout << "rank: " << rank << " :bu: " << bu_ << "\n";
	}


	const MPI_Fint commf_; // for savety resason this should be mpi_fint, but that is const int (not good)
	MPI_Comm comm_;

	Ordinal dim_;

	TO3 nGlo_;
	TO3 nLoc_;
	TO3 sInd_;
	TO3 eInd_;
	TO3 bl_;
	TO3 bu_;


}; // class FieldSpace

extern "C" {
	void SVS_get_comm(MPI_Fint&);
	void FS_get_dim(int&);
	void SVS_get_nGlo(int&,int&,int&);
	void SVS_get_nLoc(int&,int&,int&);
	void SVS_get_sInd(int&,int&,int&);
	void SVS_get_eInd(int&,int&,int&);
	void SVS_get_bl(int&,int&,int&);
	void SVS_get_bu(int&,int&,int&);
}

/** \brief function that creates \c Pimpact:FieldSpace
 * by getting values from \c IMPACT
 */
template<class Ordinal>
const Teuchos::RCP<const FieldSpace<Ordinal> > createFieldSpace(){

	typedef typename FieldSpace<Ordinal>::TO3 TO3;
	MPI_Fint comm;
	SVS_get_comm( comm );

	int dim;
	FS_get_dim( dim );

	TO3 nGlo;
	SVS_get_nGlo( nGlo[0], nGlo[1], nGlo[2] );

	TO3 nLoc;
	SVS_get_nLoc( nLoc[0], nLoc[1], nLoc[2] );

	TO3 sInd;
	SVS_get_sInd( sInd[0], sInd[1], sInd[2] );

	TO3 eInd;
	SVS_get_eInd( eInd[0], eInd[1], eInd[2] );

	TO3 bl;
	SVS_get_bl( bl[0], bl[1], bl[2] );

	TO3 bu;
	SVS_get_bu( bu[0], bu[1], bu[2] );

	return Teuchos::RCP<const FieldSpace<Ordinal> > (
			new FieldSpace<Ordinal>( comm, MPI_Comm_f2c(comm), dim, nGlo, nLoc, sInd, eInd, bl, bu ) );
}

} // namespace Pimpact

#endif // PIMPACT_FIELDSPACE_HPP
