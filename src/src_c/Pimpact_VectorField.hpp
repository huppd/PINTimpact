#pragma once
#ifndef PIMPACT_VECTORFIELD_HPP
#define PIMPACT_VECTORFIELD_HPP

#include <vector>
#include <iostream>
#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"
//#include "Pimpact_Operator.hpp"



namespace Pimpact {

extern "C" {

  void SV_add(
			const int& N1,  const int& N2,  const int& N3,
			const int& SS1, const int& SS2, const int& SS3,
			const int& NN1, const int& NN2, const int& NN3,
			const int& b1L, const int& b2L, const int& b3L,
			const int& b1U, const int& b2U, const int& b3U,
  		double* phi, const double* const  phi1, const double* const  phi2,
  		const double& scalar1, const double& scalar2);

	void VF_compNorm(
			double* phi1, double* phi2, double*  phi3,
//			MPI_Fint& comm,
//			const int& N1,  const int& N2,  const int& N3,
//			const int& SS1, const int& SS2, const int& SS3,
//			const int& NN1, const int& NN2, const int& NN3,
//			const int& b1L, const int& b2L, const int& b3L,
//			const int& b1U, const int& b2U, const int& b3U,
//			double* phi,
			const bool& inf_yes, const bool& two_yes,
			double& normInf, double& normTwo);

  void VF_dot(
//  		MPI_Fint& comm,
//			const int& N1,  const int& N2,  const int& N3,
//			const int& SS1, const int& SS2, const int& SS3,
//			const int& NN1, const int& NN2, const int& NN3,
//			const int& b1L, const int& b2L, const int& b3L,
//			const int& b1U, const int& b2U, const int& b3U,
  		const double* const phi1U, const double* const phi1V, const double* const phi1W,
  		const double* const phi2U, const double* const phi2V, const double* const phi2W,
  		double& scalar);

  void SV_scale(
			const int& N1,  const int& N2,  const int& N3,
			const int& SS1, const int& SS2, const int& SS3,
			const int& NN1, const int& NN2, const int& NN3,
			const int& b1L, const int& b2L, const int& b3L,
			const int& b1U, const int& b2U, const int& b3U,
  		double* phi, const double& scalar );

  void SV_random(
			const int& N1,  const int& N2,  const int& N3,
			const int& SS1, const int& SS2, const int& SS3,
			const int& NN1, const int& NN2, const int& NN3,
			const int& b1L, const int& b2L, const int& b3L,
			const int& b1U, const int& b2U, const int& b3U,
  		double* phi );

  void SV_init(
			const int& N1,  const int& N2,  const int& N3,
			const int& SS1, const int& SS2, const int& SS3,
			const int& NN1, const int& NN2, const int& NN3,
			const int& b1L, const int& b2L, const int& b3L,
			const int& b1U, const int& b2U, const int& b3U,
  		double* phi, const double& scalar );

  void SV_print(
			const int& N1,  const int& N2,  const int& N3,
			const int& SS1, const int& SS2, const int& SS3,
			const int& NN1, const int& NN2, const int& NN3,
			const int& b1L, const int& b2L, const int& b3L,
			const int& b1U, const int& b2U, const int& b3U,
  		const double* phi );

  void VF_write( double* phiU, double* phiV, double* phiW, const int& count );

  void VF_init_Zero( double* phiU, double* phiV, double* phiW );
  void VF_init_2DPoiseuilleX( double* phiU, double* phiV, double* phiW );
  void VF_init_2DPoiseuilleY( double* phiU, double* phiV, double* phiW );
  void VF_init_2DPulsatileXC( double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );
  void VF_init_2DPulsatileYC( double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );
  void VF_init_2DPulsatileXS( double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );
  void VF_init_2DPulsatileYS( double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );
  void VF_init_Streaming( double* phiU, double* phiV, double* phiW );
}


/** \brief important basic Vector class
 * vector for a vector field, e.g.: velocity,
 * here also happens the fortran wrapping
 */
template<class S, class O>
class VectorField {

	template<class S1, class O1>
	friend class Grad;
	template<class S1,class O1>
	friend class Div;
	template<class S1,class O1>
	friend class Helmholtz;

public:
	typedef S Scalar;
	typedef O Ordinal;
//	friend class Operator;
//	friend class Grad<Scalar,Ordinal>;

private:

//	using Teuchos::RCP;
	typedef Scalar* ScalarArray;
	typedef VectorField<Scalar,Ordinal> MV;

public:
	typedef Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > >  IndexSpaces;

	VectorField(): fieldS_(Teuchos::null),innerIS_(Teuchos::null),fullIS_(Teuchos::null),vec_(0) {};

	VectorField(Teuchos::RCP<const FieldSpace<Ordinal> > fieldS, IndexSpaces innerIS, IndexSpaces fullIS ):fieldS_(fieldS),innerIS_(innerIS),fullIS_(fullIS) {
		Ordinal N = 1;
		for(int i=0; i<3; ++i)
			N *= nLoc(i)+bu(i)-bl(i);

		for(int i=0; i<3; ++i) {
			vec_[i] = new Scalar[N];
//#ifdef DEBUG
			for(int j=0; j<N; ++j){
				vec_[i][j] = 0.;
			}
//#endif
		}
	};

	/** copy constructor
	 * shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
	 * @param sF
	 * @param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	 */
	VectorField(const VectorField& vF, ECopyType copyType=ShallowCopy):fieldS_(vF.fieldS_),innerIS_(vF.innerIS_),fullIS_(vF.fullIS_) {
		Ordinal N = 1;
		for(int i=0; i<3; ++i)
			N *= nLoc(i)+bu(i)-bl(i);

		for(int i=0; i<3; ++i) {
			vec_[i] = new Scalar[N];

			switch( copyType ) {

				case ShallowCopy:
//	#ifdef DEBUG
					for(int j=0; j<N; ++j) {
						vec_[i][j] = 0;
					}
//	#endif
					break;
			case DeepCopy:
				for( int j=0; j<N; ++j) {
						vec_[i][j] = vF.vec_[i][j];
				}
				break;
			}
		}
	};

	~VectorField() { for(int i=0; i<3; ++i) delete[] vec_[i];}

  //! \name Attribute methods
  //@{

	Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const {return fieldS_;}


	/**
	 * \brief returns the length of Field.
	 * the vector length is withregard to the inner points such that
	 * \f[ N_u = (N_x-1)(N_y-2)(N_z-?) \]
	 * \f[ N_v = (N_x-2)(N_y-1)(N_z-?) \]
	 * \f[ N_w = (N_x-?)(N_y-?)(N_z-?) \]
	 * @return vect length \f[= N_u+N_v+N+w]f
	 */
	Ordinal getLength( bool dummy=false ) const {
  	Ordinal n = 0;
  	for( int i=0; i<dim(); ++i ) {
			Ordinal vl = 1;
			for( int j=0; j<dim(); ++j) {
				if( i==j )
					vl *= nGlo(j)-1;
				else
					vl *= nGlo(j)-2;
			}
			n += vl;
  	}
  	return( n );
  }


	/// \brief get number of stored Field's
	int getNumberVecs() const { return( 1 ); }


	//@}
	/// @name Update methods
	//@{

	/**
	 * \brief Replace \c this with \f$\alpha A + \beta B\f$.
	 */
	void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {
		// add test for consistent VectorSpaces in debug mode
		for( int i=0; i<dim(); ++i )
			SV_add(
						nLoc(0), nLoc(1), nLoc(2),
						sInd(0,i), sInd(1,i), sInd(2,i),
						eInd(0,i), eInd(1,i), eInd(2,i),
						bl(0),   bl(1),   bl(2),
						bu(0),   bu(1),   bu(2),
						vec_[i], A.vec_[i], B.vec_[i], alpha, beta);
	}


	/**
	 * \brief Scale each element of the vectors in \c this with \c alpha.
	 */
	void scale( const Scalar& alpha ) {
		for(int i=0; i<dim(); ++i)
			SV_scale(
					nLoc(0), nLoc(1), nLoc(2),
					sInd(0,i), sInd(1,i), sInd(2,i),
					eInd(0,i), eInd(1,i), eInd(2,i),
					bl(0),   bl(1),   bl(2),
					bu(0),   bu(1),   bu(2),
					vec_[i], alpha);
	}


	/**
	 * \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
	 */
	Scalar dot ( const MV& a ) const {
		/// \todo add test in debuging mode for testing equality of VectorSpaces
		Scalar b;
		VF_dot(
//				commf(),
//				nLoc(0), nLoc(1), nLoc(2),
//				sInd(0), sInd(1), sInd(2),
//				eInd(0), eInd(1), eInd(2),
//				bl(0),   bl(1),   bl(2),
//				bu(0),   bu(1),   bu(2),
				vec_[0],     vec_[1],   vec_[2],
				a.vec_[0], a.vec_[1], a.vec_[2],
				b);
		return( b );
	}


	//@}
  /// @name Norm method
  //@{

  /**
   * \brief Compute the norm of each individual vector.
   * Upon return, \c normvec[i] holds the value of \f$||this_i||_2^2\f$, the \c i-th column of \c this.
   * \attention the two norm is not the real two norm but its square
   * \todo implement OneNorm
   * \todo implement in fortran
   */
	Scalar norm(  Belos::NormType type = Belos::TwoNorm ) const {
		bool twoNorm_yes = false;
  	bool infNorm_yes = false;

  	switch(type) {
  	case Belos::TwoNorm: twoNorm_yes = true; break;
  	case Belos::InfNorm: infNorm_yes = true; break;
  	case Belos::OneNorm: std::cout << "norm: not implemented"; return 0.;
  	default: std::cout << "unkown norm"; return 0.;
  	}

  	Scalar normvec;
  	VF_compNorm(
  			vec_[0], vec_[1], vec_[2],
//  				commf(),
//  				nLoc(0), nLoc(1), nLoc(2),
//  				sInd(0), sInd(1), sInd(2),
//  				eInd(0), eInd(1), eInd(2),
//  				bl(0),   bl(1),   bl(2),
//  				bu(0),   bu(1),   bu(2),
//  				s_, // mabyeb .data() or getPTr()
  				infNorm_yes, twoNorm_yes,
  				normvec, normvec );
  	return( normvec );
  }


  //@}
  //! @name Initialization methods
  //@{


  /**
   * \brief mv := A
   * Assign (deep copy) A into mv.
   */
  void assign( const MV& a ) {
  	#ifdef DEBUG
  	for(int i=0; i<3; ++i) {
  		TEST_EQUALITY( nLoc(i), a.Nloc(i) )
			TEST_EQUALITY( bu(i), a.bu(i) )
			TEST_EQUALITY( bl(i), a.bl(i) )
  	}
  	#endif

  	Ordinal N = 1;
		for(int i=0; i<3; ++i)
			N *= nLoc(i)+bu(i)-bl(i);

		for( int d=0; d<3; ++d)
			for(int i=0; i<N; ++i) {
				vec_[d][i] = a.vec_[d][i];
			}
  }


  /**
   * \brief Replace the vectors with a random vectors.
   * depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is requiered
   */
  void random(bool useSeed = false, int seed = 1) {
  	for( int i=0; i<dim(); ++i )
			SV_random(
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0,i), sInd(1,i), sInd(2,i),
				eInd(0,i), eInd(1,i), eInd(2,i),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				vec_[i] );
  }

  /*! \brief Replace each element of the vector  with \c alpha.
   */
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
  	for( int i=0; i<dim(); ++i )
			SV_init(
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0,i), sInd(1,i), sInd(2,i),
				eInd(0,i), eInd(1,i), eInd(2,i),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				vec_[i], alpha);
  }

  /** \brief Replace each element of the vector \c vec[i] with \c alpha[i].
    */
   void init( const Teuchos::Tuple<Scalar,3>& alpha ) {
   	for( int i=0; i<dim(); ++i )
 			SV_init(
 				nLoc(0), nLoc(1), nLoc(2),
 				sInd(0,i), sInd(1,i), sInd(2,i),
 				eInd(0,i), eInd(1,i), eInd(2,i),
 				bl(0),   bl(1),   bl(2),
 				bu(0),   bu(1),   bu(2),
 				vec_[i], alpha[i]);
   }


  /**
   *  \brief initializes VectorField with the initial field defined in Fortran
   */
  void init_field( EFlowProfile flowType = Poiseuille2D_inX, double re=1., double om=1., double px = 1. ) {
  	switch(flowType) {
  				case ZeroProf :
						VF_init_Zero( vec_[0], vec_[1], vec_[2] );
						break;
  				case Poiseuille2D_inX :
						VF_init_2DPoiseuilleX( vec_[0], vec_[1], vec_[2] );
						break;
  				case Poiseuille2D_inY :
						VF_init_2DPoiseuilleY( vec_[0], vec_[1], vec_[2] );
						break;
  				case Pulsatile2D_inXC :
						VF_init_2DPulsatileXC( vec_[0], vec_[1], vec_[2], re, om, px);
						break;
  				case Pulsatile2D_inYC :
						VF_init_2DPulsatileYC( vec_[0], vec_[1], vec_[2], re, om, px);
						break;
  				case Pulsatile2D_inXS :
						VF_init_2DPulsatileXS( vec_[0], vec_[1], vec_[2], re, om, px);
						break;
  				case Pulsatile2D_inYS:
						VF_init_2DPulsatileYS( vec_[0], vec_[1], vec_[2], re, om, px);
						break;
  				case Streaming2D:
						VF_init_Streaming( vec_[0], vec_[1], vec_[2] );
						break;

  	  	}
  }


  //@}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os )  {
  		int rank;
			MPI_Comm_rank(comm(),&rank);
			for(int i=0; i<dim(); ++i) {
				os << "rank: " << rank << " :dir: " << i << "\n";
				os << "rank: " << rank << " :nGlo: " << nGlo(i) << "\n";
				os << "rank: " << rank << " :nLoc: " << nLoc(i) << "\n";
				for( int j=0; j<3; ++j ) {
					os << "rank: " << rank << "field: " << j << " :sInd: " << sInd(i,j) << "\n";
					os << "rank: " << rank << "field: " << j << " :eInd: " << eInd(i,j) << "\n";
				}
				os << "rank: " << rank << " :bl: " << bl(i) << "\n";
				os << "rank: " << rank << " :bu: " << bu(i) << "\n\n";
			}
//		Ordinal N = 1;
//		for(int i=0; i<3; ++i)
//			N *= nLoc(i)+bu(i)-bl(i);
//		Ordinal Nx = nLoc(0)+bu(0)-bl(0);
//		Ordinal Ny = nLoc(1)+bu(1)-bl(1);
//		Ordinal Nz = nLoc(2)+bu(2)-bl(2);

//		Scalar bla = s_[0];
//					std::cout << "s_[0]: "<< bla <<"\n" ;
//		for(int ix=0; ix<Nx; ++ix) {
//			for(int iy=0; iy<Ny; ++iy) {
//				for(int iz=0; iz<Nz;++iz) {
//					std::cout << "rank: " << rank << " " <<
//							"ind: (" << ix+bl(0) << ", " << iy+bl(1) << ", " << iz+bl(2) <<  ") u(ind): " << s_[iz + Nz*iy + Nz*Ny*ix] << "\n" ;
//				}
//				std::cout << "\n";
//			}
//			std::cout << "\n";
//		}
		std::cout << "rank: " << rank << "\n";
		for( int i=0; i<dim(); ++i ) {
			std::cout << "field: " << i << "\n";
			SV_print(
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0,i), sInd(1,i), sInd(2,i),
				eInd(0,i), eInd(1,i), eInd(2,i),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				vec_[i] );
		}

  }

  void write( int count=0 ) {
  	VF_write( vec_[0], vec_[1], vec_[2], count );
  }

protected:
	Teuchos::RCP<const FieldSpace<Ordinal> > fieldS_;
	IndexSpaces innerIS_;
	IndexSpaces fullIS_;
	Teuchos::Tuple<ScalarArray,3> vec_;
//	ScalarArray vec_[3];

	/**
	 * \todo add good documetnation here
	 * @return
	 */
	MPI_Fint& commf() const { return  fieldS_->commf_ ; }
	MPI_Comm comm() const { return  fieldS_->comm_ ; }
	const int& dim() const { return fieldS_->dim_; }
	const Ordinal& nGlo(int i) const { return fieldS_->nGlo_[i]; }
	const Ordinal& nLoc(int i) const { return  fieldS_->nLoc_[i]; }
	const Ordinal& sInd(int i, int fieldType) const { return innerIS_[fieldType]->sInd_[i]; }
	const Ordinal& eInd(int i, int fieldType) const { return innerIS_[fieldType]->eInd_[i]; }
	const Ordinal& sIndB(int i, int fieldType) const { return fullIS_[fieldType]->sInd_[i]; }
	const Ordinal& eIndB(int i, int fieldType) const { return fullIS_[fieldType]->eInd_[i]; }
	const Ordinal& bl(int i) const { return fieldS_->bl_[i]; }
	const Ordinal& bu(int i) const { return fieldS_->bu_[i]; }

}; //class VectorField


/** @brief creates a vector field(vector) belonging to a \c FieldSpace and two \c IndexSpaces
 *
 * @param sVS scalar Vector Space to which returned vector belongs
 * @return field vector
 */
template<class Scalar, class Ordinal>
Teuchos::RCP< VectorField<Scalar,Ordinal> > createVectorField(
		Teuchos::RCP<const FieldSpace<Ordinal> > fieldS,
		typename VectorField<Scalar,Ordinal>::IndexSpaces innerIS,
		typename VectorField<Scalar,Ordinal>::IndexSpaces fullIS ) {
	return Teuchos::RCP<VectorField<Scalar,Ordinal> > (
				new VectorField<Scalar,Ordinal>( fieldS, innerIS, fullIS ) );
}

} // namespace Pimpact

#endif // PIMPACT_VECTORFIELD_HPP
