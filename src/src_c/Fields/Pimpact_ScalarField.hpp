#pragma once
#ifndef PIMPACT_SCALARFIELD_HPP
#define PIMPACT_SCALARFIELD_HPP

#include <vector>
#include <iostream>
#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Pimpact_Types.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"



namespace Pimpact {

extern "C" {


void SF_add(
	const int& N1,  const int& N2,  const int& N3,
	const int& SS1, const int& SS2, const int& SS3,
	const int& NN1, const int& NN2, const int& NN3,
	const int& b1L, const int& b2L, const int& b3L,
	const int& b1U, const int& b2U, const int& b3U,
	double* phi, const double* const  phi1, const double* const  phi2,
	const double& scalar1, const double& scalar2);


void SF_abs(
    const int& N1,  const int& N2,  const int& N3,
    const int& SS1, const int& SS2, const int& SS3,
    const int& NN1, const int& NN2, const int& NN3,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phi, const double* const  phi1 );


void SF_reciprocal(
    const int& N1,  const int& N2,  const int& N3,
    const int& SS1, const int& SS2, const int& SS3,
    const int& NN1, const int& NN2, const int& NN3,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phi, const double* const  phi1 );


void SF_compNorm(const MPI_Fint& comm,
	const int& N1,  const int& N2,  const int& N3,
	const int& SS1, const int& SS2, const int& SS3,
	const int& NN1, const int& NN2, const int& NN3,
	const int& b1L, const int& b2L, const int& b3L,
	const int& b1U, const int& b2U, const int& b3U,
	double* phi,
	const bool& inf_yes, const bool& two_yes,
	double& normInf, double& normTwo);


void SF_weightedNorm(
  const MPI_Fint& comm,
  const int& N1,  const int& N2,  const int& N3,
  const int& SS1, const int& SS2, const int& SS3,
  const int& NN1, const int& NN2, const int& NN3,
  const int& b1L, const int& b2L, const int& b3L,
  const int& b1U, const int& b2U, const int& b3U,
  double* phi,
  const double* const weights,
  double& norm);


void SF_dot( const MPI_Fint& comm,
	const int& N1,  const int& N2,  const int& N3,
	const int& SS1, const int& SS2, const int& SS3,
	const int& NN1, const int& NN2, const int& NN3,
	const int& b1L, const int& b2L, const int& b3L,
	const int& b1U, const int& b2U, const int& b3U,
	const double* const phi1, const double* const phi2, double& scalar);


void SF_scale(
	const int& N1,  const int& N2,  const int& N3,
	const int& SS1, const int& SS2, const int& SS3,
	const int& NN1, const int& NN2, const int& NN3,
	const int& b1L, const int& b2L, const int& b3L,
	const int& b1U, const int& b2U, const int& b3U,
	double* phi, const double& scalar );


void SF_scale2(
    const int& N1,  const int& N2,  const int& N3,
    const int& SS1, const int& SS2, const int& SS3,
    const int& NN1, const int& NN2, const int& NN3,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phi, const double* const  phi1 );


void SF_random(
	const int& N1,  const int& N2,  const int& N3,
	const int& SS1, const int& SS2, const int& SS3,
	const int& NN1, const int& NN2, const int& NN3,
	const int& b1L, const int& b2L, const int& b3L,
	const int& b1U, const int& b2U, const int& b3U,
	double* phi );


void SF_init(
	const int& N1,  const int& N2,  const int& N3,
	const int& SS1, const int& SS2, const int& SS3,
	const int& NN1, const int& NN2, const int& NN3,
	const int& b1L, const int& b2L, const int& b3L,
	const int& b1U, const int& b2U, const int& b3U,
	double* phi, const double& scalar );


void SF_print(
	const int& N1,  const int& N2,  const int& N3,
	const int& SS1, const int& SS2, const int& SS3,
	const int& NN1, const int& NN2, const int& NN3,
	const int& b1L, const int& b2L, const int& b3L,
	const int& b1U, const int& b2U, const int& b3U,
	const double* phi );


void SF_write( double* phi, const int& count );


} // end of extern 'C'


/// \brief important basic Vector class
/// vector for a scalar field, e.g.: pressure,
/// \note all indexing is done in Fortran
/// \ingroup Field
template<class S, class O>
class ScalarField {


	template<class S1,class O1>
	friend class Grad;
	template<class S1,class O1>
	friend class Div;
	template<class S1,class O1>
	friend class Div_Grad;


public:
	typedef S Scalar;
	typedef O Ordinal;


protected:
//	using Teuchos::RCP;
	typedef Scalar* array;
	typedef ScalarField<Scalar,Ordinal> MV;


public:
	ScalarField():s_(0),fieldSpace_(Teuchos::null) {};

	ScalarField( const Teuchos::RCP<const FieldSpace<Ordinal> >& sVS):fieldSpace_(sVS) {
		Ordinal N = 1;
		for(int i=0; i<3; ++i)
			N *= nLoc(i)+bu(i)-bl(i);

		s_ = new Scalar[N];
//#ifdef DEBUG
		for(int i=0; i<N; ++i){
			s_[i] = 0.;
		}
//#endif
	};


	/**
	 * \brief copy constructor.
	 *
	 * shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
	 * @param sF
	 * @param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	 */
	ScalarField(const ScalarField& sF, ECopyType copyType=ShallowCopy):fieldSpace_(sF.fieldSpace_) {
		Ordinal N = 1;
		for(int i=0; i<3; ++i)
			N *= nLoc(i)+bu(i)-bl(i);

		s_ = new Scalar[N];

		switch( copyType ) {

			case ShallowCopy:
//#ifdef DEBUG
				for(int i=0; i<N; ++i) {
					s_[i] = 0;
				}
//#endif
				break;
		case DeepCopy:
			for( int i=0; i<N; ++i) {
					s_[i] = sF.s_[i];
			}
			break;
		}
	};


	~ScalarField() { delete[] s_;}

	Teuchos::RCP<MV> clone( ECopyType ctype=DeepCopy ) const {
	  return( Teuchos::rcp( new MV(*this, ctype) ) );
	}

	/// \name Attribute methods
	///@{

	Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const { return( fieldSpace_ ); }


	/// \brief returns the length of Field.
	Ordinal getLength( bool dummy=false ) const {
		Ordinal vl = 1;
		for(int i = 0; i<dim(); ++i)
			vl *= nGlo(i);
		return( vl );
	}


	/// \brief get number of stored Field's
	int getNumberVecs() const { return( 1 ); }


	//@}
	/// @name Update methods
	//@{

	/// \brief Replace \c this with \f$\alpha A + \beta B\f$.
	void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {
		// add test for consistent VectorSpaces in debug mode
		SF_add(
						nLoc(0), nLoc(1), nLoc(2),
						sInd(0), sInd(1), sInd(2),
						eInd(0), eInd(1), eInd(2),
						bl(0),   bl(1),   bl(2),
						bu(0),   bu(1),   bu(2),
						s_, A.s_, B.s_, alpha, beta);
	}


  /**
   * \brief Put element-wise absolute values of source vector \c y into this
   * vector.
   *
   * Here x represents this vector, and we update it as
   * \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
   * \return Reference to this object
   * \todo implement me
   */
  void abs(const MV& y) {
		// add test for consistent VectorSpaces in debug mode
		SF_abs(
						nLoc(0), nLoc(1), nLoc(2),
						sInd(0), sInd(1), sInd(2),
						eInd(0), eInd(1), eInd(2),
						bl(0),   bl(1),   bl(2),
						bu(0),   bu(1),   bu(2),
						s_, y.s_ );
  }


  /**
    * \brief Put element-wise reciprocal of source vector \c y into this vector.
    *
    * Here x represents this vector, and we update it as
    * \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
    * \return Reference to this object
    * \todo implement me
    */
   void reciprocal(const MV& y){
     // add test for consistent VectorSpaces in debug mode
     SF_reciprocal(
         nLoc(0), nLoc(1), nLoc(2),
         sInd(0), sInd(1), sInd(2),
         eInd(0), eInd(1), eInd(2),
         bl(0),   bl(1),   bl(2),
         bu(0),   bu(1),   bu(2),
         s_, y.s_ );
   }


	/**
	 * \brief Scale each element of the vector with \c alpha.
	 */
	void scale( const Scalar& alpha ) {
		SF_scale(
					nLoc(0), nLoc(1), nLoc(2),
					sInd(0), sInd(1), sInd(2),
					eInd(0), eInd(1), eInd(2),
					bl(0),   bl(1),   bl(2),
					bu(0),   bu(1),   bu(2),
					s_, alpha);
	}


  /**
   * \brief Scale this vector <em>element-by-element</em> by the vector a.
   *
   * Here x represents this vector, and we update it as
   * \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
   * \return Reference to this object
   * \todo implement me
   */
  void scale(const MV& a) {
    // add test for consistent VectorSpaces in debug mode
    SF_scale2(
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
        s_, a.s_ );
  }


	/**
	 * \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
	*/
	Scalar dot ( const MV& a ) const {
		/// \todo add test in debuging mode for testing equality of VectorSpaces
		Scalar b;
		SF_dot( commf(),
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0), sInd(1), sInd(2),
				eInd(0), eInd(1), eInd(2),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				s_, a.s_, b);
		return b;
	}


	///@}
  /// @name Norm method
  ///@{

	/**
	 * \brief compute the norm
	 * \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
	 * \todo implement OneNorm
	 */
	Scalar norm(  Belos::NormType type = Belos::TwoNorm ) const {
		bool twoNorm_yes = false;
		bool infNorm_yes = false;

		switch(type) {
		case Belos::TwoNorm: twoNorm_yes = true; break;
		case Belos::InfNorm: infNorm_yes = true; break;
		case Belos::OneNorm: std::cout << "!!! Warning Belos::OneNorm not implemented \n"; return(0.);
  	default: std::cout << "!!! Warning unknown Belos::NormType:\t" << type << "\n"; return(0.);
		}

		Scalar normvec;
		SF_compNorm(
				commf(),
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0), sInd(1), sInd(2),
				eInd(0), eInd(1), eInd(2),
				  bl(0),   bl(1),   bl(2),
				  bu(0),   bu(1),   bu(2),
				s_,
				infNorm_yes, twoNorm_yes,
				normvec, normvec );
		return( normvec );
	}


  /**
   * \brief Weighted 2-Norm.
   *
   * Here x represents this vector, and we compute its weighted norm as follows:
   * \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
   * \return \f$ \|x\|_w \f$
   */
  double norm(const MV& weights) const {
    Scalar normvec;
    SF_weightedNorm(
        commf(),
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
        s_, weights.s_,
        normvec );
    return( normvec );
  }


  //@}
  //! @name Initialization methods
  //@{

	/** \brief mv := A
	 * Assign (deep copy) \c a into \c this.
	 * total deep, boundaries and everythin.
	 * \note the \c FieldSpace is not take care of assuming every field is generated with one
	 * \note "indexing" is done c++
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

  		for(int i=0; i<N; ++i) {
  			s_[i] = a.s_[i];
  		}
  }


	/**
   * \brief Replace the vectors with a random vectors.
   * depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is requiered
   */
	void random(bool useSeed = false, int seed = 1) {
		SF_random(
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0), sInd(1), sInd(2),
				eInd(0), eInd(1), eInd(2),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				s_);
	}


  /*! \brief Replace each element of the vector  with \c alpha.
   */
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
  	SF_init(
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0), sInd(1), sInd(2),
				eInd(0), eInd(1), eInd(2),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				s_, alpha);
  }

  //@}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os )  {
  		int rank;
			MPI_Comm_rank(comm(),&rank);
			for(int i=0; i<3; ++i) {
				os << "rank: " << rank << " :dir: " << i << "\n";
				os << "rank: " << rank << " :nGlo: " << nGlo(i) << "\n";
				os << "rank: " << rank << " :nLoc: " << nLoc(i) << "\n";
				os << "rank: " << rank << " :sInd: " << sInd(i) << "\n";
				os << "rank: " << rank << " :eInd: " << eInd(i) << "\n";
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
  	SF_print(
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0), sInd(1), sInd(2),
				eInd(0), eInd(1), eInd(2),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				s_ );

  }


  void write( int count=0 ) {
  	SF_write( s_, count );
  }

protected:
	const Teuchos::RCP<const FieldSpace<Ordinal> > fieldSpace_;
	array s_;

	const MPI_Fint& commf() const { return(  fieldSpace_->commf_ ); }
	const MPI_Comm& comm() const { return( fieldSpace_->comm_ ); }
	const int& dim() const {return( fieldSpace_->dim_ ); }
	const Ordinal& nGlo(int i) const { return( fieldSpace_->nGlo_[i] ); }
	const Ordinal& nLoc(int i) const { return(  fieldSpace_->nLoc_[i] ); }
	const Ordinal& sInd(int i) const { return( fieldSpace_->sInd_[i] ); }
	const Ordinal& eInd(int i) const { return( fieldSpace_->eInd_[i] ); }
	const Ordinal& bl(int i) const { return( fieldSpace_->bl_[i] ); }
	const Ordinal& bu(int i) const { return( fieldSpace_->bu_[i] ); }

}; // end of class ScalarField



/// \brief creates a scalar field(vector) belonging to a FieldSpace
///
/// \param fS scalar Vector Space to which returned vector belongs
/// \return scalar vector
/// \relates ScalarField
template<class Scalar, class Ordinal>
Teuchos::RCP< ScalarField<Scalar,Ordinal> > createScalarField( const Teuchos::RCP<const FieldSpace<Ordinal> >& fS) {
	return( Teuchos::RCP<ScalarField<Scalar,Ordinal> > (
			new ScalarField<Scalar,Ordinal>( fS ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_SCALARFIELD_HPP
