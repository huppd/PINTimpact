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


void VF_compNorm(
    const MPI_Fint& comm, const int& dimens,
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
		const int& b1U, const int& b2U, const int& b3U,
    const double* const phi1, const double* const phi2, const double* const phi3,
    const bool& inf_yes, const bool& two_yes,
    double& normInf, double& normTwo );


void VF_weightedNorm(
    const MPI_Fint& comm, const int& dimens,
		const int& N1,  const int& N2,  const int& N3,
		const int& S1U, const int& S2U, const int& S3U,
		const int& N1U, const int& N2U, const int& N3U,
		const int& S1V, const int& S2V, const int& S3V,
		const int& N1V, const int& N2V, const int& N3V,
		const int& S1W, const int& S2W, const int& S3W,
		const int& N1W, const int& N2W, const int& N3W,
		const int& b1L, const int& b2L, const int& b3L,
		const int& b1U, const int& b2U, const int& b3U,
    const double* const phi1, const double* const phi2, const double* const phi3,
    const double* const wei1, const double* const wei2, const double* const wei3,
    double& norm );


void VF_dot(
  		const MPI_Fint& comm,
  		const int& dimens,
			const int& N1,  const int& N2,  const int& N3,
			const int& S1U, const int& S2U, const int& S3U,
			const int& N1U, const int& N2U, const int& N3U,
			const int& S1V, const int& S2V, const int& S3V,
			const int& N1V, const int& N2V, const int& N3V,
			const int& S1W, const int& S2W, const int& S3W,
			const int& N1W, const int& N2W, const int& N3W,
			const int& b1L, const int& b2L, const int& b3L,
			const int& b1U, const int& b2U, const int& b3U,
			const double* const phi1U, const double* const phi1V, const double* const phi1W,
			const double* const phi2U, const double* const phi2V, const double* const phi2W,
			double& scalar);


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


void VF_write( double* phiU, double* phiV, double* phiW, const int& count );


void VF_init_Zero(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_2DPoiseuilleX(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
		const int& b1L, const int& b2L, const int& b3L,
		const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_2DPoiseuilleY(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_2DPulsatileXC(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );


void VF_init_2DPulsatileYC(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );


void VF_init_2DPulsatileXS(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );


void VF_init_2DPulsatileYS(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );


void VF_init_Streaming(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_Circle(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_RankineVortex(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );

}


/// \brief important basic Vector class
/// vector for a vector field, e.g.: velocity,
/// here also happens the fortran wrapping
/// \ingroup Field
template<class S, class O>
class VectorField {

	template<class S1, class O1>
	friend class Grad;
	template<class S1,class O1>
	friend class Div;
	template<class S1,class O1>
	friend class Helmholtz;
	template<class S1,class O1>
	friend class Nonlinear;
	template<class S1,class O1>
	friend class NonlinearJacobian;

public:
	typedef S Scalar;
	typedef O Ordinal;
//	friend class Operator;
//	friend class Grad<Scalar,Ordinal>;

private:

//	using Teuchos::RCP;
	typedef Scalar* ScalarArray;
	typedef VectorField<Scalar,Ordinal> VF;

public:
	typedef Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > >  IndexSpaces;

	VectorField(): fieldS_(Teuchos::null),innerIS_(Teuchos::null),fullIS_(Teuchos::null),vec_(0) {};

	VectorField( const Teuchos::RCP<const FieldSpace<Ordinal> >& fieldS, IndexSpaces innerIS, IndexSpaces fullIS ):fieldS_(fieldS),innerIS_(innerIS),fullIS_(fullIS) {
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

	/// \brief copy constructor.
	///
	/// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
	/// \param sF
	/// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	VectorField(const VectorField& vF, ECopyType copyType=DeepCopy):
	    fieldS_(vF.fieldS_),innerIS_(vF.innerIS_),fullIS_(vF.fullIS_) {

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

	Teuchos::RCP<VF> clone( ECopyType ctype=DeepCopy ) const {
	  return( Teuchos::rcp( new VF( *this, ctype ) ) );
	}

  /// \name Attribute methods
  //@{

	/// \brief get \c FieldSpace.
	Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const { return( fieldS_ ); }


	/// \brief returns the length of Field.
	///
	/// the vector length is withregard to the inner points such that
	/// \f[ N_u = (N_x-1)(N_y-2)(N_z-?) \f]
	/// \f[ N_v = (N_x-2)(N_y-1)(N_z-?) \f]
	/// \f[ N_w = (N_x-?)(N_y-?)(N_z-?) \f]
	/// \return vect length \f[= N_u+N_v+N_w\f]
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
	/// \name Update methods
	//@{

	/// \brief Replace \c this with \f$\alpha A + \beta B\f$.
	///
	/// only inner points.
	void add( const Scalar& alpha, const VF& A, const Scalar& beta, const VF& B ) {
		// add test for consistent VectorSpaces in debug mode
		for( int i=0; i<dim(); ++i )
			SF_add(
						nLoc(0), nLoc(1), nLoc(2),
						sInd(0,i), sInd(1,i), sInd(2,i),
						eInd(0,i), eInd(1,i), eInd(2,i),
						bl(0),   bl(1),   bl(2),
						bu(0),   bu(1),   bu(2),
						vec_[i], A.vec_[i], B.vec_[i], alpha, beta);
	}


   /// \brief Put element-wise absolute values of source vector \c y into this
   /// vector.
   ///
   /// Here x represents this vector, and we update it as
   /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
   /// \return Reference to this object
  void abs(const VF& y) {
    for( int i=0; i<dim(); ++i )
      SF_abs(
          nLoc(0), nLoc(1), nLoc(2),
          sInd(0,i), sInd(1,i), sInd(2,i),
          eInd(0,i), eInd(1,i), eInd(2,i),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[i], y.vec_[i] );
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal(const VF& y){
    // add test for consistent VectorSpaces in debug mode
    for( int i=0; i<dim(); ++i)
      SF_reciprocal(
          nLoc(0), nLoc(1), nLoc(2),
          sInd(0,i), sInd(1,i), sInd(2,i),
          eInd(0,i), eInd(1,i), eInd(2,i),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[i], y.vec_[i] );
  }


	/// \brief Scale each element of the vectors in \c this with \c alpha.
	void scale( const Scalar& alpha ) {
		for(int i=0; i<dim(); ++i)
			SF_scale(
					nLoc(0), nLoc(1), nLoc(2),
					sInd(0,i), sInd(1,i), sInd(2,i),
					eInd(0,i), eInd(1,i), eInd(2,i),
					bl(0),   bl(1),   bl(2),
					bu(0),   bu(1),   bu(2),
					vec_[i], alpha);
	}


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void scale(const VF& a) {
    // add test for consistent VectorSpaces in debug mode
    for(int i=0; i<dim(); ++i)
      SF_scale2(
          nLoc(0), nLoc(1), nLoc(2),
					sInd(0,i), sInd(1,i), sInd(2,i),
					eInd(0,i), eInd(1,i), eInd(2,i),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[i], a.vec_[i] );
  }


	/// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
	/// \todo add test in debuging mode for testing equality of VectorSpaces
	Scalar dot ( const VF& a ) const {
		Scalar b;
		VF_dot(
				commf(), dim(),
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0,0), sInd(1,0), sInd(2,0),
				eInd(0,0), eInd(1,0), eInd(2,0),
				sInd(0,1), sInd(1,1), sInd(2,1),
				eInd(0,1), eInd(1,1), eInd(2,1),
				sInd(0,2), sInd(1,2), sInd(2,2),
				eInd(0,2), eInd(1,2), eInd(2,2),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				vec_[0],     vec_[1],   vec_[2],
				a.vec_[0], a.vec_[1], a.vec_[2],
				b);
		return( b );
	}


	//@}
  /// @name Norm method
  //@{

  /// \brief Compute the norm of each individual vector.
  ///
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2^2\f$, the \c i-th column of \c this.
  /// \attention the two norm is not the real two norm but its square
  /// \todo implement OneNorm
	Scalar norm(  Belos::NormType type = Belos::TwoNorm ) const {
		bool twoNorm_yes = false;
  	bool infNorm_yes = false;

  	switch(type) {
  	case Belos::TwoNorm: twoNorm_yes = true; break;
  	case Belos::InfNorm: infNorm_yes = true; break;
  	case Belos::OneNorm: std::cout << "norm: not implemented"; return( 0. );
  	default: std::cout << "unkown norm"; return( 0. );
  	}

  	Scalar normvec;
  	VF_compNorm(
  	    commf(), dim(),
  	    nLoc(0), nLoc(1), nLoc(2),
				sInd(0,0), sInd(1,0), sInd(2,0),
				eInd(0,0), eInd(1,0), eInd(2,0),
				sInd(0,1), sInd(1,1), sInd(2,1),
				eInd(0,1), eInd(1,1), eInd(2,1),
				sInd(0,2), sInd(1,2), sInd(2,2),
				eInd(0,2), eInd(1,2), eInd(2,2),
  	    bl(0),   bl(1),   bl(2),
  			bu(0),   bu(1),   bu(2),
        vec_[0], vec_[1], vec_[2],
  			infNorm_yes, twoNorm_yes,
  				normvec, normvec );
  	return( normvec );
  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  double norm(const VF& weights) const {
    Scalar normvec;
    VF_weightedNorm(
  	    commf(), dim(),
  	    nLoc(0), nLoc(1), nLoc(2),
        sInd(0,0), sInd(1,0), sInd(2,0),
        eInd(0,0), eInd(1,0), eInd(2,0),
        sInd(0,1), sInd(1,1), sInd(2,1),
        eInd(0,1), eInd(1,1), eInd(2,1),
        sInd(0,2), sInd(1,2), sInd(2,2),
        eInd(0,2), eInd(1,2), eInd(2,2),
  	    bl(0),   bl(1),   bl(2),
  			bu(0),   bu(1),   bu(2),
        vec_[0], vec_[1], vec_[2],
        weights.vec_[0], weights.vec_[1], weights.vec_[2],
        normvec);
    return( normvec );
  }


  //@}
  /// \name Initialization methods
  //@{


  /// \brief mv := A
  ///
  /// Assign (deep copy) A into mv.
  void assign( const VF& a ) {
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


  /// \brief Replace the vectors with a random vectors.
  ///
  /// depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is requiered
  void random(bool useSeed = false, int seed = 1) {
  	for( int i=0; i<dim(); ++i )
			SF_random(
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0,i), sInd(1,i), sInd(2,i),
				eInd(0,i), eInd(1,i), eInd(2,i),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				vec_[i] );
  }

  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
  	for( int i=0; i<dim(); ++i )
			SF_init(
				nLoc(0), nLoc(1), nLoc(2),
				sInd(0,i), sInd(1,i), sInd(2,i),
				eInd(0,i), eInd(1,i), eInd(2,i),
				bl(0),   bl(1),   bl(2),
				bu(0),   bu(1),   bu(2),
				vec_[i], alpha);
  }


  /// \brief Replace each element of the vector \c vec[i] with \c alpha[i].
  void init( const Teuchos::Tuple<Scalar,3>& alpha ) {
   	for( int i=0; i<dim(); ++i )
 			SF_init(
 				nLoc(0), nLoc(1), nLoc(2),
 				sInd(0,i), sInd(1,i), sInd(2,i),
 				eInd(0,i), eInd(1,i), eInd(2,i),
 				bl(0),   bl(1),   bl(2),
 				bu(0),   bu(1),   bu(2),
 				vec_[i], alpha[i]);
   }


  ///  \brief initializes VectorField with the initial field defined in Fortran
  void initField( EFlowProfile flowType = Poiseuille2D_inX, double re=1., double om=1., double px = 1. ) {
    switch(flowType) {
    case ZeroProf :
      VF_init_Zero(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
					vec_[0], vec_[1], vec_[2] );
      break;
    case Poiseuille2D_inX :
      VF_init_2DPoiseuilleX(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
					vec_[0], vec_[1], vec_[2] );
      break;
    case Poiseuille2D_inY :
      VF_init_2DPoiseuilleY(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2] );
      break;
    case Pulsatile2D_inXC :
      VF_init_2DPulsatileXC(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2], re, om, px);
      break;
    case Pulsatile2D_inYC :
      VF_init_2DPulsatileYC(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2], re, om, px);
      break;
    case Pulsatile2D_inXS :
      VF_init_2DPulsatileXS(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2], re, om, px);
      break;
    case Pulsatile2D_inYS:
      VF_init_2DPulsatileYS(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2], re, om, px);
      break;
    case Streaming2D:
      VF_init_Streaming(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2] );
      break;
    case Circle2D:
      VF_init_Circle(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2] );
      break;
    case RankineVortex2D:
      VF_init_RankineVortex(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2] );
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
			SF_print(
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

	/// \todo add good documetnation here
	/// @return
	const MPI_Fint& commf()                     const { return( fieldS_->commf_  ); }
	MPI_Comm        comm()                      const { return( fieldS_->comm_   ); }
	const int&      dim()                       const { return( fieldS_->dim_    ); }
	const Ordinal&  nGlo(int i)                 const { return( fieldS_->nGlo_[i] ); }
	const Ordinal&  nLoc(int i)                 const { return( fieldS_->nLoc_[i]) ; }
	const Ordinal&  sInd(int i, int fieldType)  const { return( innerIS_[fieldType]->sInd_[i] ); }
	const Ordinal&  eInd(int i, int fieldType)  const { return( innerIS_[fieldType]->eInd_[i] ); }
	const Ordinal&  sIndB(int i, int fieldType) const { return( fullIS_[fieldType]->sInd_[i] ); }
	const Ordinal&  eIndB(int i, int fieldType) const { return( fullIS_[fieldType]->eInd_[i] ); }
	const Ordinal&  bl(int i)                   const { return( fieldS_->bl_[i] ); }
	const Ordinal&  bu(int i)                   const { return( fieldS_->bu_[i] ); }

}; // end of class VectorField



/// \brief creates a vector field belonging to a \c FieldSpace and two \c IndexSpaces
template<class Scalar, class Ordinal>
Teuchos::RCP< VectorField<Scalar,Ordinal> > createVectorField(
		const Teuchos::RCP<const FieldSpace<Ordinal> >& fieldS,
		typename VectorField<Scalar,Ordinal>::IndexSpaces innerIS,
		const typename VectorField<Scalar,Ordinal>::IndexSpaces& fullIS ) {
	return( Teuchos::RCP<VectorField<Scalar,Ordinal> > (
				new VectorField<Scalar,Ordinal>( fieldS, innerIS, fullIS ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_VECTORFIELD_HPP
