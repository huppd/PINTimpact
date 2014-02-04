#pragma once
#ifndef PIMPACT_COMPOUNDFIELD_HPP
#define PIMPACT_COMPOUNDFIELD_HPP

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


/** \brief important basic Vector class
 * vector for wrapping 2 fields into one mode
 */
template<class VField, class SField>
class CompoundField {

//	template<class S1, class O1>
//	friend class Grad;
//	template<class S1,class O1>
//	friend class Div;
//	template<class S1,class O1>
//	friend class Helmholtz;

public:
	typedef typename VField::Scalar Scalar;
	typedef typename VField::Ordinal Ordinal;

private:

//	using Teuchos::RCP;
//	typedef Scalar* ScalarArray;
	typedef CompoundField<VField,SField> MV;

public:
//	typedef Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > >  IndexSpaces;

	CompoundField(Teuchos::RCP<VField> vfield, Teuchos::RCP<SField> sfield):vfield_(vfield),sfield_(sfield) {};

	/** copy constructor
	 * shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
	 * @param sF
	 * @param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	 */
	CompoundField(const CompoundField& vF, ECopyType copyType=ShallowCopy):
		vfield_( Teuchos::rcp( new VField(*vF.vfield_,copyType) ) ),
		sfield_( Teuchos::rcp( new SField(*vF.sfield_,copyType) ) )
	{};

//	~CompoundField() { for(int i=0; i<3; ++i) delete[] vec_[i];}

  //! \name Attribute methods
  //@{

	/**
	 * \warning it is assumed that both fields have the same \c FieldSpace
	 * @return field space of \c cfield_
	 */
	Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const {return vfield_->getFieldSpace();}

	Teuchos::RCP<VField> getVField() { return vfield_; }
	Teuchos::RCP<SField> getSField() { return sfield_; }

	Teuchos::RCP<const VField> getConstVField() const { return vfield_; }
	Teuchos::RCP<const SField> getConstSField() const { return sfield_; }

	/** \brief get Vect length
	 * shoud be the same as 2*vfield_->getVecLength()
	 * the vector length is withregard to the inner points
	 * @return vector length
	 */
  Ordinal getVecLength() const {
  	return vfield_->getVecLength() + sfield_->getVecLength();
  }



	//@}

	//! @name Update methods
	//@{


	/*! \brief Replace \c this with \f$\alpha A + \beta B\f$.
	 */
	void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {
		// add test for consistent VectorSpaces in debug mode
		vfield_->add(alpha, *A.vfield_, beta, *B.vfield_);
		sfield_->add(alpha, *A.sfield_, beta, *B.sfield_);
	}

	/** \brief Scale each element of the vectors in \c this with \c alpha.
	 */
	void scale( const Scalar& alpha ) {
		vfield_->scale(alpha);
		sfield_->scale(alpha);
	}


	/** \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
	 * \todo has to be redone
	 */
	Scalar dot ( const MV& a ) const {
		return vfield_->dot( *a.vfield_ ) + sfield_->dot( *a.sfield_ );
	}

	//@}

  //! @name Norm method
  //@{

  /** \brief Compute the norm of each individual vector.
    Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
    \todo implement OneNorm
    \todo implement in fortran
  */
  Scalar norm(  Belos::NormType type = Belos::TwoNorm ) const {
  	switch(type) {
  			case Belos::TwoNorm: return vfield_->norm(type) + sfield_->norm(type);
  			case Belos::InfNorm: return std::max(vfield_->norm(type), sfield_->norm(type) );
  			case Belos::OneNorm: std::cout << "norm: not implemented"; return 0.;
  			default: std::cout << "unkown norm"; return 0.;
    	}
  }


  //@}
  //! @name Initialization methods
  //@{


  /** \brief mv := A
   * Assign (deep copy) A into mv.
   */
  void assign( const MV& a ) {
  	vfield_->assign(*a.vfield_);
  	sfield_->assign(*a.sfield_);
  }

  /** \brief Replace the vectors with a random vectors.
   */
  void random() {
  	vfield_->random();
  	sfield_->random();
  }

  /*! \brief Replace each element of the vector  with \c alpha.
   */
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
  	vfield_->init(alpha);
  	sfield_->init(alpha);
  }

//  /** \brief Replace each element of the vector \c vec[i] with \c alpha[i].
//    */
//   void init( const Teuchos::Tuple<Scalar,3>& alpha ) {
//   	for( int i=0; i<dim(); ++i )
// 			SV_init(
// 				nLoc(0), nLoc(1), nLoc(2),
// 				sInd(0,i), sInd(1,i), sInd(2,i),
// 				eInd(0,i), eInd(1,i), eInd(2,i),
// 				bl(0),   bl(1),   bl(2),
// 				bu(0),   bu(1),   bu(2),
// 				vec_[i], alpha[i]);
//   }


//  /**
//   *  \brief initializes CompoundField with the initial field defined in Fortran
//   */
//  void init_field() {
//  	for( int i=0; i<dim(); ++i )
//  	VF_init_field( vec_[0], vec_[1], vec_[2] );
//  }

  //@}

  void print()  {
  	vfield_->print();
  	sfield_->print();
  }

  void write( int count=0 ) {
  	vfield_->write(count);
  	sfield_->write(count+1);
  }

protected:
	Teuchos::RCP<VField> vfield_;
	Teuchos::RCP<SField> sfield_;

//	/**
//	 * \todo add good documetnation here
//	 * @return
//	 */
//	MPI_Fint& commf() const { return  fieldS_->commf_ ; }
//	MPI_Comm comm() const { return  fieldS_->comm_ ; }
//	const int& dim() const { return fieldS_->dim_; }
//	const Ordinal& nGlo(int i) const { return fieldS_->nGlo_[i]; }
//	const Ordinal& nLoc(int i) const { return  fieldS_->nLoc_[i]; }
//	const Ordinal& sInd(int i, int fieldType) const { return innerIS_[fieldType]->sInd_[i]; }
//	const Ordinal& eInd(int i, int fieldType) const { return innerIS_[fieldType]->eInd_[i]; }
//	const Ordinal& sIndB(int i, int fieldType) const { return fullIS_[fieldType]->sInd_[i]; }
//	const Ordinal& eIndB(int i, int fieldType) const { return fullIS_[fieldType]->eInd_[i]; }
//	const Ordinal& bl(int i) const { return fieldS_->bl_[i]; }
//	const Ordinal& bu(int i) const { return fieldS_->bu_[i]; }

}; //class CompoundField


/** @brief creates a compound vector+scalar field(vector)
 *
 * @param sVS scalar Vector Space to which returned vector belongs
 * @return field vector
 */
template<class VField, class SField>
Teuchos::RCP< CompoundField<VField,SField> > createCompoundField( Teuchos::RCP<VField>  vfield, Teuchos::RCP<SField> sfield ) {
	return Teuchos::RCP<CompoundField<VField,SField> > (
				new CompoundField<VField,SField>( vfield, sfield ) );
}


} // namespace Pimpact

#endif // PIMPACT_COMPOUNDFIELD_HPP
