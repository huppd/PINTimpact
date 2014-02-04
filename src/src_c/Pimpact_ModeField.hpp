#pragma once
#ifndef PIMPACT_MODEFIELD_HPP
#define PIMPACT_MODEFIELD_HPP

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
template<class Field>
class ModeField {

//	template<class S1, class O1>
//	friend class Grad;
//	template<class S1,class O1>
//	friend class Div;
//	template<class S1,class O1>
//	friend class Helmholtz;

public:
	typedef typename Field::Scalar Scalar;
	typedef typename Field::Ordinal Ordinal;

private:

//	using Teuchos::RCP;
//	typedef Scalar* ScalarArray;
	typedef ModeField<Field> MV;

public:
//	typedef Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > >  IndexSpaces;

	ModeField(Teuchos::RCP<Field> fieldc, Teuchos::RCP<Field> fields):fieldc_(fieldc),fields_(fields) {};

	/** copy constructor
	 * shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
	 * @param sF
	 * @param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	 */
	ModeField(const ModeField& vF, ECopyType copyType=ShallowCopy):
		fieldc_( Teuchos::rcp( new Field(*vF.fieldc_,copyType) ) ),
		fields_( Teuchos::rcp( new Field(*vF.fields_,copyType) ) )
	{};

//	~ModeField() { for(int i=0; i<3; ++i) delete[] vec_[i];}

  //! \name Attribute methods
  //@{

	/**
	 * \warning it is assumed that both fields have the same \c FieldSpace
	 * @return field space of \c fieldc_
	 */
	Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const {return fieldc_->getFieldSpace();}

	Teuchos::RCP<Field> getFieldC() { return fieldc_; }
	Teuchos::RCP<Field> getFieldS() { return fields_; }

	Teuchos::RCP<const Field> getConstFieldC() const { return fieldc_; }
	Teuchos::RCP<const Field> getConstFieldS() const { return fields_; }

	/** \brief get Vect length
	 * shoud be the same as 2*fieldc_->getVecLength()
	 * the vector length is withregard to the inner points
	 * @return vector length
	 */
  Ordinal getVecLength() const {
  	return fieldc_->getVecLength() + fields_->getVecLength();
  }



	//@}

	//! @name Update methods
	//@{


	/*! \brief Replace \c this with \f$\alpha A + \beta B\f$.
	 */
	void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {
		// add test for consistent VectorSpaces in debug mode
		fieldc_->add(alpha, *A.fieldc_, beta, *B.fieldc_);
		fields_->add(alpha, *A.fields_, beta, *B.fields_);
	}

	/** \brief Scale each element of the vectors in \c this with \c alpha.
	 */
	void scale( const Scalar& alpha ) {
		fieldc_->scale(alpha);
		fields_->scale(alpha);
	}


	/** \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
	 * \todo has to be redone
	 */
	Scalar dot ( const MV& a ) const {
		return fieldc_->dot( *a.fieldc_ ) + fields_->dot( *a.fields_ );
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
  			case Belos::TwoNorm: return fieldc_->norm(type) + fields_->norm(type);
  			case Belos::InfNorm: return std::max(fieldc_->norm(type), fields_->norm(type) );
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
  	fieldc_->assign(*a.fieldc_);
  	fields_->assign(*a.fields_);
  }

  /** \brief Replace the vectors with a random vectors.
   */
  void random() {
  	fieldc_->random();
  	fields_->random();
  }

  /*! \brief Replace each element of the vector  with \c alpha.
   */
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
  	fieldc_->init(alpha);
  	fields_->init(alpha);
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
//   *  \brief initializes ModeField with the initial field defined in Fortran
//   */
//  void init_field() {
//  	for( int i=0; i<dim(); ++i )
//  	VF_init_field( vec_[0], vec_[1], vec_[2] );
//  }

  //@}

  void print()  {
  	fieldc_->print();
  	fields_->print();
  }

  void write( int count=0 ) {
  	fieldc_->write(count);
  	fields_->write(count+1);
  }

protected:
	Teuchos::RCP<Field> fieldc_;
	Teuchos::RCP<Field> fields_;

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

}; //class ModeField


/** @brief creates a scalar/vector mode field(vector)
 *
 * @param sVS scalar Vector Space to which returned vector belongs
 * @return field vector
 */
template<class Field>
Teuchos::RCP< ModeField<Field> > createModeField( Teuchos::RCP<Field>  fieldc, Teuchos::RCP<Field> fields ) {
	return Teuchos::RCP<ModeField<Field> > (
				new ModeField<Field>( fieldc, fields ) );
}


} // namespace Pimpact

#endif // PIMPACT_MODEFIELD_HPP
