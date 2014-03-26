#pragma once
#ifndef PIMPACT_MULTIHARMONICFIELD_HPP
#define PIMPACT_MULTIHARMONICFIELD_HPP

#include <vector>
#include <iostream>
//#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"

#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"



namespace Pimpact {


/// \brief important basic Vector class.
///
/// vector for wrapping many fields into one multiharmonic field
/// \ingroup Field
template<class Field>
class MultiHarmonicField {

public:

	typedef typename Field::Scalar Scalar;
	typedef typename Field::Ordinal Ordinal;

private:

	typedef MultiHarmonicField<Field> MV;

public:

	MultiHarmonicField(): field0_(Teuchos::null), fields_(Teuchos::null) {};

	MultiHarmonicField( const Teuchos::RCP<Field>& field0,
	    const Teuchos::RCP< MultiField< ModeField<Field> > >& fields):field0_(field0),fields_(fields) {};


	/// \brief copy constructor.
	///
	/// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
	/// @param vF
	/// @param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	MultiHarmonicField(const MultiHarmonicField& vF, ECopyType copyType=ShallowCopy):
		field0_( vF.field0_->clone(copyType) ),
		fields_( vF.fields_->clone(copyType) )
	{};


  Teuchos::RCP<MV> clone( ECopyType ctype=DeepCopy ) const {
    return( Teuchos::rcp( new MV( field0_->clone(ctype), fields_->clone(ctype) ) ) );
  }


  /// \name Attribute methods
  //@{

	/// \warning it is assumed that both fields have the same \c FieldSpace
	/// @return field space of \c field0_
	Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const {return( field0_->getFieldSpace() );}


	Field&                    get0Field()               { return( *field0_ ); }
	const Field&              getConst0Field()    const { return( *field0_ ); }
	Teuchos::RCP<Field>       get0FieldPtr()            { return( field0_ ); }
	Teuchos::RCP<const Field> getConst0FieldPtr() const { return( field0_ ); }


	ModeField<Field>&                     getField( int i )              { return( fields_->getField(i) ); }
	const ModeField<Field>&               getConstField( int i )   const { return( fields_->getConstField(i) ); }
	Teuchos::RCP<ModeField<Field> >       getFieldPtr( int i )           { return( fields_->getFieldPtr(i) ); }
	Teuchos::RCP<const ModeField<Field> > getConstFieldPtr( int i) const { return( fields_->getConstFieldPtr(i) ); }


	Field&              getCField( int i )            { return( fields_->getFieldPtr(i)->getCField() ); }
	const Field&        getConstCField( int i ) const { return( fields_->getConstFieldPtr(i)->getConstCField() ); }
	Teuchos::RCP<Field> getCFieldPtr( int i )         { return( fields_->getFieldPtr(i)->getCFieldPtr() ); }


	Field&              getSField( int i )            { return( fields_->getFieldPtr(i)->getSField() ); }
	const Field&        getConstSField( int i ) const { return( fields_->getConstFieldPtr(i)->getConstSField() ); }
	Teuchos::RCP<Field> getSFieldPtr( int i )         { return( fields_->getFieldPtr(i)->getSFieldPtr() ); }


	/// \brief returns the length of Field.
	///
	/// the vector length is with regard to the inner points
	Ordinal getLength( bool nox_vec=false ) const {
  	return( field0_->getLength(nox_vec) + fields_->getLength(true) );
  }


	/// \brief get number of stored Field's
	/// \todo what makes sense here?
	int getNumberVecs() const { return( 1 ); }

	/// \brief get number of mode Field's
	/// \todo what makes sense here?
	int getNumberModes() const { return( fields_->getNumberVecs() ); }


	//@}
	/// \name Update methods
	//@{

	/// \brief Replace \c this with \f$\alpha A + \beta B\f$.
	void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {
		// add test for consistent VectorSpaces in debug mode
		field0_->add(alpha, *A.field0_, beta, *B.field0_);
		fields_->add(alpha, *A.fields_, beta, *B.fields_);
	}


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void abs(const MV& y) {
      field0_->abs( *y.field0_ );
      fields_->abs( *y.fields_ );
  }


   /// \brief Put element-wise reciprocal of source vector \c y into this vector.
   ///
   /// Here x represents this vector, and we update it as
   /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
   /// \return Reference to this object
   void reciprocal(const MV& y){
      field0_->reciprocal( *y.field0_ );
      fields_->reciprocal( *y.fields_ );
   }


	/// \brief Scale each element of the vectors in \c this with \c alpha.
	void scale( const Scalar& alpha ) {
		field0_->scale(alpha);
		fields_->scale(alpha);
	}


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  /// \todo implement me
  void scale(const MV& a) {
    field0_->scale( *a.field0_ );
    fields_->scale( *a.fields_ );
  }


	/// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
	/// \todo has to be redone
	Scalar dot ( const MV& a ) const {
		return( field0_->dot( *a.field0_ ) + fields_->dot( *a.fields_ ) );
	}


	//@}
  /// \name Norm method
  //@{

  /// \brief Compute the norm of Field.
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  /// \todo implement OneNorm
  /// \todo implement in fortran
	Scalar norm(  Belos::NormType type = Belos::TwoNorm ) const {
		switch(type) {
		case Belos::TwoNorm:
		  return( std::sqrt( std::pow(field0_->norm(type),2) + std::pow(fields_->norm(type),2) ) );
		case Belos::InfNorm:
		  return( std::max(field0_->norm(type), fields_->norm(type) ) );
		case Belos::OneNorm:
		  std::cout << "!!! Warning Belos::OneNorm not implemented \n"; return(0.);
  	default: std::cout << "!!! Warning unknown Belos::NormType:\t" << type << "\n"; return(0.);
		}
	}


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  double norm(const MV& weights) const {
    return( field0_->norm( *weights.field0_)+fields_->norm( *weights.fields_ ) );
  }


  //@}
  /// \name Initialization methods
  //@{


  /// \brief mv := A
  /// Assign (deep copy) A into mv.
  void assign( const MV& a ) {
  	field0_->assign(*a.field0_);
  	fields_->assign(*a.fields_);
  }


  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed = false, int seed = 1) {
  	field0_->random();
  	fields_->random();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
  	field0_->init(alpha);
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

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os )  {
  	field0_->print( os );
  	fields_->print( os );
  }

  void write( int count=0 ) {
  	field0_->write(count);
  	for( int i=0; i<getNumberModes(); ++i ) {
  	  getCFieldPtr(i)->write( count+2*i+1 );
  	  getSFieldPtr(i)->write( count+2*i+2 );
  	}
  }

protected:
	Teuchos::RCP<Field> field0_;
	Teuchos::RCP< MultiField< ModeField<Field> > > fields_;

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

}; //class MultiHarmonicField


/// \brief creates a scalar/vector mode field(vector)
///
/// \param sVS scalar Vector Space to which returned vector belongs
/// \return field vector
/// \relates MultiHarmonicField
template<class Field>
Teuchos::RCP< MultiHarmonicField<Field> > createMultiHarmonicField(
    const Teuchos::RCP<Field>&  field0,
    const Teuchos::RCP<MultiField<ModeField<Field> > >& fields ) {
	return Teuchos::rcp(
				new MultiHarmonicField<Field>( field0, fields ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MULTIHARMONICFIELD_HPP
