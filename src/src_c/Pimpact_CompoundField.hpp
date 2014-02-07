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
	CompoundField(): vfield_(Teuchos::null),sfield_(Teuchos::null) {};

	CompoundField(Teuchos::RCP<VField> vfield, Teuchos::RCP<SField> sfield):vfield_(vfield),sfield_(sfield) {};

	/**
	 * \brief copy constructor.
	 *
	 * shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
	 * @param sF
	 * @param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	 */
	CompoundField(const CompoundField& vF, ECopyType copyType=ShallowCopy):
		vfield_( Teuchos::rcp( new VField(*vF.vfield_,copyType) ) ),
		sfield_( Teuchos::rcp( new SField(*vF.sfield_,copyType) ) )
	{};


	Teuchos::RCP<MV> clone( ECopyType ctype=DeepCopy ){
	  return( Teuchos::rcp( new MV( vfield_->clone(ctype), sfield_->clone(ctype) ) ) );
	}
//	~CompoundField() { for(int i=0; i<3; ++i) delete[] vec_[i];}

  //! \name Attribute methods
  //@{

	/**
	 * \warning it is assumed that both fields have the same \c FieldSpace
	 * @return field space of \c cfield_
	 */
	Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const { return( vfield_->getFieldSpace() );}

	Teuchos::RCP<VField> getVField() { return( vfield_ ); }
	Teuchos::RCP<SField> getSField() { return( sfield_ ); }

	Teuchos::RCP<const VField> getConstVField() const { return( vfield_ ); }
	Teuchos::RCP<const SField> getConstSField() const { return( sfield_ ); }

	/** \brief get Vect length
	 * shoud be the same as 2*vfield_->getVecLength()
	 * the vector length is withregard to the inner points
	 * @return vector length
	 */
	/// \brief returns the length of Field.
	Ordinal getLength( bool nox_vec=false ) const {
		return( vfield_->getLength(nox_vec) + sfield_->getLength(nox_vec) );
  }


	/// \brief get number of stored Field's
	int getNumberVecs() const { return( 1 ); }


	//@}
	/// \name Update methods
	//@{

	/**
	 * \brief Replace \c this with \f$\alpha A + \beta B\f$.
	 */
	void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {
		// add test for consistent VectorSpaces in debug mode
		vfield_->add(alpha, *A.vfield_, beta, *B.vfield_);
		sfield_->add(alpha, *A.sfield_, beta, *B.sfield_);
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
      vfield_->abs( *y.vfield_ );
      sfield_->abs( *y.sfield_ );
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
      vfield_->reciprocal( *y.vfield_ );
      sfield_->reciprocal( *y.sfield_ );
   }


	/**
	 * \brief Scale each element of the vectors in \c this with \c alpha.
	 */
	void scale( const Scalar& alpha ) {
		vfield_->scale(alpha);
		sfield_->scale(alpha);
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
		vfield_->scale( *a.vfield_ );
		sfield_->scale( *a.sfield_ );
  }


	/**
	 * \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
	 * \todo has to be redone
	 */
	Scalar dot ( const MV& a ) const {
		return( vfield_->dot( *a.vfield_ ) + sfield_->dot( *a.sfield_ ) );
	}


	//@}
  /// @name Norm method
  //@{

  /**
   * \brief Compute the norm of the field.
   * \todo implement OneNorm
   */
	Scalar norm(  Belos::NormType type = Belos::TwoNorm ) const {
		switch(type) {
		case Belos::TwoNorm: return( vfield_->norm(type) + sfield_->norm(type) );
		case Belos::InfNorm: return( std::max(vfield_->norm(type), sfield_->norm(type) ) );
		case Belos::OneNorm: std::cout << "!!! Warning Belos::OneNorm not implemented \n"; return(0.);
  	default: std::cout << "!!! Warning unknown Belos::NormType:\t" << type << "\n"; return(0.);
		}
	}


  /**
   * \brief Weighted 2-Norm.
   *
   * Here x represents this vector, and we compute its weighted norm as follows:
   * \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
   * \return \f$ \|x\|_w \f$
   */
  double norm(const MV& weights) const {
    return( vfield_->norm( *weights.vfield_)+sfield_->norm( *weights.sfield_ ) );
  }


  //@}
  /// \name Initialization methods
  //@{

  /**
   * \brief mv := A
   * Assign (deep copy) A into mv.
   */
  void assign( const MV& a ) {
  	vfield_->assign(*a.vfield_);
  	sfield_->assign(*a.sfield_);
  }


  /**
   * \brief Replace the vectors with a random vectors.
   */
  void random(bool useSeed = false, int seed = 1) {
  	vfield_->random();
  	sfield_->random();
  }


  /**
   * \brief Replace each element of the vector  with \c alpha.
   */
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
  	vfield_->init(alpha);
  	sfield_->init(alpha);
  }


  //@}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os )  {
  	vfield_->print( os );
  	sfield_->print( os );
  }


  void write( int count=0 ) {
  	vfield_->write(count);
  	sfield_->write(count+1);
  }

protected:
	Teuchos::RCP<VField> vfield_;
	Teuchos::RCP<SField> sfield_;

}; // end of class CompoundField


/**
 * \brief creates a compound vector+scalar field(vector).
 * \param sVS scalar Vector Space to which returned vector belongs
 * \return Field vector
 */
template<class VField, class SField>
Teuchos::RCP< CompoundField<VField,SField> > createCompoundField(
		Teuchos::RCP<VField>  vfield, Teuchos::RCP<SField> sfield ) {

	return( Teuchos::RCP<CompoundField<VField,SField> > (
				new CompoundField<VField,SField>( vfield, sfield ) ) );

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_COMPOUNDFIELD_HPP
