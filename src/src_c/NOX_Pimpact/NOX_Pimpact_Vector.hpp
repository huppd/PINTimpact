#pragma once
#ifndef NOX_PIMPACT_VECTOR_HPP
#define NOX_PIMPACT_VECTOR_HPP

#include "NOX_Common.H"  // for NOX_Config.h
#include "NOX_Abstract_Vector.H"

#include "Teuchos_RCP.hpp"

//#include "BelosTypes.hpp"

#include "Pimpact_Types.hpp"
#include "Pimpact_MultiField.hpp"



// Forward declaration of MultiVector class
namespace NOX{
namespace Pimpact{
class MultiVector;
}
}



//! Nonlinear solvers package namespace
namespace NOX{

/// \brief %NOX PIMPACT interface for Vector
namespace Pimpact {

/// \brief %NOX's pure PIMPACT vector interface for vectors that are used by the
/// nonlinear solver.
///
/// This class is a member of the namespace NOX::PIMPACT.
/// \note actually this class wrapps \c MultiField, where \c MultiField is
/// considered as one Vector, this is relevent for norms, dots and length
template<class Field>
class Vector : public virtual NOX::Abstract::Vector {

public:

	/// %PIMPACT %Vector constructor (does nothing)
	Vector() { field_ = Teuchos::null; };

	/// constructor from \c Field
	Vector( const Teuchos::RCP<Field>& field):field_(field) {};

	/// %PIMPACT %Vector destructor
	virtual ~Vector() { field_ = Teuchos::null; };


	//@{ \name Initialization methods.

  /// \brief Initialize every element of this vector with \c gamma.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = \gamma \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  virtual NOX::Abstract::Vector& init( double gamma=0 ) {
  	field_->init( static_cast<typename Field::Scalar>(gamma) );
  	return( *this );
  }


  /// \brief Initialize each element of this vector with a random value.
  ///
  /// If \c useSeed is true, uses the value of \c seed to seed the
  /// random number generator before filling the entries of this
  /// vector. So, if two calls are made where \c useSeed is true and \c
  /// seed is the same, then the vectors returned should be the same.
  /// Default implementation throw an error. Only referenced by LOCA methods.
  /// \return Reference to this object
  virtual NOX::Abstract::Vector& random(bool useSeed = false, int seed = 1) {
  	field_->random( useSeed, seed );
  	return( *this );
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  virtual NOX::Abstract::Vector& abs(const Vector<Field>& y) {
    field_->abs( *y.field_ );
  	return( *this );
  }
  virtual NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y) {
  	return( abs( dynamic_cast<const Vector<Field>& >(y) ) );
  }


  /// \brief Copy source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = y_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  virtual NOX::Abstract::Vector& operator=(const Vector<Field>& y) {
  	field_->assign( *y.field_ );
  	return( *this );
  }
  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y) {
  	return( operator=( dynamic_cast<const Vector<Field>& >(y) ) );
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  virtual NOX::Abstract::Vector& reciprocal(const Vector<Field>& y){
    field_->reciprocal( *y.field_ );
  	return( *this );
  }
  virtual NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y){
  	return( reciprocal( dynamic_cast<const Vector<Field>& >(y) ) );
  }


  //@}
  //@{ \name Update methods.

  /// \brief Scale each element of this vector by \c gamma.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = \gamma x_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  virtual NOX::Abstract::Vector& scale(double gamma) {
  	field_->scale( static_cast<typename Field::Scalar>(gamma) );
  	return( *this );
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  virtual NOX::Abstract::Vector& scale(const Vector<Field>& a) {
    field_->scale( *a.field_ );
  	return( *this );
  }
  virtual NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a) {
  	return( scale( dynamic_cast<const Vector<Field>& >(a) ) );
  }


  /// \brief Compute x = (alpha * a) + (gamma * x) where x is this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = \alpha \; a_i + \gamma \; x_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  /// \todo test me good
  virtual NOX::Abstract::Vector& update(double alpha, const Vector<Field>& a, double gamma = 0.0) {
  	field_->add( alpha, *a.field_, gamma, *field_);
  	return( *this );
  }
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, double gamma = 0.0) {
  	return( update( alpha, dynamic_cast<const Vector<Field>& >(a), gamma ) );
  }


  /// \brief Compute x = (alpha * a) + (beta * b) + (gamma * x) where x is this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = \alpha \; a_i + \beta \; b_i + \gamma \; x_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  /// \todo test me
  virtual NOX::Abstract::Vector& update(
  		double alpha, const Vector<Field>& a,
			double beta, const Vector<Field>& b,
			double gamma = 0.0) {
  	field_->add( alpha, *a.field_, gamma,* field_);
  	field_->add( beta,  *b.field_, 1.,* field_);
  	return( *this);
  }
  virtual NOX::Abstract::Vector& update(
  		double alpha, const NOX::Abstract::Vector& a,
			double beta, const NOX::Abstract::Vector& b,
			double gamma = 0.0) {
  	return(
  			update(
  					alpha, dynamic_cast<const Vector<Field>& >(a),
  					beta,  dynamic_cast<const Vector<Field>& >(b),
  					gamma ) );
  }


  //@}
  //@{ \name Creating new Vectors.

  /// \brief Create a new %Vector of the same underlying type by
  /// cloning "this", and return a pointer to the new vector.
  ///
  /// If type is NOX::DeepCopy, then we need to create an exact replica of
  /// "this". Otherwise, if type is NOX::ShapeCopy, we need only replicate the
  /// shape of "this" (the memory is allocated for the objects, but the current
  /// values are not copied into the vector). Note that there is <em>no
  /// assumption</em> that a vector created by ShapeCopy is initialized to zeros.
  /// \return Pointer to newly created vector or NULL if clone is not supported.
  /// \todo implement
  /// \warning only Deep copy is working
  /// \todo fix schallow copy
  virtual Teuchos::RCP<NOX::Abstract::Vector>
  clone(NOX::CopyType type = NOX::DeepCopy) const {
  	switch(type) {
  	case NOX::DeepCopy:
  		return( Teuchos::rcp(new Vector<Field>( field_->clone( ::Pimpact::DeepCopy) ) ) );
  	case NOX::ShapeCopy:
  		return( Teuchos::rcp(new Vector<Field>( field_->clone( ::Pimpact::ShallowCopy) ) ) );
  	default: // just to make the compliler happy
  	  return( Teuchos::null );
  	}

  }


  /// \brief Create a MultiVector with \c numVecs+1 columns out of an array of
  /// Vectors.  The vector stored under \c this will be the first column with
  /// the remaining \c numVecs columns given by \c vecs.
  ///
  /// The default implementation creates a generic NOX::MultiVector with
  /// either Shape or Deep copies of the supplied vectors.
//  virtual Teuchos::RCP<NOX::Pimpact::MultiVector>
//  createMultiVector(const Vector<Field>* const* vecs,
//		    int numVecs, NOX::CopyType type = NOX::DeepCopy) const;


  /// \brief Create a MultiVector with \c numVecs columns.
  ///
  /// The default implementation creates a generic NOX::MultiVector with
  /// either Shape or Deep copies of the supplied vector.
//  virtual Teuchos::RCP<NOX::Pimpact::MultiVector>
//  createMultiVector(int numVecs, NOX::CopyType type = NOX::DeepCopy) const;


  //@}
  //@{ \name Norms.

  /// \brief Norm.
  ///
  /// Here x represents this vector, and we compute its norm as follows:
  /// for each NOX::PIMPACT::Vector::NormType:
  /// <ul>
  /// <li>  NOX::PIMPACT::Vector::TwoNorm  \f[ \|x\| = \sqrt{\sum_{i=1}^{n} x_i^2} \f]
  /// <li>  NOX::PIMPACT::Vector::OneNorm  \f[ \|x\| = \sum_{i=1}^{n} |x_i| \f]
  /// <li>  NOX::PIMPACT::Vector::MaxNorm  \f[ \|x\| = \max_{i} |x_i| \f]
  /// </uL>
  /// \return \f$\|x\|\f$
  virtual double norm( NOX::Abstract::Vector::NormType type=NOX::Abstract::Vector::TwoNorm) const {
  	switch( type ) {
  	case TwoNorm: return( field_->norm( Belos::TwoNorm ) );
  	case OneNorm: return( field_->norm( Belos::OneNorm ) );
  	case MaxNorm: return( field_->norm( Belos::InfNorm ) );
  	default: std::cout << "!!! Warning unknown NOX::Pimpact::Vector::NormType:\t" << type << "\n"; return(0.);
  	}
  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  virtual double norm(const Vector<Field>& weights) const {
  	return( field_->norm( *weights.field_) );
  }
  virtual double norm(const NOX::Abstract::Vector& weights) const {
  	return( norm( dynamic_cast<const Vector<Field>& >(weights) ) );
  }


  //@}
  //@{ \name Inner product.

  /// \brief Inner product with \c y.
  ///
  /// Here x represents this vector, and we compute its inner product with y as
  /// follows:
  /// \f[ \langle x,y \rangle = \sum_{i=1}^n x_i y_i \f]
  /// \return \f$\langle x,y \rangle\f$
  virtual double innerProduct( const Vector<Field>& y ) const {
  	return( field_->dot( *y.field_ ) );
  }
  virtual double innerProduct( const NOX::Abstract::Vector& y ) const {
  	return( innerProduct( dynamic_cast<const Vector<Field>& >(y) ) );
  }
  

  //@}

  /// \brief Return the length of vector.
  ///
  /// \return The length of this vector
  /// \note Even if the vector is distributed across processors, this
  /// should return the <em> global length </em> of the vector.
  virtual int length() const {
  	return( field_->getLength(true) );
  }


  /// Print the vector.  To be used for debugging only.
  virtual void print( std::ostream& os ) const {
  	field_->print( os );
  }


  Teuchos::RCP<Field> getFieldPtr() { return( field_ ); };
  Teuchos::RCP<const Field> getConstFieldPtr() const { return( field_ ); };
  Field& getField() { return( *field_ ); };
  const Field& getConstField() const { return( *field_ ); };

protected:
  Teuchos::RCP<Field> field_;

}; // end of class Vector



/// \relates Vector
template<class Field>
Teuchos::RCP< Vector<Field> > createVector( const Teuchos::RCP<Field>& field ){
  return( Teuchos::rcp( new Vector<Field>(field) ) );
}


} // end of namespace Pimpact
} // end of namespace NOX

#endif //end of #ifndef NOX_PIMPACT_VECTOR_HPP

