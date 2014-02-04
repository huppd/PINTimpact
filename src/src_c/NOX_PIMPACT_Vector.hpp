#ifndef NOX_PIMPACT_VECTOR_H
#define NOX_PIMPACT_VECTOR_H

#include "NOX_Common.H"  // for NOX_Config.h
#include "Teuchos_RCP.hpp"

// Forward declaration of MultiVector class
namespace NOX {
  namespace PIMPACT {
    class MultiVector;
  }
}

//! Nonlinear solvers package namespace
namespace NOX {

  //! Specify whether to copy using deep copy or just copy by shape.
  enum CopyType {
    //! Copy object including all data
    DeepCopy,
    //! Copy the shape of the object only
    ShapeCopy
  };

/*! 
  \brief %NOX PIMPACT interface for vector and group

  The user should implement their own concrete implementation of this
  class or use one of the implementations provided by us, as defined
  in the NOX::Epetra and NOX::LAPACK namespaces.
*/
namespace PIMPACT {

/*! 
  \brief %NOX's pure PIMPACT vector interface for vectors that are
  used by the nonlinear solver. 

  This class is a member of the namespace NOX::PIMPACT.

  The user should implement their own concrete implementation of this
  class or use one of the implementations provided by us.

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

template<class ScalarType, int dim>
class Vector {

public:
  
  //! Norm types used in norm() calculations
  enum NormType {
    //! Use the 2-norm
    TwoNorm, 
    //! Use the 1-norm
    OneNorm, 
    //! Use the max-norm, a.k.a. the \f$\infty\f$-norm
    MaxNorm
  };

  //! %PIMPACT %Vector constructor (does nothing)
  Vector() {};

  //! %PIMPACT %Vector destructor (does nothing)
  virtual ~Vector() {};

  //@{ \name Initialization methods.

  //! Initialize every element of this vector with \c gamma.
  /*! 
    Here x represents this vector, and we update it as
    \f[ x_i = \gamma \quad \mbox{for } i=1,\dots,n \f] 
    \return Reference to this object
  */
  virtual NOX::PIMPACT::Vector& init(double gamma) = 0;

  //! Initialize each element of this vector with a random value
  /*!
    If \c useSeed is true, uses the value of \c seed to seed the
    random number generator before filling the entries of this
    vector. So, if two calls are made where \c useSeed is true and \c
    seed is the same, then the vectors returned should be the same.

    Default implementation throw an error. Only referenced by LOCA methods.

    \return Reference to this object
   */
  virtual NOX::PIMPACT::Vector& random(bool useSeed = false, int seed = 1);

  //! Put element-wise absolute values of source vector \c y into this vector.
  /*! 
    Here x represents this vector, and we update it as
    \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to this object
  */
  virtual NOX::PIMPACT::Vector& abs(const NOX::PIMPACT::Vector& y) = 0;

  //! Copy source vector \c y into this vector.
  /*! 
    Here x represents this vector, and we update it as
    \f[ x_i = y_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to this object
  */
  virtual NOX::PIMPACT::Vector& operator=(const NOX::PIMPACT::Vector& y) = 0;

  //! Put element-wise reciprocal of source vector \c y into this vector.
  /*! 
    Here x represents this vector, and we update it as
    \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f] 

    \return Reference to this object
  */
  virtual NOX::PIMPACT::Vector& reciprocal(const NOX::PIMPACT::Vector& y) = 0;

  //@}

  //@{ \name Update methods.

  //! Scale each element of this vector by \c gamma.
  /*! 
    Here x represents this vector, and we update it as
    \f[ x_i = \gamma x_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to this object
  */
  virtual NOX::PIMPACT::Vector& scale(double gamma) = 0;

  //! Scale this vector <em>element-by-element</em> by the vector a.
  /*! 
    Here x represents this vector, and we update it as
    \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to this object
  */
  virtual NOX::PIMPACT::Vector& scale(const NOX::PIMPACT::Vector& a) = 0;

  //! Compute x = (alpha * a) + (gamma * x) where x is this vector.
  /*! 
    Here x represents this vector, and we update it as
    \f[ x_i = \alpha \; a_i + \gamma \; x_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to this object
  */
  virtual NOX::PIMPACT::Vector& update(double alpha, const NOX::PIMPACT::Vector& a, double gamma = 0.0) = 0;

  //! Compute x = (alpha * a) + (beta * b) + (gamma * x) where x is this vector.
  /*! 
    Here x represents this vector, and we update it as
    \f[ x_i = \alpha \; a_i + \beta \; b_i + \gamma \; x_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to this object
  */
  virtual NOX::PIMPACT::Vector& update(double alpha, const NOX::PIMPACT::Vector& a,
			 double beta, const NOX::PIMPACT::Vector& b,
			 double gamma = 0.0) = 0;

  //@}

  //@{ \name Creating new Vectors.

  /*! 
    \brief Create a new %Vector of the same underlying type by
    cloning "this", and return a pointer to the new vector.  

    If type is NOX::DeepCopy, then we need to create an exact replica
    of "this". Otherwise, if type is NOX::ShapeCopy, we need only
    replicate the shape of "this" (the memory is allocated for the
    objects, but the current values are not copied into the
    vector). Note that there is <em>no assumption</em> that a vector
    created by ShapeCopy is initialized to zeros.

    \return Pointer to newly created vector or NULL if clone is not supported. 
  */
  virtual Teuchos::RCP<NOX::PIMPACT::Vector>
  clone(NOX::CopyType type = NOX::DeepCopy) const = 0;

  /*! 
   * \brief Create a MultiVector with \c numVecs+1 columns out of an array of 
   * Vectors.  The vector stored under \c this will be the first column with
   * the remaining \c numVecs columns given by \c vecs.
   *
   * The default implementation creates a generic NOX::MultiVector with
   * either Shape or Deep copies of the supplied vectors.
   */
  virtual Teuchos::RCP<NOX::PIMPACT::MultiVector>
  createMultiVector(const NOX::PIMPACT::Vector* const* vecs,
		    int numVecs, NOX::CopyType type = NOX::DeepCopy) const;

  /*! 
   * \brief Create a MultiVector with \c numVecs columns.  
   *
   * The default implementation creates a generic NOX::MultiVector with
   * either Shape or Deep copies of the supplied vector.
   */
  virtual Teuchos::RCP<NOX::PIMPACT::MultiVector>
  createMultiVector(int numVecs, NOX::CopyType type = NOX::DeepCopy) const;

  //@}

  //@{ \name Norms.

  //! Norm.
  /*! 
    Here x represents this vector, and we compute its norm as follows:
    for each NOX::PIMPACT::Vector::NormType:
    <ul>
    <li>  NOX::PIMPACT::Vector::TwoNorm  \f[ \|x\| = \sqrt{\sum_{i=1}^{n} x_i^2} \f]
    <li>  NOX::PIMPACT::Vector::OneNorm  \f[ \|x\| = \sum_{i=1}^{n} |x_i| \f]
    <li>  NOX::PIMPACT::Vector::MaxNorm  \f[ \|x\| = \max_{i} |x_i| \f]
    </uL>

    \return \f$\|x\|\f$
  */
  virtual double norm(NOX::PIMPACT::Vector::NormType type = NOX::PIMPACT::Vector::TwoNorm) const = 0;

  //! Weighted 2-Norm.
  /*! 
    Here x represents this vector, and we compute its weighted norm as follows:
    \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f] 
    \return \f$ \|x\|_w \f$
  */
  virtual double norm(const NOX::PIMPACT::Vector& weights) const = 0;

  //@}

  //@{ \name Inner product.

  //! Inner product with \c y.
  /*! 
    Here x represents this vector, and we compute its inner product with y as follows:
    \f[ \langle x,y \rangle = \sum_{i=1}^n x_i y_i \f] 
    \return \f$\langle x,y \rangle\f$
   */
  virtual double innerProduct(const NOX::PIMPACT::Vector& y) const = 0;
  
  //@}

  //! Return the length of vector.
  /*! 
    \return The length of this vector
    \note Even if the vector is distributed across processors, this
    should return the <em> global length </em> of the vector.
   */
  virtual int length() const = 0;

  //! Print the vector.  To be used for debugging only.  
  virtual void print(std::ostream& stream) const;

protected:
  RCP<FieldSpace> sVS;

  ScalarType* field;

}; // class Vector
} // namespace PIMPACT
} // namespace NOX

#endif
