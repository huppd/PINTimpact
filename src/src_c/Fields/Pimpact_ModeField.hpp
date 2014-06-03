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



namespace Pimpact {

/// \brief important basic Vector class
/// vector for wrapping 2 fields into one mode
/// \ingroup Field
template<class Field>
class ModeField {

public:
  typedef typename Field::Scalar Scalar;
  typedef typename Field::Ordinal Ordinal;

private:

  typedef ModeField<Field> MV;

public:


  ModeField(
      const Teuchos::RCP<Field>& fieldc=Teuchos::null,
      const Teuchos::RCP<Field>& fields=Teuchos::null):
        fieldc_(fieldc),fields_(fields) {};


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
  ModeField(const ModeField& vF, ECopyType copyType=DeepCopy):
    fieldc_( vF.fieldc_->clone(copyType) ),
    fields_( vF.fields_->clone(copyType) )
  {};


  Teuchos::RCP<MV> clone( ECopyType ctype=DeepCopy ) const {
    return( Teuchos::rcp( new MV( fieldc_->clone(ctype), fields_->clone(ctype) ) ) );
  }



  /// \name Attribute methods
  //@{

  /// \warning it is assumed that both fields have the same \c FieldSpace
  /// @return field space of \c fieldc_
  Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const {return( fieldc_->getFieldSpace() );}

  Teuchos::RCP<Field> getCFieldPtr() { return( fieldc_ ); }
  Teuchos::RCP<Field> getSFieldPtr() { return( fields_ ); }

  Teuchos::RCP<const Field> getConstCFieldPtr() const { return( fieldc_ ); }
  Teuchos::RCP<const Field> getConstSFieldPtr() const { return( fields_ ); }

  Field& getCField() { return( *fieldc_ ); }
  Field& getSField() { return( *fields_ ); }

  const Field& getConstCField() const { return( *fieldc_ ); }
  const Field& getConstSField() const { return( *fields_ ); }


  /// \brief returns the length of Field.
  ///
  /// should be the same as 2*fieldc_->getVecLength()
  /// the vector length is with regard to the inner points
  Ordinal getLength( bool nox_vec=false ) const {
    return( fieldc_->getLength(nox_vec) + fields_->getLength(nox_vec) );
  }


  /// \brief get number of stored Field's
  int getNumberVecs() const { return( 1 ); }


  //@}
  /// \name Update methods
  //@{

  /// \brief Replace \c this with \f$\alpha A + \beta B\f$.
  void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {
    // add test for consistent VectorSpaces in debug mode
    fieldc_->add(alpha, *A.fieldc_, beta, *B.fieldc_);
    fields_->add(alpha, *A.fields_, beta, *B.fields_);
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  /// \todo implement me
  void abs(const MV& y) {
    fieldc_->abs( *y.fieldc_ );
    fields_->abs( *y.fields_ );
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  /// \todo implement me
  void reciprocal(const MV& y){
    fieldc_->reciprocal( *y.fieldc_ );
    fields_->reciprocal( *y.fields_ );
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale( const Scalar& alpha ) {
    fieldc_->scale(alpha);
    fields_->scale(alpha);
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  /// \todo implement me
  void scale(const MV& a) {
    fieldc_->scale( *a.fieldc_ );
    fields_->scale( *a.fields_ );
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  Scalar dot ( const MV& a, bool global=true ) const {
    Scalar b=0.;
    b = fieldc_->dot( *a.fieldc_, false ) + fields_->dot( *a.fields_, false );
    if( global ) {
      Scalar b_global=0.;
      MPI_Allreduce( &b, &b_global, 1, MPI_REAL8, MPI_SUM, comm() );
      b = b_global;
    }
    return( b );
  }


  //@}
  /// \name Norm method
  //@{

  /// \brief Compute the norm of Field.
  ///
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  /// \todo implement OneNorm
  /// \todo implement in fortran
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {
    Scalar normvec=0;
    switch(type) {
    case Belos::TwoNorm:
      normvec = std::pow(fieldc_->norm(type,false),2) + std::pow(fields_->norm(type,false),2);
      if( global ) {
        Scalar normvec_global;
        MPI_Allreduce( &normvec, &normvec_global, 1, MPI_REAL8, MPI_SUM, comm() );
        normvec = normvec_global;
      }
      return( std::sqrt(normvec) );
    case Belos::InfNorm:
      normvec = std::max(fieldc_->norm(type,false), fields_->norm(type,false) ) ;
      if( global ) {
        Scalar normvec_global;
        MPI_Allreduce( &normvec, &normvec_global, 1, MPI_REAL8, MPI_MAX, comm() );
        normvec = normvec_global;
      }
      return( normvec );
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
    return( fieldc_->norm( *weights.fieldc_)+fields_->norm( *weights.fields_ ) );
  }


  //@}
  /// \name Initialization methods
  //@{

  /// \brief mv := A
  /// Assign (deep copy) A into mv.
  void assign( const MV& a ) {
    fieldc_->assign(*a.fieldc_);
    fields_->assign(*a.fields_);
  }

  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed = false, int seed = 1) {
    fieldc_->random();
    fields_->random();
  }

  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    fieldc_->init(alpha);
    fields_->init(alpha);
  }

  //@}
  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os )  {
    fieldc_->print( os );
    fields_->print( os );
  }

  void write( int count=0 ) {
    fieldc_->write(count);
    fields_->write(count+1);
  }

  MPI_Comm comm() const { return( fieldc_->comm() ); }

protected:
  Teuchos::RCP<Field> fieldc_;
  Teuchos::RCP<Field> fields_;


}; // end of class ModeField



/// \brief creates a scalar/vector mode field(vector)
///
/// \param sVS scalar Vector Space to which returned vector belongs
/// \return field vector
/// \relates ModeField
template<class Field>
Teuchos::RCP< ModeField<Field> > createModeField( const Teuchos::RCP<Field>&  fieldc, const Teuchos::RCP<Field>& fields ) {
  return( Teuchos::rcp(
      new ModeField<Field>( fieldc, fields ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MODEFIELD_HPP
