#pragma once
#ifndef PIMPACT_MULTIHARMONICFIELD_HPP
#define PIMPACT_MULTIHARMONICFIELD_HPP

#include <vector>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_AbstractField.hpp"



namespace Pimpact {


/// \brief important basic Vector class.
///
/// vector for wrapping many fields into one multiharmonic field
/// \todo SpaceT constructor
/// \todo continous memory
/// \ingroup Field
template<class Field>
class MultiHarmonicField : private AbstractField<typename Field::SpaceT> {

public:

  typedef typename Field::SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

protected:

  typedef MultiHarmonicField<Field> MV;

  typedef AbstractField<SpaceT> AF;

  Teuchos::RCP<Field> field0_;
  Teuchos::RCP< MultiField< ModeField<Field> > > fields_;

public:

  MultiHarmonicField(
      const Teuchos::RCP<const SpaceT>& space ):
        AF( space ),
        field0_( create<Field>(space) ),
        fields_( create<MultiField< ModeField<Field> > >(space) ) {};

  /// \deprecated
  MultiHarmonicField(
      const Teuchos::RCP<Field>& field0,
      const Teuchos::RCP< MultiField< ModeField<Field> > >& fields):
        AF( field0->space() ),
        field0_(field0),
        fields_(fields) {};

protected:

  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
  MultiHarmonicField( const MultiHarmonicField& vF, ECopyType copyType=DeepCopy ):
    AF( vF.space() ),
    field0_( vF.field0_->clone(copyType) ),
    fields_( vF.fields_->clone(copyType) )
{};

public:

  Teuchos::RCP<MV> clone( ECopyType ctype=DeepCopy ) const {
    return( Teuchos::rcp( new MV( field0_->clone(ctype), fields_->clone(ctype) ) ) );
  }

  void push_back( const Teuchos::RCP< ModeField<Field> >& modeField=Teuchos::null ) {
    fields_->push_back( modeField );
  }

  /// \name Attribute methods
  /// \{

  /// \warning it is assumed that both fields have the same \c StencilWidths


  Field&                    get0Field()               { return( *field0_ ); }
  const Field&              getConst0Field()    const { return( *field0_ ); }
  Teuchos::RCP<      Field> get0FieldPtr()            { return(  field0_ ); }
  Teuchos::RCP<const Field> getConst0FieldPtr() const { return(  field0_ ); }


  ModeField<Field>&                     getField        ( int i )       { return( fields_->getField        (i) ); }
  const ModeField<Field>&               getConstField   ( int i ) const { return( fields_->getConstField   (i) ); }
  Teuchos::RCP<      ModeField<Field> > getFieldPtr     ( int i )       { return( fields_->getFieldPtr     (i) ); }
  Teuchos::RCP<const ModeField<Field> > getConstFieldPtr( int i ) const { return( fields_->getConstFieldPtr(i) ); }


  Field&                    getCField        ( int i )       { return( fields_->getFieldPtr     (i)->getCField()      ); }
  const Field&              getConstCField   ( int i ) const { return( fields_->getConstFieldPtr(i)->getConstCField() ); }
  Teuchos::RCP<      Field> getCFieldPtr     ( int i )       { return( fields_->getFieldPtr     (i)->getCFieldPtr()   ); }
  Teuchos::RCP<const Field> getConstCFieldPtr( int i ) const { return( fields_->getConstFieldPtr(i)->getConstCFieldPtr() ); }


  Field&                    getSField        ( int i )       { return( fields_->getFieldPtr     (i)->getSField()         ); }
  const Field&              getConstSField   ( int i ) const { return( fields_->getConstFieldPtr(i)->getConstSField()    ); }
  Teuchos::RCP<      Field> getSFieldPtr     ( int i )       { return( fields_->getFieldPtr     (i)->getSFieldPtr()      ); }
  Teuchos::RCP<const Field> getConstSFieldPtr( int i ) const { return( fields_->getConstFieldPtr(i)->getConstSFieldPtr() ); }


  Teuchos::RCP<const SpaceT> space() const { return( AF::space_ ); }

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


  /// \}
  /// \name Update methods
  /// \{

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
  void scale(const MV& a) {
    field0_->scale( *a.field0_ );
    fields_->scale( *a.fields_ );
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  Scalar dot ( const MV& a, bool global=true ) const {

    Scalar b = 0.;

    b = field0_->dot( *a.field0_, false ) + fields_->dot( *a.fields_, false );

    if( global ) this->reduceNorm( comm(), b );

    return( b );

  }


  /// \}
  /// \name Norm method
  /// \{

  /// \brief Compute the norm of Field.
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = 0.;

    switch(type) {
    case Belos::OneNorm:
      normvec = field0_->norm(type,false) + fields_->norm(type,false);
      break;
    case Belos::TwoNorm:
      normvec = field0_->norm(type,false) + fields_->norm(type,false);
      break;
    case Belos::InfNorm:
      normvec = std::max( field0_->norm(type,false), fields_->norm(type,false) );
      break;
    }

    if( global ) this->reduceNorm( comm(), normvec, type );

    return( normvec );

  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  double norm( const MV& weights, bool global=true ) const {

    double normvec=field0_->norm(*weights.field0_,false)+fields_->norm(*weights.fields_,false);

    if( global ) this->reduceNorm( comm(), normvec, Belos::TwoNorm );

    return( normvec );

  }


  /// \}
  /// \name Initialization methods
  /// \{


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

  /// \}

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

  MPI_Comm comm() const { return( field0_->space()->comm() ); }




}; // end of class MultiHarmonicField



/// \brief creates a scalar/vector mode field(vector)
///
/// \param field0 scalar Vector Space to which returned vector belongs
/// \param fields vector
/// \relates MultiHarmonicField
template<class Field>
Teuchos::RCP< MultiHarmonicField<Field> > createMultiHarmonicField(
    const Teuchos::RCP<Field>&  field0,
    const Teuchos::RCP<MultiField<ModeField<Field> > >& fields ) {
  return(
      Teuchos::rcp(
          new MultiHarmonicField<Field>( field0, fields ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIHARMONICFIELD_HPP
