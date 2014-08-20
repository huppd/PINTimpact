#pragma once
#ifndef PIMPACT_MULTIFIELD_HPP
#define PIMPACT_MULTIFIELD_HPP


/// \file Pimpact_MultiField.hpp

#include <vector>
#include <Teuchos_Array.hpp>
#include "Teuchos_RCP.hpp"
#include <Teuchos_Range1D.hpp>

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_AbstractField.hpp"



namespace Pimpact {


/// \brief templated class which is the interface to \c Belos and \c NOX
///
/// has multiple \c Field's, where \c Field can be a \c Pimpact:ScalarField, \c
/// Pimpact:VectorField or some combination with \c Pimpact::ModeField or \c
/// Pimpact::CompoundField \note if this is heavily used for many Field's, then
/// the implementation should be improved such that communication is done such
/// that only done once per MV not per Field
///
/// \note for better documentation, look at the equivalent documentation in the \c Belos::...
/// \ingroup Field
template<class Field>
class MultiField : private AbstractField<typename Field::Scalar, typename Field::Ordinal> {

public:

  typedef typename Field::Scalar Scalar;
  typedef typename Field::Ordinal Ordinal;

private:

  typedef Pimpact::MultiField<Field> MV;
  Teuchos::Array<Teuchos::RCP<Field> > mfs_;

public:
  MultiField():mfs_(Teuchos::null) {};

  /// \brief constructor taking a \c Field constructing multiple shallow copys.
  /// \note maybe hide and make it private
  MultiField( const Field& field, const int numvecs, ECopyType ctyp=ShallowCopy ):mfs_(numvecs) {
    for( int i=0; i<numvecs; ++i )
      mfs_[i] = field.clone(ctyp);
  }


  /// \brief cheap constructor from one Field.
  ///
  /// creates simple wrapper from Field(no coppying).
  MultiField( const Teuchos::RCP<Field>& field ):mfs_(1) {
    mfs_[0] = field;
  }


  /// \brief copy constructor creating a view
  /// note clear if here only referencing or copy is happening
  MultiField( const MV& mv ):mfs_(mv.mfs_) {}


  /// \brief  constructor, creates \c numvecs  empty Fields
  /// \param numvecs
  /// @return
  MultiField( int numvecs ):mfs_(numvecs) {}


  /// \brief Create a new \c MultiField with \c numvecs columns.
  ///
  /// The returned Pimpact::MultiField has the same (distribution over one or
  /// more parallel processes) Its entries are not initialized and have undefined
  /// values.
  Teuchos::RCP< MV > clone( const int numvecs ) const {
    return( Teuchos::rcp( new MV( *mfs_[0], numvecs ) ) );
  }


  /// \brief Create a new \c MultiField with \c numvecs columns.
  ///
  /// The returned Pimpact::MultiField has the same (distribution over one or
  /// more parallel processes) Its entries are not initialized and have undefined
  /// values.
  Teuchos::RCP< MV > clone( ECopyType ctype = DeepCopy ) const {
    auto mv_ = Teuchos::rcp( new MV(getNumberVecs()) );
    for( int i=0; i<getNumberVecs(); ++i ) {
      mv_->mfs_[i] = mfs_[i]->clone(ctype);
    }
    return( mv_ );
  }


  /// \return deep copy of \c MultiField
  Teuchos::RCP<MV> CloneCopy() const {
    return( clone() );
  }


  /// \brief deep copy of \c index fields, the new \c MultiFields stores \c index.size() vector.
  ///
  /// \param index
  /// \return
  Teuchos::RCP<MV>  CloneCopy( const std::vector<int>& index ) const {
    auto mv_ = Teuchos::rcp( new MV(index.size()) );
    for( unsigned int i=0; i<index.size(); ++i ) {
      mv_->mfs_[i] = Teuchos::rcp( new Field( *mfs_[ index[i] ], DeepCopy ) );
    }
    return( mv_ );
  }


  /// deep copy of \c index fields, the new \c MultiFields stores \c index.size()
  /// \param index here index means an interval
  /// @return
  Teuchos::RCP<MV> CloneCopy( const Teuchos::Range1D& index) const {
    auto mv_ = Teuchos::rcp( new MV(index.size()) );
    int j = 0;
    for( int i=index.lbound(); i<=index.ubound(); ++i ) {
      mv_->mfs_[j] = Teuchos::rcp( new Field( *mfs_[ i ], DeepCopy ) );
      ++j;
    }
    return( mv_ );
  }


  /// \param index
  /// @return nonConst View
  Teuchos::RCP<MV> CloneViewNonConst( const std::vector<int>& index) {
    auto mv_ = Teuchos::rcp( new MV( index.size() ) );
    for( unsigned int i=0; i<index.size(); ++i ) {
      mv_->mfs_[i] =  mfs_[ index[i] ];
    }
    return( mv_ );
  }


  Teuchos::RCP<MV> CloneViewNonConst( const Teuchos::Range1D& index ) {
    auto mv_ = Teuchos::rcp( new MV(index.size()) );
    int j=0;
    for( int i=index.lbound(); i<=index.ubound(); ++i ) {
      mv_->mfs_[j++] =  mfs_[i];
    }
    return( mv_ );
  }


  Teuchos::RCP<const MV > CloneView( const std::vector<int>& index ) const {
    auto mv_ = Teuchos::rcp( new MV(index.size()) );
    for( unsigned int i=0; i<index.size(); ++i ) {
      mv_->mfs_[i] =  mfs_[ index[i] ];
    }
    return( mv_ );
  }


  Teuchos::RCP<const MV > CloneView( const Teuchos::Range1D& index ) const {
    auto mv_ = Teuchos::rcp( new MV(index.size()) );
    int j=0;
    for( int i=index.lbound(); i<=index.ubound(); ++i ) {
      mv_->mfs_[j++] =  mfs_[i];
    }
    return( mv_ );
  }


  /// \brief returns the length of Field.
  ///
  /// \param nox_vec if \c MultiField is used for NOX the Vector length is
  /// considered for all Fields
  Ordinal getLength( bool nox_vec=false ) const {
    if( nox_vec ){
      Ordinal n = 0;
      for( int i=0; i<getNumberVecs(); ++i )
        n += mfs_[i]->getLength(nox_vec);
      return( n );
    }
    else
      return( mfs_[0]->getLength(nox_vec) );
  }


  /// \brief get number of stored Field's
  int getNumberVecs() const { return( mfs_.size() ); }


  /// \brief is true
  bool HasConstantStride() const { return( true ); }


  //@}
  /// \name Update methods
  //@{

  /// \brief addes new field at end
  void push_back( const Teuchos::RCP<Field>& field=Teuchos::null ) {
    if( Teuchos::is_null(field) )
      mfs_.push_back( mfs_.back()->clone(ShallowCopy) );
    else
      mfs_.push_back( field );
  }



  /// \brief \f[ *this = \alpha A B + \beta *this \f]
  ///
  /// \param alpha scalar
  /// \param A Vector
  /// \param B Matrix
  /// \param beta
  void TimesMatAdd( const Scalar& alpha, const MV& A,
      const Teuchos::SerialDenseMatrix<int,Scalar>& B,
      const Scalar& beta ) {

    int m1 = A.getNumberVecs(); ///< is assumed to be equal to number vecs of this and ncolumns and nrows of b
    int m2 = getNumberVecs();   ///< is assumed to be equal to number vecs of this and ncolumns and nrows of b

    scale( beta );

    for( int i=0; i<m2; ++i ) {
      for( int j=0; j<m1; ++j ) {
        mfs_[i]->add( 1., *mfs_[i], alpha*B(j,i), *A.mfs_[j] );
      }
    }
  }



  /// \brief <tt>mv := alpha*A + beta*B</tt>
  ///
  ///	The Tpetra specialization of this method ignores and completely
  ///	overwrites any NaN or Inf entries in A.  Thus, it does <i>not</i> mean
  ///	the same thing as <tt>mv := 0*mv + alpha*A + beta*B</tt> in IEEE 754
  ///	floating-point arithmetic. (Remember that NaN*0 = NaN.)
  void add( Scalar alpha, const MV& A, Scalar beta, const MV& B ) {
    for( int i=0; i<getNumberVecs(); ++i )
      mfs_[i]->add( alpha, *A.mfs_[i], beta, *B.mfs_[i] );
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void abs(const MV& y) {
    for( int i=0; i<getNumberVecs(); ++i )
      mfs_[i]->abs( *y.mfs_[i] );
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal(const MV& y){
    for( int i=0; i<getNumberVecs(); ++i )
      mfs_[i]->reciprocal( *y.mfs_[i] );
  }


  /// \brief Scale each element of every \c Field by \c gamma.
  ///
  /// Here x represents on \c Field, and we update it as
  /// \f[ x_i = \alpha x_i \quad \mbox{for } i=1,\dots,n \f]
  void scale( const Scalar& alpha ) {
    for( int i=0; i<getNumberVecs(); ++i )
      mfs_[i]->scale(alpha);
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void scale(const MV& a) {
    for( int i=0; i<getNumberVecs(); ++i )
      mfs_[i]->scale( *a.mfs_[i] );
  }


  /// \brief Scale each element of a \c Field by a \c gamma.
  ///
  /// Here x_j represents the j'th field, and we update it as
  /// \f[ x_j[i] = \alpha_j x_j[i] \quad \mbox{for } i=1,\dots,n \f]
  void scale( const std::vector<Scalar>& alphas ) {
#ifdef DEBUG
    TEST_EQUALITY( alphas.size(), getNumberVecs() );
#endif
    for( unsigned int i=0; i<alphas.size(); ++i )
      mfs_[i]->scale( alphas[i] );
  }


  /// \f[ C_{ij} = \alpha A_i^T \c this_j \f]
  /// \param alpha
  /// \param A
  /// \param C
  void Trans( Scalar alpha, const MV& A,
      Teuchos::SerialDenseMatrix<int,Scalar>& C) const {

    const int n1 = getNumberVecs();
    const int n2 = A.getNumberVecs();
    //#ifdef DEBUG
    //    	if( n2!=C.numCols() )
    //    		std::cout << "n1!=C.numCols() n1= " << n1 << ", C.numcols() = " << C.numRows() << "\n";
    //    	if( n1!=C.numRows() )
    //    		std::cout << "n1!=C.numCols() n1= " << n2 << ", C.numcols() = " << C.numCols() << "\n";
    //    	TEST_EQUALITY( n2, C.numRows() );
    //    	TEST_EQUALITY( n1, C.numCols() );
    //#endif
    for( int i=0; i<n1; ++i)
      for( int j=0; j<n2; ++j )
        C(i,j) = alpha*mfs_[i]->dot( *A.mfs_[j] );
  }


  //@}

  /// For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
  /// \todo make reduction over vector
  void dot( const MV& A, std::vector<Scalar>& dots) const {

    const int n = getNumberVecs();
    Scalar* temp = new Scalar[n];

    for( int i=0; i<n; ++i)
      temp[i] = A.mfs_[i]->dot( *mfs_[i], false );

    MPI_Allreduce( temp, dots.data(), n, MPI_REAL8, MPI_SUM, comm() );
    delete[] temp;

  }


  /// \brief Compute the inner product for the \c MultiField considering it as one Vector.
  Scalar dot( const MV& A, bool global=true ) const {
    int n = getNumberVecs();

    Scalar b = 0.;

    for( int i=0; i<n; ++i )
      b+= mfs_[i]->dot( *A.mfs_[i], false );

    if( global ) this->reduceNorm( comm(), b );

    return( b );
  }


  /// \brief Compute the norm of each individual vector.
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  /// \todo implement OneNorm
  void norm(
      std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec,
      Belos::NormType type=Belos::TwoNorm ) const {

    const int n = getNumberVecs();
    Scalar* temp = new Scalar[n];

    switch(type) {
    case Belos::OneNorm:
      for( int i=0; i<n; ++i )
        temp[i] = mfs_[i]->norm(type,false);
      MPI_Allreduce( temp, normvec.data(), n, MPI_REAL8, MPI_SUM, comm() );
      break;
    case Belos::TwoNorm:
      for( int i=0; i<n; ++i )
        temp[i] = mfs_[i]->norm(type,false);
      MPI_Allreduce( temp, normvec.data(), n, MPI_REAL8, MPI_SUM, comm() );
      for( int i=0; i<n; ++i )
        normvec[i] = std::sqrt( normvec[i] );
      break;
    case Belos::InfNorm:
      for( int i=0; i<n; ++i )
        temp[i] = mfs_[i]->norm(type,false);
      MPI_Allreduce( temp, normvec.data(), n, MPI_REAL8, MPI_MAX, comm() );
      break;
    }
    delete[] temp;
  }


  /// \brief Compute the norm for the \c MultiField as it is considered as one Vector .
  /// \todo implement OneNorm
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    int n = getNumberVecs();
    Scalar normvec = 0.;

    for( int i=0; i<n; ++i ) {
      switch(type) {
          case Belos::OneNorm:
            normvec += mfs_[i]->norm(type,false);
            break;
          case Belos::TwoNorm:
            normvec += mfs_[i]->norm(type,false);
            break;
          case Belos::InfNorm:
            normvec = std::max( mfs_[i]->norm(type,false), normvec ) ;
            break;
          }
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

    double nor=0.;

    for( int i=0; i<getNumberVecs(); ++i )
      nor += mfs_[i]->norm( *weights.mfs_[i], false );

    if( global ) this->reduceNorm( comm(), nor, Belos::TwoNorm );

    return( nor );
  }

  /// \test this
  /// \param A
  /// \param index
  void SetBlock( const MV& A, const std::vector<int>& index ) {
    const int n = index.size();
    for( int i=0; i<n; ++i )
      mfs_[index[i]]->assign(*A.mfs_[i]);
  }

  /// \param A
  /// \param index
  void SetBlock( const MV& A, const Teuchos::Range1D& index) {
    for( int i=index.lbound(); i<=index.ubound(); ++i )
      mfs_[i]->assign(*A.mfs_[i-index.lbound()]);
  }

  /// \brief mv := A.
  ///
  /// assign (deep copy) A into mv.
  void assign( const MV& A ) {
    const int n = getNumberVecs();
    for( int i=0; i<n; ++i )
      mfs_[i]->assign( *A.mfs_[i] );
  }

  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed = false, int seed = 1) {
    const int n = getNumberVecs();
    for( int i=0; i<n; ++i )
      mfs_[i]->random();
  }

  /// \brief \f[ *this = \alpha \f]
  void init( Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    const int n = getNumberVecs();
    for( int i=0; i<n; ++i )
      mfs_[i]->init(alpha);
  }

  /// \todo add os to \c ScalarField and \c VectorField
  /// \param os
  void print( std::ostream& os ) {
    const int n = getNumberVecs();
    for( int i=0; i<n; ++i )
      mfs_[i]->print( os );
  }


  void write( int count=0 ) const {
    const int n = getNumberVecs();
    for( int i=0; i<n; ++i )
      mfs_[i]->write(count);
  }


  Field&       getField     (int i)       { return( *mfs_[i] ); }
  const Field& getConstField(int i) const { return( *mfs_[i] ); }

  Teuchos::RCP<Field>       getFieldPtr     (int i)       { return( mfs_[i] ); }
  Teuchos::RCP<const Field> getConstFieldPtr(int i) const { return( mfs_[i] ); }

  MPI_Comm comm() const { return( mfs_[0]->comm() ); }

}; // end of class MultiField


/// \brief factory for \c MultiField
/// \relates MultiField
template<class Field>
Teuchos::RCP< MultiField<Field> > createMultiField( const Field& field,
    const int numvecs, ECopyType ctype = ShallowCopy ) {

  return( Teuchos::rcp( new MultiField<Field>( field, numvecs, ctype ) ) );
}


/// \brief factory for \c MultiField.
///
/// simple wrapper.
/// \relates MultiField
template<class Field>
Teuchos::RCP< MultiField<Field> > createMultiField( const Teuchos::RCP<Field>& field ) {
  return( Teuchos::rcp( new MultiField<Field>( field ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIFIELD_HPP
