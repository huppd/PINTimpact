#pragma once
#ifndef PIMPACT_MULTIFIELD_HPP
#define PIMPACT_MULTIFIELD_HPP


/*! \file Pimpact_MultiField.hpp
    \brief Provides several interfaces between Belos virtual classes and Tpetra concrete classes.
*/

#include <vector>
//#include <Teuchos_Assert.hpp>
//#include <Teuchos_ScalarTraits.hpp>
//#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Array.hpp>
#include "Teuchos_RCP.hpp"
#include <Teuchos_Range1D.hpp>
//#include "Teuchos_LocalTestingHelpers.hpp"

//#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
//#include <BelosMultiVecTraits.hpp>
//#include <BelosOperatorTraits.hpp>

//#ifdef HAVE_BELOS_TSQR
//#  include <Tpetra_TsqrAdaptor.hpp>
//#endif // HAVE_BELOS_TSQR

#include "Pimpact_Types.hpp"
//#include "Pimpact_ModeField.hpp"
//#include "Pimpact_VectorField.hpp"
//#include "Pimpact_ScalarField.hpp"


namespace Pimpact {


/**
 * \brief templated class which is the interface to \c Belos and \c NOX
 *
 * has multiple \c Field's, where \c Field can be a \c Pimpact:ScalarField, \c
 * Pimpact:VectorField or some combination with \c Pimpact::ModeField or \c
 * Pimpact::CompoundField \note if this is heavily used for many Field's, then
 * the implementation should be improved such that communication is done such
 * that only done once per MV not per Field
 *
 * \note for better documentation, look at the equivalent documentation in the \c Belos::...
 */
template<class Field>
class MultiField {

public:

	typedef typename Field::Scalar Scalar;
	typedef typename Field::Ordinal Ordinal;

private:

	typedef Pimpact::MultiField<Field> MV;
	Teuchos::Array<Teuchos::RCP<Field> > mfs_;

public:
	MultiField():mfs_(Teuchos::null) {};

	/// \brief constructor taking a \c Field constructing multiple shallow copys
  /// \note maybe hide and make it private
  MultiField( const Field& field, const int numvecs, Pimpact::ECopyType ctyp = Pimpact::ShallowCopy):mfs_(numvecs) {
  	for( int i=0; i<numvecs; ++i )
  		mfs_[i] = Teuchos::rcp( new Field(field, ctyp ) );
  }


  /**
   * \brief copy constructor creating a view
   */
  MultiField( const MV& mv ):mfs_(mv.mfs_) {}


	/**
   * \brief  constructor, creates \c numvecs  empty Fields
   * @param numvecs
   * @return
  */
  MultiField( int numvecs ):mfs_(numvecs) {}


  /** \brief Create a new \c MultiField with \c numvecs columns.
   *
   * The returned Pimpact::MultiField has the same (distribution over one or
   * more parallel processes) Its entries are not initialized and have undefined
   * values.
  */
  Teuchos::RCP< MV > clone( const int numvecs ) const {
  	return( Teuchos::rcp( new MV( *mfs_[0], numvecs ) ) );
  }


  /** \brief Create a new \c MultiField with \c numvecs columns.
   *
   * The returned Pimpact::MultiField has the same (distribution over one or
   * more parallel processes) Its entries are not initialized and have undefined
   * values.
  */
  Teuchos::RCP< MV > clone( ECopyType ctype = DeepCopy ) const {
    auto mv_ = Teuchos::rcp( new MV(getNumberVecs()) );
      for( int i=0; i<getNumberVecs(); ++i ) {
//        mv_->mfs_[i] = Teuchos::rcp( new Field( *mfs_[i], type ) );
        mv_->mfs_[i] = mfs_[i]->clone(ctype);
      }
      return( mv_ );
  }


  /**
   * @return deep copy of \c MultiField
   */
  Teuchos::RCP<MV> CloneCopy() const {
    return( clone() );
  }


  /**
   * deep copy of \c index fields, the new \c MultiFields stores \c index.size() vector
   * @param index
   * @return
  */
  Teuchos::RCP<MV>  CloneCopy( const std::vector<int>& index ) const {
  	auto mv_ = Teuchos::rcp( new MV(index.size()) );
  	for( unsigned int i=0; i<index.size(); ++i ) {
  		mv_->mfs_[i] = Teuchos::rcp( new Field( *mfs_[ index[i] ], DeepCopy ) );
  	}
  	return( mv_ );
  }


  /**
   * deep copy of \c index fields, the new \c MultiFields stores \c index.size()
   * @param index here index means an interval
   * @return
  */
  Teuchos::RCP<MV> CloneCopy( const Teuchos::Range1D& index) const {
  	auto mv_ = Teuchos::rcp( new MV(index.size()) );
  	int j = 0;
  	for( int i=index.lbound(); i<=index.ubound(); ++i ) {
  		mv_->mfs_[j] = Teuchos::rcp( new Field( *mfs_[ i ], DeepCopy ) );
  		++j;
  	}
  	return( mv_ );
  }


  /**
   *
   * @param index
   * @return nonConst View
  */
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


  /** \brief returns the length of Field.
   * \param nox_vec if \c MultiField is used for NOX the Vector length is
   * considered for all Fields
   */
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


  /**
   * \brief \f[ *this = \alpha A B + \beta *this \f]
   * \todo add debugtests
   * @param alpha
   * @param A
   * @param B
   * @param beta
   */
  void TimesMatAdd( const Scalar& alpha, const MV& A,
  		const Teuchos::SerialDenseMatrix<int,Scalar>& B,
  		const Scalar& beta ) {
  	auto AB = CloneCopy(); /// costly should be avoided
  	int m1 = A.getNumberVecs(); ///< is assumed to be equal to number vecs of this and ncolumns and nrows of b
  	int m2 = getNumberVecs(); ///< is assumed to be equal to number vecs of this and ncolumns and nrows of b

  	AB->init(0.);
  	for( int j=0; j<m1; ++j ) {
  		for( int i=0; i<m2; ++i ) {
  			AB->mfs_[i]->add(1., *(AB->mfs_[i]), B(j,i), *A.mfs_[j] );
  		}
  	}
  	for( int i=0; i<m2; ++i ) {
  		mfs_[i]->add( alpha, *AB->mfs_[i], beta, *mfs_[i] );
  	}
  }


  /**
   * \brief <tt>mv := alpha*A + beta*B</tt>
   *
   *	The Tpetra specialization of this method ignores and completely
   *	overwrites any NaN or Inf entries in A.  Thus, it does <i>not</i> mean
   *	the same thing as <tt>mv := 0*mv + alpha*A + beta*B</tt> in IEEE 754
   *	floating-point arithmetic. (Remember that NaN*0 = NaN.)
   */
  void add( Scalar alpha, const MV& A, Scalar beta, const MV& B ) {
  	for( int i=0; i<A.getNumberVecs(); ++i )
  		mfs_[i]->add( alpha, *A.mfs_[i], beta, *B.mfs_[i] );
  }


  /**
   * \brief Scale each element of every \c Field by \c gamma.
   *
   * Here x represents on \c Field, and we update it as
   * \f[ x_i = \alpha x_i \quad \mbox{for } i=1,\dots,n \f]
   */
  void scale( const Scalar& alpha ) {
  	for( int i=0; i<getNumberVecs(); ++i )
  		mfs_[i]->scale(alpha);
  }


  /**
   * \brief Scale each element of a \c Field by a \c gamma.
   *
   * Here x_j represents the j'th field, and we update it as
   * \f[ x_j[i] = \alpha_j x_j[i] \quad \mbox{for } i=1,\dots,n \f]
   */
  void scale( const std::vector<Scalar>& alphas ) {
#ifdef DEBUG
  	TEST_EQUALITY( alphas.size(), getNumberVecs() );
#endif
  	for( unsigned int i=0; i<alphas.size(); ++i )
  		mfs_[i]->scale( alphas[i] );
  }


  /**
   * \f[ C_{ij} = \alpha A_i^T \c this_j \f]
   * @param alpha
   * @param A
   * @param C
  */
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


  /// For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
  void dot( const MV& A, std::vector<Scalar>& dots) const {
  	const int n = getNumberVecs();
  	for( int i=0; i<n; ++i)
  		dots[i] = A.mfs_[i]->dot( *mfs_[i] );
  }


  /**
   * \brief Compute the inner product for the \c MultiField considering it as one Vector.
   * \todo test this
   */
	Scalar dot( const MV& A ) const {
		int n = getNumberVecs();
		Scalar innpro=0.;
		std::vector<Scalar> innprovec( n );
		dot( A, innprovec );

		for( int i=0; i<n; ++i ) {
			innpro += innprovec[i];
		}
		return( innpro );
	}


  /**
   * \brief Compute the norm of each individual vector.
   * Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
   * \todo implement OneNorm
   */
  void norm( std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec,
  		Belos::NormType type=Belos::TwoNorm) const {
  	const int n = getNumberVecs();
  	for( int i=0; i<n; ++i )
  		normvec[i] = mfs_[i]->norm(type);
  }


  /**
   * \brief Compute the norm for the \c MultiField as it is considered as one Vector .
   * \todo implement OneNorm
   */
	Scalar norm(  Belos::NormType type = Belos::TwoNorm ) const {
		int n = getNumberVecs();
		Scalar nor=0.;
		std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> normvec( n );
		norm( normvec, type );

		for( int i=0; i<n; ++i ) {
			nor += normvec[i];
		}
		return( nor );
	}


  /**
   * \test this
   * @param A
   * @param index
  */
  void SetBlock( const MV& A, const std::vector<int>& index ) {
  	const int n = index.size();
  	for( int i=0; i<n; ++i )
  		mfs_[index[i]]->assign(*A.mfs_[i]);
  }


  /**
   * @param A
   * @param index
   */
  void SetBlock( const MV& A, const Teuchos::Range1D& index) {
  	for( int i=index.lbound(); i<=index.ubound(); ++i )
  		mfs_[i]->assign(*A.mfs_[i-index.lbound()]);
  }

  /**
    * \brief mv := A
    * assign (deep copy) A into mv.
    */
  void assign( const MV& A ) {
  	const int n = getNumberVecs();
  	for( int i=0; i<n; ++i )
  		mfs_[i]->assign( *A.mfs_[i] );
  }


  /**
   * \brief Replace the vectors with a random vectors.
   */
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


  /**
   * \todo add os to \c ScalarField and \c VectorField
   * @param os
  */
  void print( std::ostream& os ) {
  	const int n = getNumberVecs();
  	for( int i=0; i<n; ++i )
  		mfs_[i]->print( os );
  }


  void write( int count=0 ) {
  	const int n = getNumberVecs();
  	for( int i=0; i<n; ++i )
  		mfs_[i]->write(count);
  }


  Field& GetVec(int i) { return( *mfs_[i] ); }


  const Field& GetConstVec(int i) const { return( *mfs_[i] ); }

}; // end of class MultiField


template<class Field, class Scalar, class Ordinal>
Teuchos::RCP< MultiField<Field> > createMultiField( const Field& field,
		const int numvecs, ECopyType ctype = ShallowCopy ) {

	return( Teuchos::rcp( new MultiField<Field>( field, numvecs, ctype ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MULTIFIELD_HPP
