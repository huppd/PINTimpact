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
	 * has multiple \c Field types, where \c Field can be a \c Pimpact:ScalarField or \c Pimpact:VectorField
	 * \note if this is heavily used for many Field's, then the implementation should be improved such that communication
   * is done such that only done once per MV not per Field
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
//#ifdef HAVE_BELOS_TPETRA_TIMERS
//    static Teuchos::RCP<Teuchos::Time> mvTimesMatAddMvTimer_, mvTransMvTimer_;
//#endif

  	/// \brief constructor taking a \c Field constructing multiple shallow copys
  	/// \todo maybe hide and make it private
  	MultiField( const Field& field, const int numvecs, Pimpact::ECopyType ctyp = Pimpact::ShallowCopy):mfs_(numvecs) {
  		for( int i=0; i<numvecs; ++i )
  			mfs_[i] = Teuchos::rcp( new Field(field, ctyp ) );
  	}


  	/**
  	 * copy constructor creating a view
  	 */
  	MultiField( const MV& mv ):mfs_(mv.mfs_) {}


  	/**
  	 * constructor, creates \c numvecs  empty Fields
  	 * @param numvecs
  	 * @return
  	 */
  	MultiField( int numvecs ):mfs_(numvecs) {}


    /// \brief Create a new multivector with \c numvecs columns.
    ///
    /// The returned Pimpact::MultiVector has the same
    /// (distribution over one or more parallel processes)
    /// Its entries are not initialized and have undefined values.
    Teuchos::RCP< MV >
    Clone ( const int numvecs) const {
      return Teuchos::rcp( new MV( *mfs_[0], numvecs ) );
    }


    /**
     * @return deep copy of \c MultiField
     */
    Teuchos::RCP<MV> CloneCopy() const {
    	auto mv_ = Teuchos::rcp( new MV(GetNumberVecs()) );
  		for( int i=0; i<GetNumberVecs(); ++i ) {
  			mv_->mfs_[i] = Teuchos::rcp( new Field( *mfs_[i], DeepCopy ) );
  		}
      return mv_;
    }


    /**
     * deep copy of \c index fields, the new \c MultiFields stores \c index.size() vectors
     * @param index
     * @return
     */
    Teuchos::RCP<MV>  CloneCopy( const std::vector<int>& index ) const {
    	auto mv_ = Teuchos::rcp( new MV(index.size()) );
  		for( unsigned int i=0; i<index.size(); ++i ) {
  			mv_->mfs_[i] = Teuchos::rcp( new Field( *mfs_[ index[i] ], DeepCopy ) );
  		}
      return mv_;
    }


    /**
     * deep copy of \c index fields, the new \c MultiFields stores \c index.size()
     * @param index here index means an interval
     * @return
     */
    Teuchos::RCP<MV>
    CloneCopy( const Teuchos::Range1D& index) const {
    	auto mv_ = Teuchos::rcp( new MV(index.size()) );
    	int j = 0;
  		for( int i=index.lbound(); i<=index.ubound(); ++i ) {
  			mv_->mfs_[j] = Teuchos::rcp( new Field( *mfs_[ i ], DeepCopy ) );
  			++j;
  		}
      return mv_;
    }


    /**
     *
     * @param index
     * @return nonConst View
     */
    Teuchos::RCP<MV>
    CloneViewNonConst( const std::vector<int>& index) {
    	auto mv_ = Teuchos::rcp( new MV( index.size() ) );
  		for( unsigned int i=0; i<index.size(); ++i ) {
  			mv_->mfs_[i] =  mfs_[ index[i] ];
  		}
      return mv_;
    }


    Teuchos::RCP<MV>
    CloneViewNonConst( const Teuchos::Range1D& index ) {
    	auto mv_ = Teuchos::rcp( new MV(index.size()) );
    	int j=0;
  		for( int i=index.lbound(); i<=index.ubound(); ++i ) {
  			mv_->mfs_[j++] =  mfs_[i];
  		}
      return mv_;
    }


    Teuchos::RCP<const MV >
    CloneView( const std::vector<int>& index ) const {
    	auto mv_ = Teuchos::rcp( new MV(index.size()) );
  		for( unsigned int i=0; i<index.size(); ++i ) {
  			mv_->mfs_[i] =  mfs_[ index[i] ];
  		}
      return mv_;
    }

    Teuchos::RCP<const MV >
    CloneView( const Teuchos::Range1D& index ) const {
    	auto mv_ = Teuchos::rcp( new MV(index.size()) );
    	int j=0;
  		for( int i=index.lbound(); i<=index.ubound(); ++i ) {
  			mv_->mfs_[j++] =  mfs_[i];
  		}
      return mv_;
    }

    int GetVecLength() const
    { return mfs_[0]->getVecLength(); }

    int GetNumberVecs() const
    { return mfs_.size(); }

    bool HasConstantStride() const
    { return true; }

    /**
     * \todo add debugtests
     * \warning alot can go wrong here
     * @param alpha
     * @param A
     * @param B
     * @param beta
     */
    void
    TimesMatAdd( const Scalar& alpha,
                     const MV& A,
                     const Teuchos::SerialDenseMatrix<int,Scalar>& B,
                     const Scalar& beta ) {
    	auto AB = CloneCopy(); /// costly should be avoided
    	int m1 = A.GetNumberVecs(); ///< is assumed to be equal to number vecs of this and ncolumns and nrows of b
    	int m2 = GetNumberVecs(); ///< is assumed to be equal to number vecs of this and ncolumns and nrows of b

//    	std::cout << "A.GetNumberVecs(): " <<  A.GetNumberVecs() << ", GetNumberVecs(): " << GetNumberVecs() << ", B.numRows(): " << B.numRows()<< ", B.numCols(): " << B.numCols() << "\n";
    	AB->Init(0.);
    	for( int j=0; j<m1; ++j ) {
				for( int i=0; i<m2; ++i ) {
					AB->mfs_[i]->add(1., *(AB->mfs_[i]), B(j,i), *A.mfs_[j] );
				}
    	}
    	for( int i=0; i<m2; ++i ) {
    		mfs_[i]->add( alpha, *AB->mfs_[i], beta, *mfs_[i] );
    	}
    }

    /// \brief <tt>mv := alpha*A + beta*B</tt>
    ///
    /// The Tpetra specialization of this method ignores and
    /// completely overwrites any NaN or Inf entries in A.  Thus, it
    /// does <i>not</i> mean the same thing as <tt>mv := 0*mv +
    /// alpha*A + beta*B</tt> in IEEE 754 floating-point arithmetic.
    /// (Remember that NaN*0 = NaN.)
    void
    Add( Scalar alpha,
             const MV& A,
             Scalar beta,
             const MV& B ) {
    	for( int i=0; i<A.GetNumberVecs(); ++i )
    		mfs_[i]->add( alpha, *A.mfs_[i], beta, *B.mfs_[i] );
    }

    void
    Scale( const Scalar& alpha ) {
    	for( int i=0; i<GetNumberVecs(); ++i )
				mfs_[i]->scale(alpha);
    }

    void
    Scale( const std::vector<Scalar>& alphas ) {
#ifdef DEBUG
    	TEST_EQUALITY( alphas.size(), GetNumberVecs() );
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
    void
    Trans( Scalar alpha,
               const MV& A,
//               const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
               Teuchos::SerialDenseMatrix<int,Scalar>& C) const {

    	const int n1 = GetNumberVecs();
    	const int n2 = A.GetNumberVecs();
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

    //! For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
    void
    Dot( const MV& A,
//           const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
           std::vector<Scalar>& dots) const {

    	const int n = GetNumberVecs();
			for( int i=0; i<n; ++i)
				dots[i] = A.mfs_[i]->dot( *mfs_[i] );
    }

    //! For all columns j of mv, set <tt>normvec[j] = norm(mv[j])</tt>.
    void Norm( std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec,
            Belos::NormType type=Belos::TwoNorm) const {
    	const int n = GetNumberVecs();
			for( int i=0; i<n; ++i )
				normvec[i] = mfs_[i]->norm(type);
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
     * \test this
     * @param A
     * @param index
     */
    void
    SetBlock( const MV& A,
              const Teuchos::Range1D& index) {
			for( int i=index.lbound(); i<=index.ubound(); ++i )
				mfs_[i]->assign(*A.mfs_[i-index.lbound()]);
    }

    void
    Assign( const MV& A ) {
    	const int n = GetNumberVecs();
			for( int i=0; i<n; ++i )
				mfs_[i]->assign( *A.mfs_[i] );
    }


    void Random() {
    	const int n = GetNumberVecs();
			for( int i=0; i<n; ++i )
				mfs_[i]->random();
    }

    void Init( Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    	const int n = GetNumberVecs();
			for( int i=0; i<n; ++i )
				mfs_[i]->init(alpha);
    }

    /**
     * \todo add os to \c ScalarField and \c VectorField
     * @param os
     */
    void Print( std::ostream& os )
    {
    	const int n = GetNumberVecs();
			for( int i=0; i<n; ++i )
				mfs_[i]->print();
    }

    void write( int count=0 ) {
    	const int n = GetNumberVecs();
			for( int i=0; i<n; ++i )
				mfs_[i]->write(count);
    }

    Field& GetVec(int i) { return *mfs_[i]; }
    const Field& GetConstVec(int i) const { return *mfs_[i]; }

  };


template<class Field, class Scalar, class Ordinal>
Teuchos::RCP< MultiField<Field> > createMultiField(const Field& field, const int numvecs, Pimpact::ECopyType ctype = Pimpact::ShallowCopy ) {
	return Teuchos::rcp( new MultiField<Field>( field, numvecs, ctype ) );
//         Teuchos::rcp( new Pimpact::MultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10) );
//	Teuchos::RCP<Pimpact::MultiField<Pimpact::ScalarField<double,int>,double,int> > mv =

}

} // end of Pimpact namespace

#endif // end of file PIMPACT_MULTIFIELD_HPP
