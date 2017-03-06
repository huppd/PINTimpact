#pragma once
#ifndef PIMPACT_MULTIFIELD_HPP
#define PIMPACT_MULTIFIELD_HPP


#include <vector>

#include "Teuchos_Array.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_AbstractField.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief templated class which is the interface to \c Belos and \c NOX
///
/// has multiple \c Field's, where \c Field can be a \c Pimpact:ScalarField, \c
/// Pimpact:VectorField or some combination with \c Pimpact::ModeField or \c
/// Pimpact::CompoundField
///
/// \note for better documentation, look at the equivalent documentation in the \c Belos::...
/// \note continous memory is not necessarily desirable because also views are used
/// \ingroup Field
template<class IFT>
class MultiField : private AbstractField<typename IFT::SpaceT> {

public:

  using InnerFieldT = IFT;

  using SpaceT = typename InnerFieldT::SpaceT;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;


  using FieldT = Pimpact::MultiField<InnerFieldT>;

  using AF = AbstractField<SpaceT>;

  Teuchos::Array<Teuchos::RCP<InnerFieldT> > mfs_;

public:

  /// \brief constructor taking a \c FieldT constructing multiple shallow copys.
  /// \note maybe hide and make it private
	MultiField( const InnerFieldT& field, const int& numvecs, ECopy ctyp=ECopy::Shallow ):
		AF( field.space() ), mfs_(numvecs) {
			for( int i=0; i<numvecs; ++i )
				mfs_[i] = field.clone(ctyp);
	}


  /// \brief cheap constructor from one InnerFieldT.
  ///
  /// creates simple wrapper from field(no coppying).
	MultiField( const Teuchos::RCP<InnerFieldT>& field ):
		AF( field->space() ), mfs_(1) {
			mfs_[0] = field;
	}


	MultiField( const FieldT& mv, ECopy ctype ):
		AF( mv.space() ), mfs_( mv.getNumberVecs() ) {

			for( int i=0; i<getNumberVecs(); ++i )
				mfs_[i] = mv.mfs_[i]->clone(ctype);
	}


	/// \brief  constructor, creates \c numvecs  empty fields
	///
	/// \param space
	/// \param numvecs
	/// \return
	MultiField( const Teuchos::RCP<const SpaceT>& space, const int& numvecs=1 ):
		AF( space ), mfs_(numvecs) {

			for( int i=0; i<numvecs; ++i )
				mfs_[i] = create<InnerFieldT>( space );
	}


  /// \brief Create a new \c MultiField with \c numvecs columns.
  ///
  /// The returned Pimpact::MultiField has the same (distribution over one or
  /// more parallel processes) Its entries are not initialized and have undefined
  /// values.
	Teuchos::RCP< FieldT > clone( const int& numvecs ) const {
		return( Teuchos::rcp( new FieldT( *mfs_[0], numvecs ) ) );
	}


  /// \brief Create a new \c MultiField with \c numvecs columns.
  ///
  /// The returned Pimpact::MultiField has the same (distribution over one or
  /// more parallel processes) Its entries are not initialized and have undefined
  /// values.
	Teuchos::RCP< FieldT > clone( ECopy ctype = ECopy::Deep ) const {
		return( Teuchos::rcp( new FieldT(*this,ctype) ) );
	}


  /// \return deep copy of \c MultiField
	Teuchos::RCP<FieldT> CloneCopy() const {
		return( clone() );
	}


  /// \brief deep copy of \c index fields, the new \c MultiFields stores \c index.size() vector.
  ///
  /// \param index
  /// \return
	Teuchos::RCP<FieldT> CloneCopy( const std::vector<int>& index ) const {
		Teuchos::RCP<FieldT> mv_ = Teuchos::rcp( new FieldT( space(), index.size() ) );
		for( unsigned int i=0; i<index.size(); ++i ) {
			mv_->mfs_[i] = mfs_[ index[i] ]->clone( ECopy::Deep );
		}
		return( mv_ );
	}


  /// deep copy of \c index fields, the new \c MultiFields stores \c index.size()
  /// \param index here index means an interval
  /// \return
  Teuchos::RCP<FieldT> CloneCopy( const Teuchos::Range1D& index) const {
    Teuchos::RCP<FieldT> mv_ = Teuchos::rcp( new FieldT(space(), index.size()) );
    int j = 0;
    for( int i=index.lbound(); i<=index.ubound(); ++i ) {
      mv_->mfs_[j] = mfs_[i]->clone( ECopy::Deep );
      ++j;
    }
    return( mv_ );
  }


  /// \param index
  /// \return nonConst View
  Teuchos::RCP<FieldT> CloneViewNonConst( const std::vector<int>& index) {
    Teuchos::RCP<FieldT> mv_ = Teuchos::rcp( new FieldT( space(), index.size() ) );
    for( unsigned int i=0; i<index.size(); ++i ) {
      mv_->mfs_[i] =  mfs_[ index[i] ];
    }
    return( mv_ );
  }


  Teuchos::RCP<FieldT> CloneViewNonConst( const Teuchos::Range1D& index ) {
    Teuchos::RCP<FieldT> mv_ = Teuchos::rcp( new FieldT( space(), index.size()) );
    int j=0;
    for( int i=index.lbound(); i<=index.ubound(); ++i ) {
      mv_->mfs_[j++] =  mfs_[i];
    }
    return( mv_ );
  }


  Teuchos::RCP<const FieldT > CloneView( const std::vector<int>& index ) const {
    Teuchos::RCP<FieldT> mv_ = Teuchos::rcp( new FieldT( space(), index.size()) );
    for( unsigned int i=0; i<index.size(); ++i ) {
      mv_->mfs_[i] =  mfs_[ index[i] ];
    }
    return( mv_ );
  }


  Teuchos::RCP<const FieldT > CloneView( const Teuchos::Range1D& index ) const {
    Teuchos::RCP<FieldT> mv_ = Teuchos::rcp( new FieldT( space(), index.size()) );
    int j=0;
    for( int i=index.lbound(); i<=index.ubound(); ++i ) {
      mv_->mfs_[j++] =  mfs_[i];
    }
    return( mv_ );
  }


  /// \brief returns the length of MultiField.
  constexpr Ordinal getLength() const {
		return( mfs_[0]->getLength() );
	}


  /// \brief get number of stored Field's
  constexpr int getNumberVecs() const { return( mfs_.size() ); }


  /// \brief is true
  constexpr bool HasConstantStride() const { return( true ); }


  /// \}
  /// \name Update methods
  /// \{

  /// \brief addes new field at end
  void push_back( const Teuchos::RCP<InnerFieldT>& field=Teuchos::null ) {
    if( Teuchos::is_null(field) )
      mfs_.push_back( mfs_.back()->clone(ECopy::Shallow) );
    else
      mfs_.push_back( field );
  }



  /// \brief \f[ *this = \alpha a b + \beta *this \f]
  ///
  /// \param alpha scalar
  /// \param a Vector
  /// \param b Matrix
  /// \param beta
	void TimesMatAdd( const Scalar& alpha, const FieldT& a,
			const Teuchos::SerialDenseMatrix<int,Scalar>& b,
			const Scalar& beta ) {

		int m1 = a.getNumberVecs(); ///< is assumed to be equal to number vecs of this and ncolumns and nrows of b
		int m2 = getNumberVecs();   ///< is assumed to be equal to number vecs of this and ncolumns and nrows of b

		//assert( m1==m2 ); // not necessary 

		for( int i=0; i<m2; ++i ) {
			mfs_[i]->add( beta, *mfs_[i], alpha*b(0,i), *a.mfs_[0] );
			for( int j=1; j<m1; ++j ) {
				mfs_[i]->add( 1., *mfs_[i], alpha*b(j,i), *a.mfs_[j] );
			}
		}
	}



  /// \brief <tt>mv := alpha*a + beta*b</tt>
	void add( Scalar alpha, const FieldT& a, Scalar beta, const FieldT& b, const B& wb=B::Y ) {
		for( int i=0; i<getNumberVecs(); ++i )
			mfs_[i]->add( alpha, *a.mfs_[i], beta, *b.mfs_[i], wb );
	}


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
	void abs( const FieldT& y) {
		for( int i=0; i<getNumberVecs(); ++i )
			mfs_[i]->abs( *y.mfs_[i] );
	}


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
	void reciprocal( const FieldT& y){
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
	void scale( const FieldT& a) {
		for( int i=0; i<getNumberVecs(); ++i )
			mfs_[i]->scale( *a.mfs_[i] );
	}


  /// \brief Scale each element of a \c Field by a \c gamma.
  ///
  /// Here x_j represents the j'th field, and we update it as
  /// \f[ x_j[i] = \alpha_j x_j[i] \quad \mbox{for } i=1,\dots,n \f]
	void scale( const std::vector<Scalar>& alphas ) {

		assert( alphas.size()==getNumberVecs() );

		for( unsigned int i=0; i<alphas.size(); ++i )
			mfs_[i]->scale( alphas[i] );
	}


  /// \f[ C_{ij} = \alpha A_i^T \c this_j \f]
  /// \param alpha
  /// \param A
  /// \param C
	/// \todo make dot local and reduce over C (cf dot)
	void Trans( Scalar alpha, const FieldT& A,
			Teuchos::SerialDenseMatrix<int,Scalar>& C) const {

		const int& n1 = getNumberVecs();
		const int& n2 = A.getNumberVecs();

		assert( n1==C.numRows() );
		assert( n2==C.numCols() );

		for( int i=0; i<n1; ++i )
			for( int j=0; j<n2; ++j )
				C(i,j) = alpha*mfs_[i]->dot( *A.mfs_[j] );
	}


  /// \}

  /// For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
	void dot( const FieldT& A, std::vector<Scalar>& dots) const {

		const int n = getNumberVecs();
		Scalar* temp = new Scalar[n];

		for( int i=0; i<n; ++i )
			temp[i] = A.mfs_[i]->dotLoc( *mfs_[i] );

		MPI_Allreduce( temp, dots.data(), n, MPI_REAL8, MPI_SUM, comm() );
		delete[] temp;
	}


  /// \brief Compute the inner product for the \c MultiField considering it as one Vector.
	constexpr Scalar dotLoc( const FieldT& A ) const {
		int n = getNumberVecs();

		Scalar b = 0.;

		for( int i=0; i<n; ++i )
			b+= mfs_[i]->dotLoc( *A.mfs_[i] );

		return( b );
	}

	/// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
	constexpr Scalar dot( const FieldT& y ) const {
		return( this->reduce( comm(), dotLoc( y ) ) );
	}

  /// \brief Compute the norm of each individual vector.
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
	void norm(
			std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec,
			Belos::NormType type=Belos::TwoNorm ) const {

		const int n = getNumberVecs();
		Scalar* temp = new Scalar[n];

		switch( type ) {
			case Belos::OneNorm:
				for( int i=0; i<n; ++i )
					temp[i] = mfs_[i]->normLoc(type);
				MPI_Allreduce( temp, normvec.data(), n, MPI_REAL8, MPI_SUM, comm() );
				break;
			case Belos::TwoNorm:
				for( int i=0; i<n; ++i )
					temp[i] = mfs_[i]->normLoc(type);
				MPI_Allreduce( temp, normvec.data(), n, MPI_REAL8, MPI_SUM, comm() );
				for( int i=0; i<n; ++i )
					normvec[i] = std::sqrt( normvec[i] );
				break;
			case Belos::InfNorm:
				for( int i=0; i<n; ++i )
					temp[i] = mfs_[i]->normLoc(type);
				MPI_Allreduce( temp, normvec.data(), n, MPI_REAL8, MPI_MAX, comm() );
				break;
		}
		delete[] temp;
	}


  /// \brief Compute the norm for the \c MultiField as it is considered as one Vector .
	constexpr Scalar normLoc(  Belos::NormType type = Belos::TwoNorm ) const {

		Scalar normvec = 0.;

		for( int i=0; i<getNumberVecs(); ++i )
			normvec = 
				(Belos::InfNorm==type)?
				(std::max( mfs_[i]->normLoc(type), normvec )):
				(normvec+mfs_[i]->normLoc(type));

		return( normvec );
	}

 /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  constexpr Scalar norm( Belos::NormType type = Belos::TwoNorm ) const {

		Scalar normvec = this->reduce(
				comm(),
				normLoc( type ),
				(Belos::InfNorm==type)?MPI_MAX:MPI_SUM );

		normvec =
			(Belos::TwoNorm==type) ?
				std::sqrt(normvec) :
				normvec;

    return( normvec );
  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
	constexpr Scalar normLoc( const FieldT& weights ) const {

		Scalar nor = Teuchos::ScalarTraits<Scalar>::zero();

		for( int i=0; i<getNumberVecs(); ++i )
			nor += mfs_[i]->normLoc( *weights.mfs_[i] );

		return( nor );
	}


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr Scalar norm( const FieldT& weights ) const {
		return( std::sqrt( this->reduce( comm(), normLoc( weights ) ) ) );
	}


  /// \test this
  /// \param A
  /// \param index
	void SetBlock( const FieldT& A, const std::vector<int>& index ) {
		const int n = index.size();
		for( int i=0; i<n; ++i )
			*mfs_[index[i]] = *A.mfs_[i];
	}

  /// \param A
  /// \param index
	void SetBlock( const FieldT& A, const Teuchos::Range1D& index) {
		for( int i=index.lbound(); i<=index.ubound(); ++i )
			*mfs_[i] = *A.mfs_[i-index.lbound()];
	}

  /// \brief *this := a.
  ///
  /// assign (deep copy) a into *this.
	MultiField& operator=( const MultiField& a ) {

		const int n = getNumberVecs();
		for( int i=0; i<n; ++i )
			*mfs_[i] = *a.mfs_[i];

		return *this;
	}

  /// \brief Replace the vectors with a random vectors.
	void random(bool useSeed = false, int seed = 1) {
		const int n = getNumberVecs();
		for( int i=0; i<n; ++i )
			mfs_[i]->random();
	}

  /// \brief \f[ *this = \alpha \f]
	void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero(), const B& wB=B::Y ) {
		const int n = getNumberVecs();
		for( int i=0; i<n; ++i )
			mfs_[i]->init(alpha,wB);
	}

	void extrapolateBC( const Belos::ETrans& trans=Belos::NOTRANS ) {
		const int n = getNumberVecs();
		for( int i=0; i<n; ++i )
			mfs_[i]->extrapolateBC( trans );
	}

	void level() const {
		const int n = getNumberVecs();
		for( int i=0; i<n; ++i )
			mfs_[i]->level();
	}

  /// \param out out stream 
  void print( std::ostream& out=std::cout ) const {
		const int n = getNumberVecs();
		for( int i=0; i<n; ++i )
			mfs_[i]->print( out );
	}


	void write( int count=0 ) const {
		const int n = getNumberVecs();
		for( int i=0; i<n; ++i )
			mfs_[i]->write(count + i);
	}

	/// \name Attribute methods
	/// @{ 

	InnerFieldT& getField( const int& i ) { return( *mfs_[i] ); }
	constexpr const InnerFieldT& getField( const int& i ) { return( *mfs_[i] ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return( AF::space_ ); }

	constexpr const MPI_Comm& comm() const { return(mfs_[0]->comm()); }

	///  @} 

}; // end of class MultiField



/// \brief factory for \c MultiField.
///
/// simple wrapper.
/// \relates MultiField
template<class FieldT>
Teuchos::RCP< MultiField<FieldT> > wrapMultiField( const Teuchos::RCP<FieldT>& field ) {

	return( Teuchos::rcp( new MultiField<FieldT>( field ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIFIELD_HPP
