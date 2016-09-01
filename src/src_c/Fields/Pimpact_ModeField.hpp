#pragma once
#ifndef PIMPACT_MODEFIELD_HPP
#define PIMPACT_MODEFIELD_HPP


#include <vector>
#include <iostream>

#include "mpi.h"

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_AbstractField.hpp"




namespace Pimpact {



/// \brief important basic Vector class
/// vector for wrapping 2 fields into one mode
/// \ingroup Field
template<class IFT>
class ModeField : private AbstractField<typename IFT::SpaceT> {

public:

  using SpaceT = typename IFT::SpaceT;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal =typename SpaceT::Ordinal;

	using ScalarArray =  Scalar*;
	using FieldT = ModeField<IFT>;

	using AF =  AbstractField<SpaceT>;

	bool owning_;

	Teuchos::RCP<IFT> fieldc_;
  Teuchos::RCP<IFT> fields_;

	ScalarArray s_;

private:

	void allocate() {
		Ordinal n = fieldc_->getStorageSize();
		s_ = new Scalar[2*n];
		fieldc_->setStoragePtr( s_   );
		fields_->setStoragePtr( s_+n );
	}

public:

	constexpr Ordinal getStorageSize() const { return( 2*fieldc_->getStorageSize() ); }

	constexpr Scalar* getRawPtr() const { return( s_ ); }

  void setStoragePtr( Scalar*  array ) {
    s_ = array;
		fieldc_->setStoragePtr( s_                             );
		fields_->setStoragePtr( s_ + fieldc_->getStorageSize() );
  }


  ModeField( const Teuchos::RCP<const SpaceT>& space, bool owning=true ):
    AF( space ),
		owning_(owning),
    fieldc_( Teuchos::rcp( new IFT(space,false) ) ),
    fields_( Teuchos::rcp( new IFT(space,false) ) ) {

			if( owning_ ) {
				allocate();
				initField();
			}
	};

	~ModeField() { if( owning_ ) delete[] s_; }


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
  ModeField(const ModeField& vF, ECopyType copyType=DeepCopy):
    AF( vF.space() ),
		owning_( vF.owning_ ),
    fieldc_( Teuchos::rcp( new IFT(*vF.fieldc_,copyType) ) ),
		fields_( Teuchos::rcp( new IFT(*vF.fields_,copyType) ) ) {

			if( owning_ ) {

				allocate();

				switch( copyType ) {
					case ShallowCopy:
						initField();
						break;
					case DeepCopy:
						for( int i=0; i<getStorageSize(); ++i )
							s_[i] = vF.s_[i];
						break;
				}
			}
	};


  Teuchos::RCP<FieldT> clone( ECopyType cType=DeepCopy ) const {

		Teuchos::RCP<FieldT> mv = Teuchos::rcp( new FieldT( space() ) );

		switch( cType ) {
			case ShallowCopy:
				break;
			case DeepCopy:
				for( int i=0; i<getStorageSize(); ++i )
					mv->getRawPtr()[i] = this->s_[i];
				break;
		}

    return( mv );
  }

  /// \name Attribute methods
  /// \{

  constexpr Teuchos::RCP<IFT> getCFieldPtr() { return( fieldc_ ); }
  constexpr Teuchos::RCP<IFT> getSFieldPtr() { return( fields_ ); }

  constexpr Teuchos::RCP<const IFT> getConstCFieldPtr() const { return( fieldc_ ); }
  constexpr Teuchos::RCP<const IFT> getConstSFieldPtr() const { return( fields_ ); }

  constexpr IFT& getCField() { return( *fieldc_ ); }
  constexpr IFT& getSField() { return( *fields_ ); }

  constexpr const IFT& getConstCField() const { return( *fieldc_ ); }
  constexpr const IFT& getConstSField() const { return( *fields_ ); }

  constexpr const Teuchos::RCP<const SpaceT>& space() const { return( AF::space_ ); }

  constexpr const MPI_Comm& comm() const { return(fieldc_->comm()); }

  /// \brief returns the length of Field.
  ///
  /// should be the same as 2*fieldc_->getVecLength()
  /// the vector length is with regard to the inner points
  constexpr Ordinal getLength() const {
    return( fieldc_->getLength() + fields_->getLength() );
  }


  /// \brief get number of stored Field's
  constexpr int getNumberVecs() const { return( 1 ); }


  /// \}
  /// \name Update methods
  /// \{

  /// \brief Replace \c this with \f$\alpha A + \beta B\f$.
  void add( const Scalar& alpha, const FieldT& A, const Scalar& beta, const FieldT& B ) {
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
  void abs( const FieldT& y ) {
    fieldc_->abs( *y.fieldc_ );
    fields_->abs( *y.fields_ );
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal( const FieldT& y ) {
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
  void scale(const FieldT& a) {
    fieldc_->scale( *a.fieldc_ );
    fields_->scale( *a.fields_ );
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  constexpr Scalar dotLoc( const FieldT& a ) const {

    Scalar b=0.;

    b = fieldc_->dotLoc( *a.fieldc_) + fields_->dotLoc( *a.fields_ );

    return( b );
  }


	/// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
	constexpr Scalar dot( const FieldT& y ) const {

		return( this->reduce( comm(), dotLoc( y ) ) );
	}

  /// \}
  /// \name Norm method
  /// \{

  constexpr Scalar normLoc( Belos::NormType type=Belos::TwoNorm ) const {

    Scalar normvec = 
			(Belos::InfNorm==type)?
			std::max( fieldc_->normLoc(type), fields_->normLoc(type) ):
      ( fieldc_->normLoc(type) + fields_->normLoc(type) );

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
  constexpr Scalar normLoc(const FieldT& weights ) const {
		 return(
				 fieldc_->normLoc(*weights.fieldc_) +
				 fields_->normLoc(*weights.fields_)
				 );
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


  /// \}
  /// \name Initialization methods
  /// \{

  /// \brief mv := A
  /// Assign (deep copy) A into mv.
  void assign( const FieldT& a ) {
    fieldc_->assign( *a.fieldc_ );
    fields_->assign( *a.fields_ );
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

	void initField() {
    fieldc_->initField();
    fields_->initField();
	}

  void extrapolateBC() const {
    fieldc_->extrapolateBC();
    fields_->extrapolateBC();
  }

  void level() const {
    fieldc_->level();
    fields_->level();
  }

  /// \}
  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& out=std::cout ) const {
    fieldc_->print( out );
    fields_->print( out );
  }

  void write( int count=0 ) const {
    fieldc_->write(count);
    fields_->write(count+1);
  }

  /// \name comunication methods.
  /// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
  ///
  /// \{
	
  void exchange() const {
    fieldc_->exchange();
    fields_->exchange();
  }

  void setExchanged() const {
    fieldc_->setExchanged();
    fields_->setExchanged();
  }


  /// \}

}; // end of class ModeField


} // end of namespace Pimpact


#ifdef COMPILE_ETI
#include "Pimpact_VectorField.hpp"
extern template class Pimpact::ModeField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::ModeField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::ModeField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::ModeField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >;
#endif


#endif // end of #ifndef PIMPACT_MODEFIELD_HPP
