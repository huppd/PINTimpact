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
/// \todo continous memory
template<class FieldT>
class ModeField : private AbstractField<typename FieldT::SpaceT> {

public:

  using SpaceT = typename FieldT::SpaceT;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal =typename SpaceT::Ordinal;

  static const int dimension = SpaceT::dimension;


protected:

	using ScalarArray =  Scalar*;
	using MV = ModeField<FieldT>;

	using AF =  AbstractField< typename FieldT::SpaceT>;

	bool owning_;

	Teuchos::RCP<FieldT> fieldc_;
  Teuchos::RCP<FieldT> fields_;

	ScalarArray s_;

private:

	void allocate() {
		Ordinal n = fieldc_->getStorageSize();
		s_ = new Scalar[2*n];
		fieldc_->setStoragePtr( s_   );
		fields_->setStoragePtr( s_+n );
	}

public:

	Ordinal getStorageSize() const { return( 2*fieldc_->getStorageSize() ); }

	Scalar* getRawPtr() const { return( s_ ); }

  void setStoragePtr( Scalar*  array ) {
    s_ = array;
		fieldc_->setStoragePtr( s_                             );
		fields_->setStoragePtr( s_ + fieldc_->getStorageSize() );
  }


  ModeField( const Teuchos::RCP<const SpaceT>& space, bool owning=true ):
    AF( space ),
		owning_(owning),
    fieldc_( Teuchos::rcp( new FieldT(space,false) ) ),
    fields_( Teuchos::rcp( new FieldT(space,false) ) ) {

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
    fieldc_( Teuchos::rcp( new FieldT(*vF.fieldc_,copyType) ) ),
		fields_( Teuchos::rcp( new FieldT(*vF.fields_,copyType) ) ) {

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


  Teuchos::RCP<MV> clone( ECopyType cType=DeepCopy ) const {

		Teuchos::RCP<MV> mv = Teuchos::rcp( new MV( space() ) );

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

  Teuchos::RCP<FieldT> getCFieldPtr() { return( fieldc_ ); }
  Teuchos::RCP<FieldT> getSFieldPtr() { return( fields_ ); }

  Teuchos::RCP<const FieldT> getConstCFieldPtr() const { return( fieldc_ ); }
  Teuchos::RCP<const FieldT> getConstSFieldPtr() const { return( fields_ ); }

  FieldT& getCField() { return( *fieldc_ ); }
  FieldT& getSField() { return( *fields_ ); }

  const FieldT& getConstCField() const { return( *fieldc_ ); }
  const FieldT& getConstSField() const { return( *fields_ ); }

  Teuchos::RCP<const SpaceT> space() const { return( AF::space_ ); }

  const MPI_Comm& comm() const { return(fieldc_->comm()); }

  /// \brief returns the length of Field.
  ///
  /// should be the same as 2*fieldc_->getVecLength()
  /// the vector length is with regard to the inner points
  Ordinal getLength( bool nox_vec=false ) const {
    return( fieldc_->getLength(nox_vec) + fields_->getLength(nox_vec) );
  }


  /// \brief get number of stored Field's
  int getNumberVecs() const { return( 1 ); }


  /// \}
  /// \name Update methods
  /// \{

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
  void abs( const MV& y ) {
    fieldc_->abs( *y.fieldc_ );
    fields_->abs( *y.fields_ );
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal( const MV& y ) {
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
  void scale(const MV& a) {
    fieldc_->scale( *a.fieldc_ );
    fields_->scale( *a.fields_ );
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  Scalar dot ( const MV& a, bool global=true ) const {

    Scalar b=0.;

    b = fieldc_->dot( *a.fieldc_, false ) + fields_->dot( *a.fields_, false );

    if( global ) this->reduceNorm( comm(), b );

    return( b );

  }


  /// \}
  /// \name Norm method
  /// \{

  /// \brief Compute the norm of Field.
  ///
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  Scalar norm(  Belos::NormType type=Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = 0;

    switch(type) {
    default:
      normvec = fieldc_->norm(type,false) + fields_->norm(type,false);
      break;
    case Belos::InfNorm:
      normvec = std::max( fieldc_->norm(type,false), fields_->norm(type,false) ) ;
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
  double norm(const MV& weights, bool global=true) const {

    double normvec=fieldc_->norm(*weights.fieldc_,false)+fields_->norm(*weights.fields_,false);

    if( global ) this->reduceNorm( comm(), normvec, Belos::TwoNorm );

    return( normvec );

  }


  /// \}
  /// \name Initialization methods
  /// \{

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

	void initField() {
    fieldc_->initField();
    fields_->initField();
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
	
//  void changed() const {
//    fieldc_->changed( dir );
//    fields_->changed( dir );
//  }

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
