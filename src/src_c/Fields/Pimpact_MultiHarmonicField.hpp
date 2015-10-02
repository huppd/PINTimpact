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
//  Teuchos::RCP< MultiField< ModeField<Field> > > fields_;

  Teuchos::Array< Teuchos::RCP< ModeField<Field> > > fields_;

	Scalar* s_;

private:

	void allocate() {
		Ordinal n = field0_->getStorageSize();
		s_ = new Scalar[n*(1+2*space()->nGlo(3))];
		field0_->setStoragePtr( s_ );
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			fields_[i]->setStoragePtr( s_ + n + 2*n*i );
	}


public:

	Ordinal getStorageSize() const { return( (1+2*space()->nGlo(3))*field0_->getStorageSize() ); }

  MultiHarmonicField(
			const Teuchos::RCP<const SpaceT>& space ):
		AF( space ),
		field0_( Teuchos::rcp( new Field(space,false) ) ),
		fields_(space->nGlo(3)) {

			for( Ordinal i=0; i< space->nGlo(3); ++i )
				fields_[i] = Teuchos::rcp( new ModeField<Field>(space,false) );

			allocate();
			initField();

		};


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	MultiHarmonicField( const MultiHarmonicField& vF, ECopyType copyType=DeepCopy ):
		AF( vF.space() ),
		field0_( Teuchos::rcp( new Field( *vF.field0_, copyType ) ) ),
		fields_( vF.space()->nGlo(3) ) {

			for( Ordinal i=0; i< space()->nGlo(3); ++i )
				fields_[i] = Teuchos::rcp( new ModeField<Field>( vF.getConstField(i), copyType ) );

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

		};


	~MultiHarmonicField() { delete[] s_; }

  Teuchos::RCP<MV> clone( ECopyType ctype=DeepCopy ) const {

    return( Teuchos::rcp( new MV(*this,ctype) ) );

  }


  /// \name Attribute methods
  /// \{

  /// \warning it is assumed that both fields have the same \c StencilWidths


  Field&                    get0Field()               { return( *field0_ ); }
  const Field&              getConst0Field()    const { return( *field0_ ); }
  Teuchos::RCP<      Field> get0FieldPtr()            { return(  field0_ ); }
  Teuchos::RCP<const Field> getConst0FieldPtr() const { return(  field0_ ); }


  ModeField<Field>&                     getField        ( int i )       { return( *fields_[i] ); }
  const ModeField<Field>&               getConstField   ( int i ) const { return( *fields_[i] ); }
  Teuchos::RCP<      ModeField<Field> > getFieldPtr     ( int i )       { return( fields_[i]  ); }
  Teuchos::RCP<const ModeField<Field> > getConstFieldPtr( int i ) const { return( fields_[i]  ); }


  Field&                    getCField        ( int i )       { return( fields_[i]->getCField()      ); }
  const Field&              getConstCField   ( int i ) const { return( fields_[i]->getConstCField() ); }
  Teuchos::RCP<      Field> getCFieldPtr     ( int i )       { return( fields_[i]->getCFieldPtr()   ); }
  Teuchos::RCP<const Field> getConstCFieldPtr( int i ) const { return( fields_[i]->getConstCFieldPtr() ); }


  Field&                    getSField        ( int i )       { return( fields_[i]->getSField()         ); }
  const Field&              getConstSField   ( int i ) const { return( fields_[i]->getConstSField()    ); }
  Teuchos::RCP<      Field> getSFieldPtr     ( int i )       { return( fields_[i]->getSFieldPtr()      ); }
  Teuchos::RCP<const Field> getConstSFieldPtr( int i ) const { return( fields_[i]->getConstSFieldPtr() ); }


  Teuchos::RCP<const SpaceT> space() const { return( AF::space_ ); }

  const MPI_Comm& comm() const { return(field0_->comm()); }

  /// \brief returns the length of Field.
  ///
  /// the vector length is with regard to the inner points
  Ordinal getLength( bool nox_vec=false ) const {
		Ordinal len = field0_->getLength(true);
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			len += fields_[i]->getLength(true);
    return( len );
  }


  /// \brief get number of stored Field's
  /// \todo what makes sense here?
  int getNumberVecs() const { return( 1 ); }



  /// \}
  /// \name Update methods
  /// \{

  /// \brief Replace \c this with \f$\alpha A + \beta B\f$.
  void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {
    // add test for consistent VectorSpaces in debug mode
    field0_->add(alpha, *A.field0_, beta, *B.field0_);
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->add( alpha, A.getConstField(i), beta, B.getConstField(i) );
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void abs(const MV& y) {
    field0_->abs( *y.field0_ );
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->abs( y.getConstField(i) );
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal(const MV& y){
    field0_->reciprocal( *y.field0_ );
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->reciprocal( y.getConstField(i) );
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale( const Scalar& alpha ) {
    field0_->scale(alpha);
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->scale(alpha);
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void scale(const MV& a) {
    field0_->scale( *a.field0_ );
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->scale( a.getConstField(i) );
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  Scalar dot ( const MV& a, bool global=true ) const {

    Scalar b = 0.;

    b = field0_->dot( *a.field0_, false );
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			b += getConstFieldPtr(i)->dot( a.getConstField(i), false );

    if( global ) this->reduceNorm( comm(), b );

    return( b );

  }


  /// \}
  /// \name Norm method
  /// \{

  /// \brief Compute the norm of Field.
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = field0_->norm(type,false);

    switch(type) {
    case Belos::OneNorm:
			for( Ordinal i=0; i<space()->nGlo(3); ++i )
				normvec += getConstFieldPtr(i)->norm(type,false);
      break;
    case Belos::TwoNorm:
			for( Ordinal i=0; i<space()->nGlo(3); ++i )
				normvec += getConstFieldPtr(i)->norm(type,false);
      break;
    case Belos::InfNorm:
			for( Ordinal i=0; i<space()->nGlo(3); ++i )
				normvec = std::max( getConstFieldPtr(i)->norm(type,false), normvec );
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

    double normvec=field0_->norm(*weights.field0_,false);

		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			normvec += getConstFieldPtr(i)->norm(weights.getConstField(i),false);

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
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->assign( a.getConstField(i) );
  }


  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed = false, int seed = 1) {
    field0_->random();
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->random();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    field0_->init(alpha);
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->init(alpha);
  }

	void initField() {
    field0_->initField();
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->initField();
	}

	void setCornersZero() const {
    field0_->setCornersZero();
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getConstFieldPtr(i)->setCornersZero();
  }

  void level() const {
    field0_->level();
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getConstFieldPtr(i)->level();
  }

  /// \}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os )  {
    field0_->print( os );
		for( Ordinal i=0; i<space()->nGlo(3); ++i )
			getFieldPtr(i)->print( os );
  }

  void write( int count=0, bool time_evol_yes=false ) const {

		if( time_evol_yes ) {
			Scalar pi = 4.*std::atan(1.);
			Ordinal nf = space()->nGlo(3);
			Ordinal nt = 4*nf;
			auto temp = getConst0FieldPtr()->clone( Pimpact::ShallowCopy );
			for( Ordinal i=0; i<nt;  ++i ) {
				temp->assign( getConst0Field() );
//				temp->initField(  );
				for( Ordinal j=0; j<nf; ++j ) {
					temp->add(
							1., *temp,
							std::sin( 2.*pi*i*((Scalar)j+1.)/nt ), getConstSField(j) );
					temp->add(
							std::cos( 2.*pi*i*((Scalar)j+1.)/nt ), getConstCField(j),
							1., *temp );
					temp->write( count+i );
				}


			}

		}
		else{
			field0_->write(count);
			for( Ordinal i=0; i<space()->nGlo(3); ++i ) {
				getConstCFieldPtr(i)->write( count+2*i+1 );
				getConstSFieldPtr(i)->write( count+2*i+2 );
			}
		}
  }


}; // end of class MultiHarmonicField





/// \brief creates a multi-harmonic scalar field.
///
/// \relates MultiHarmonicField
/// \param space scalar Vector Space to which returned vector belongs
template<class FieldT>
Teuchos::RCP< MultiHarmonicField< FieldT > > createMultiHarmonic(
    const Teuchos::RCP<const typename FieldT::SpaceT >& space ) {

	return(
			create< MultiHarmonicField<FieldT> >( space ) );


}


} // end of namespace Pimpact


#ifdef COMPILE_ETI
#include "Pimpact_VectorField.hpp"
extern template class Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >;
#endif

#endif // end of #ifndef PIMPACT_MULTIHARMONICFIELD_HPP
