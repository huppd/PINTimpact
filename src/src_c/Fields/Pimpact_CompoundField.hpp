#pragma once
#ifndef PIMPACT_COMPOUNDFIELD_HPP
#define PIMPACT_COMPOUNDFIELD_HPP


#include <vector>
#include <iostream>

#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_AbstractField.hpp"




namespace Pimpact {


/// \brief important basic Vector class
/// vector for wrapping 2 fields into one mode
/// \ingroup Field
/// \todo use attributes methods in vectorspace functions
template<class VField, class SField>
class CompoundField : private AbstractField<typename VField::SpaceT> {

public:

  using SpaceT = typename VField::SpaceT;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal =typename SpaceT::Ordinal;

  using FieldT = CompoundField<VField,SField>;

  using AF = AbstractField<SpaceT>;

  Teuchos::RCP<VField> vfield_;
  Teuchos::RCP<SField> sfield_;

public:

  CompoundField( const Teuchos::RCP<const SpaceT>& space, EField dummy=S ):
        AF( space ),
        vfield_( create<VField>(space) ),
        sfield_( create<SField>(space) ) {};

  CompoundField(
      const Teuchos::RCP<VField>& vfield,
      const Teuchos::RCP<SField>& sfield ):
        AF( vfield->space() ),
        vfield_(vfield),
        sfield_(sfield) {};


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param field 
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
  CompoundField( const CompoundField& field, ECopyType copyType=DeepCopy ):
    AF( field.space() ),
    vfield_( Teuchos::rcp( new VField( *field.vfield_, copyType ) ) ),
    sfield_( Teuchos::rcp( new SField( *field.sfield_, copyType ) ) )
	{};


  Teuchos::RCP<FieldT> clone( ECopyType ctype=DeepCopy ) const {
    return( Teuchos::rcp( new FieldT( *this, ctype ) ) );
  }

  /// \name Attribute methods
  /// \{

  constexpr VField& getVField() { return( *vfield_ ); }
  constexpr SField& getSField() { return( *sfield_ ); }

  constexpr const VField& getConstVField() const { return( *vfield_ ); }
  constexpr const SField& getConstSField() const { return( *sfield_ ); }

  constexpr Teuchos::RCP<VField> getVFieldPtr() { return( vfield_ ); }
  constexpr Teuchos::RCP<SField> getSFieldPtr() { return( sfield_ ); }

  constexpr Teuchos::RCP<const VField> getConstVFieldPtr() const { return( vfield_ ); }
  constexpr Teuchos::RCP<const SField> getConstSFieldPtr() const { return( sfield_ ); }

  constexpr const Teuchos::RCP<const SpaceT>& space() const { return( AF::space_ ); }

  constexpr const MPI_Comm& comm() const { return(vfield_->comm()); }


  /// \brief get Vect length
  /// shoud be the same as 2*vfield_->getVecLength()
  /// the vector length is withregard to the inner points
  /// \return vector length
  /// \brief returns the length of Field.
  constexpr Ordinal getLength() const {
    return( vfield_->getLength() + sfield_->getLength() );
  }


  /// \brief get number of stored Field's
  constexpr int getNumberVecs() const { return( 1 ); }


  /// \}
  /// \name Update methods
  /// \{

  /// \brief Replace \c this with \f$\alpha A + \beta B\f$.
  void add( const Scalar& alpha, const FieldT& A, const Scalar& beta, const FieldT& B ) {
    // add test for consistent VectorSpaces in debug mode
    vfield_->add(alpha, *A.vfield_, beta, *B.vfield_);
    sfield_->add(alpha, *A.sfield_, beta, *B.sfield_);
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void abs( const FieldT& y ) {
    vfield_->abs( *y.vfield_ );
    sfield_->abs( *y.sfield_ );
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal( const FieldT& y){
    vfield_->reciprocal( *y.vfield_ );
    sfield_->reciprocal( *y.sfield_ );
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale( const Scalar& alpha ) {
    vfield_->scale(alpha);
    sfield_->scale(alpha);
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void scale( const FieldT& a) {
    vfield_->scale( *a.vfield_ );
    sfield_->scale( *a.sfield_ );
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  constexpr Scalar dotLoc( const FieldT& a ) const {

    Scalar b = 0.;

    b = vfield_->dotLoc( *a.vfield_ ) + sfield_->dotLoc( *a.sfield_ );

    return( b );
  }


	/// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
	constexpr Scalar dot( const FieldT& y ) const {

		return( this->reduce( comm(), dotLoc( y ) ) );
	}


  /// \}
  /// \name Norm method
  /// \{

  /// \brief Compute the norm of the field.
  constexpr Scalar normLoc(  Belos::NormType type = Belos::TwoNorm ) const {

		return(
			(Belos::InfNorm==type)?
			std::max(vfield_->normLoc(type), sfield_->normLoc(type) ):
			(vfield_->normLoc(type) + sfield_->normLoc(type)) );
  }

 /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
	/// \todo include scaled norm
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
		return(
				vfield_->normLoc( *weights.vfield_ ) +
				sfield_->normLoc( *weights.sfield_ ) );
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
    vfield_->assign(*a.vfield_);
    sfield_->assign(*a.sfield_);
  }


  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed = false, int seed = 1) {
    vfield_->random();
    sfield_->random();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    vfield_->init(alpha);
    sfield_->init(alpha);
  }

	void initField() {
    vfield_->initField();
    sfield_->initField();
	}

	void level() const {
    vfield_->level();
    sfield_->level();
  }




  /// \}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& out=std::cout ) const {
    vfield_->print( out );
    sfield_->print( out );
  }


  void write( int count=0 ) const {
    vfield_->write(count);
    sfield_->write(count);
  }




}; // end of class CompoundField


/// \brief creates a compound vector+scalar field(vector).
///
/// \param vfield
/// \param sfield
/// \return Field vector
/// \relates CompoundField
template<class VField, class SField>
Teuchos::RCP< CompoundField<VField,SField> > createCompoundField(
    const Teuchos::RCP<VField>&  vfield, const Teuchos::RCP<SField>& sfield ) {

  return( Teuchos::RCP<CompoundField<VField,SField> > (
      new CompoundField<VField,SField>( vfield, sfield ) ) );
}


} // end of namespace Pimpact


#ifdef COMPILE_ETI
#include "Pimpact_VectorField.hpp"
#include "Pimpact_TimeField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"
extern template class Pimpact::CompoundField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> >, Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::CompoundField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> >, Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > >;
extern template class Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > >;
extern template class Pimpact::CompoundField< Pimpact::ModeField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >, Pimpact::ModeField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > >;
extern template class Pimpact::CompoundField< Pimpact::ModeField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >, Pimpact::ModeField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > >;
extern template class Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > >;
extern template class Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > >;
#endif


#endif // end of #ifndef PIMPACT_COMPOUNDFIELD_HPP
