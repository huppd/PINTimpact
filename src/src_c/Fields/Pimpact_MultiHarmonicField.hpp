#pragma once
#ifndef PIMPACT_MULTIHARMONICFIELD_HPP
#define PIMPACT_MULTIHARMONICFIELD_HPP


#include <vector>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "BelosTypes.hpp"


#include "Pimpact_AbstractField.hpp"
#include "Pimpact_ModeField.hpp"




namespace Pimpact {


/// \brief important basic Vector class.
///
/// vector for wrapping many fields into one multiharmonic field
/// \ingroup Field
template<class Field>
class MultiHarmonicField : private AbstractField<typename Field::SpaceT> {

public:

  using SpaceT = typename Field::SpaceT;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  static const int dimension = SpaceT::dimension;

protected:

  using MV = MultiHarmonicField<Field>;

  using AF = AbstractField<SpaceT>;

  Teuchos::RCP<Field> field0_;

  Teuchos::Array< Teuchos::RCP< ModeField<Field> > > fields_;

	Scalar* s_;

  mutable bool exchangedState_;

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

	/// \todo add etst
  MultiHarmonicField(
			const Teuchos::RCP<const SpaceT>& space ):
		AF( space ),
		field0_( Teuchos::rcp( new Field(space,false) ) ),
		fields_( space->nGlo(3) ),
    exchangedState_( true ) {

			TEUCHOS_TEST_FOR_EXCEPT( 4 != dimension  );
			TEUCHOS_TEST_FOR_EXCEPT( true != space()->getStencilWidths()->spectralT() );

			for( Ordinal i=0; i<space->nGlo(3); ++i )
				fields_[i] = Teuchos::rcp( new ModeField<Field>( space, false ) );

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
		fields_( vF.space()->nGlo(3) ),
    exchangedState_(vF.exchangedState_) {

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

  const MPI_Comm& comm() const { return( space()->getProcGrid()->getCommWorld() ); }

  /// \brief returns the length of Field.
  ///
  /// the vector length is with regard to the inner points
  Ordinal getLength( bool nox_vec=false ) const {

		Ordinal len = 0;
		len += field0_->getLength(true);

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
	/// \todo add test for consistent VectorSpaces in debug mode
  void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {

		if( space()->sInd(U,3)<0 )
			field0_->add(alpha, *A.field0_, beta, *B.field0_);

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getFieldPtr(i)->add( alpha, A.getConstField(i), beta, B.getConstField(i) );

		changed();

  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void abs(const MV& y) {

		if( space()->sInd(U,3)<0 )
			field0_->abs( *y.field0_ );

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getFieldPtr(i)->abs( y.getConstField(i) );

		changed();

  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal(const MV& y){

		if( space()->sInd(U,3)<0 )
			field0_->reciprocal( *y.field0_ );

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getFieldPtr(i)->reciprocal( y.getConstField(i) );

		changed();

  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale( const Scalar& alpha ) {

		if( space()->sInd(U,3)<0 )
			field0_->scale(alpha);

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getFieldPtr(i)->scale(alpha);

		changed();

  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void scale(const MV& a) {

		if( space()->sInd(U,3)<0 )
			field0_->scale( *a.field0_ );

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getFieldPtr(i)->scale( a.getConstField(i) );

		changed();
  }




  /// \}
  /// \name Norm method and SP
  /// \{

  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  Scalar dot ( const MV& a, bool global=true ) const {

    Scalar b = 0.;

		if( space()->sInd(U,3)<0 )
			b += field0_->dot( *a.field0_, false );
		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			b += getConstFieldPtr(i)->dot( a.getConstField(i), false );

    if( global ) this->reduceNorm( comm(), b );

    return( b );

  }

  /// \brief Compute the norm of Field.
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = 0.;

		if( space()->sInd(U,3)<0 )
			normvec += field0_->norm(type,false);

    switch(type) {
    case Belos::OneNorm:
		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
				normvec += getConstFieldPtr(i)->norm(type,false);
      break;
    case Belos::TwoNorm:
		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
				normvec += getConstFieldPtr(i)->norm(type,false);
      break;
    case Belos::InfNorm:
		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
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

    double normvec= 0.;

		if( space()->sInd(U,3)<0 )
			normvec += field0_->norm(*weights.field0_,false);

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
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

		if( space()->sInd(U,3)<0 )
			field0_->assign(*a.field0_);

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getFieldPtr(i)->assign( a.getConstField(i) );

		changed();
  }


  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed = false, int seed = 1) {

		if( space()->sInd(U,3)<0 )
			field0_->random();

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getFieldPtr(i)->random();

		changed();

  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {

		if( space()->sInd(U,3)<0 )
			field0_->init( alpha );

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getFieldPtr(i)->init(alpha);

		changed();

  }


	void initField() {
		if( space()->sInd(U,3)<0 )
			field0_->initField();
		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getFieldPtr(i)->initField();
		changed();
	}


	///  \brief initializes including boundaries to zero 
	void initField( Teuchos::ParameterList& para ) {

		if( space()->sInd(U,3)<0 )
			field0_->initField( para.sublist("0 mode") );

		if( space()->sInd(U,3)<=0 && 0<space()->eInd(U,3) ) {
			getCFieldPtr(0)->initField( para.sublist("cos mode") );
			getSFieldPtr(0)->initField( para.sublist("sin mode") );
		}
		changed();
	}


	void setCornersZero() const {

		if( space()->sInd(U,3)<0 )
			field0_->setCornersZero();

		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getConstFieldPtr(i)->setCornersZero();

		changed();
	}


  void level() const {

		if( space()->sInd(U,3)<0 )
			field0_->level();
		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getConstFieldPtr(i)->level();

		changed();

  }

  /// \}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os=std::cout ) const {

		if( space()->sInd(U,3)<0 )
			field0_->print( os );
		for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
			getConstFieldPtr(i)->print( os );
  }

  void write( int count=0, bool time_evol_yes=false ) const {

		if( time_evol_yes ) {
			exchange();

			if( space()->getProcGrid()->getIB(3)==1 ) {
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

		}
		else{
			if( space()->sInd(U,3)<0 )
				field0_->write(count);
			for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i ) {
				getConstCFieldPtr(i)->write( count+2*i+1 );
				getConstSFieldPtr(i)->write( count+2*i+2 );
			}
		}
  }


	/// \name comunication methods.
	/// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
	///
	/// \{

  void changed() const {
    exchangedState_ = false;
  }

	/// \todo biggest todo
	void exchange() const {

		// check if exchange is necessary
		if( exchangedState_==false ) {

			// exchange spatial
			if( space()->sInd(U,3)<0 )
				field0_->exchange();

			for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
				getConstFieldPtr(i)->exchange();

			// mpi stuff
			int nx = field0_->getStorageSize();

			// --- sendcount ---
			int sendcount = 0;
			if( space()->sInd(U,3)<0 )
				sendcount += nx;
			for( Ordinal i=std::max(space()->sInd(U,3),0); i<space()->eInd(U,3); ++i )
				sendcount += 2*nx;

			// --- recvcount, displacement ---
			int np = space()->getProcGrid()->getNP(3);

			int* recvcounts = new int[np];
			int* displs     = new int[np];

			Ordinal nfl  = (space()->getGridSizeLocal()->get(3)+1)/np;
			Ordinal rem  = (space()->getGridSizeLocal()->get(3)+1)%np;

//			std::cout << "\trankST: " << space()->rankST() << "\tsendcount: " << sendcount << "\n";

//			MPI_Barrier( space()->getProcGrid()->getCommSlice(3) );

			for( int rank=0; rank<np; ++rank ) {
				if( 0==rank )
					recvcounts[rank] = nx + (nfl-1)*2*nx + (rank<rem?2*nx:0);
				else
					recvcounts[rank] = nfl*2*nx + (rank<rem?2*nx:0);
				if( 0==rank )
					displs[rank] = 0;
				else
					displs[rank] = displs[rank-1] + recvcounts[rank-1];
//				if( 0==space()->rankST() )
//					std::cout << "\trank: " << rank
//						<< "\trecevcount: " << recvcounts[rank]
//						<< "\tdispl: " << displs[rank] << "\n";
			}

			// exchange modes
			MPI_Allgatherv(
					MPI_IN_PLACE,                           // sendbuf
					sendcount,                              // sendcount
					MPI_REAL8,                              // sendtype
					s_,                                     // recvbuf
					recvcounts,                             // recvcounts
					displs,                                 // displs
					MPI_REAL8,                              // recvtype
					space()->getProcGrid()->getCommSlice(3) // comm
					); 

			delete[] recvcounts;
			delete[] displs;

			// set non owning spatial block as exchanged
			if( space()->sInd(U,3)>0 )
				field0_->setExchanged();

			for( Ordinal i=0; i<space()->sInd(U,3); ++i )
				getConstFieldPtr(i)->setExchanged();

			for( Ordinal i=space()->eInd(U,3); i<space()->nGlo(3); ++i )
				getConstFieldPtr(i)->setExchanged();
		}

	}

	/// \}

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
extern template class Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > >;
extern template class Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_MULTIHARMONICFIELD_HPP
