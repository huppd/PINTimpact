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
template<class IFT>
class MultiHarmonicField : private AbstractField<typename IFT::SpaceT> {

public:

  using SpaceT = typename IFT::SpaceT;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using FieldT = MultiHarmonicField<IFT>;

  using AF = AbstractField<SpaceT>;

	const bool global_;

  Teuchos::RCP<IFT> field0_;

  Teuchos::Array< Teuchos::RCP< ModeField<IFT> > > fields_;

	Scalar* s_;

  mutable bool exchangedState_;


	void allocate() {

		if( global_ ) {

			Ordinal nx = field0_->getStorageSize();

			s_ = new Scalar[ getStorageSize() ];
			
			field0_->setStoragePtr( s_ );

			for( Ordinal i=0; i<space()->nGlo(3); ++i )
				fields_[i]->setStoragePtr( s_ + nx + 2*nx*i );
		}
		else{

			Ordinal nx = field0_->getStorageSize();

			s_ = new Scalar[ getStorageSize() ];

			field0_->setStoragePtr( s_ );
			for( Ordinal i=0; i<=space()->end(U,3) - std::max(space()->begin(U,3),1); ++i )
				fields_[i]->setStoragePtr( s_ + nx + 2*nx*i );
		}
	}


public:

	constexpr Ordinal getStorageSize() const {
		return(
				( global_ )?
					( ( 1 + 2*space()->nGlo(3)                                             )*field0_->getStorageSize() ):
					( ( 1 + 2*( space()->end(U,3) - std::max(space()->begin(U,3),1) + 1 )  )*field0_->getStorageSize() )
				);
	}


	MultiHarmonicField( const Teuchos::RCP<const SpaceT>& space ):
		MultiHarmonicField( space, space->np(3)==1 ) {};

	MultiHarmonicField( const Teuchos::RCP<const SpaceT>& space, const bool global ):
		AF( space ),
		global_(global),
		field0_( Teuchos::rcp( new IFT(space,false) ) ),
		fields_( space->nGlo(3) ),
    exchangedState_( true ) {

			TEUCHOS_TEST_FOR_EXCEPT( 4 != SpaceT::dimension  );
			TEUCHOS_TEST_FOR_EXCEPT( true != space()->getStencilWidths()->spectralT() );

			if( global_ ) {
				for( Ordinal i=0; i<space->nGlo(3); ++i )
					fields_[i] = Teuchos::rcp( new ModeField<IFT>( space, false ) );
			}
			else{
				for( Ordinal i=0; i<space()->end(U,3) - std::max(space()->begin(U,3),1) + 1; ++i )
					fields_[i] = Teuchos::rcp( new ModeField<IFT>( space, false ) );
			}

			allocate();
			initField();
	};


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ECopy::Shallow is done but allows also to deepcopy the field
	MultiHarmonicField( const MultiHarmonicField& vF, ECopy copyType=ECopy::Deep ):
		AF( vF.space() ),
		//global_( vF.global_ ),
		global_( vF.space()->np(3)==1 ),
		field0_( Teuchos::rcp( new IFT( *vF.field0_, copyType ) ) ),
		fields_( vF.space()->nGlo(3) ),
    exchangedState_(vF.exchangedState_) {

			if( global_ ) {
				for( Ordinal i=1; i<=space()->nGlo(3); ++i )
					fields_[i-1] = Teuchos::rcp( new ModeField<IFT>( vF.getConstField(i), copyType ) );
			}
			else{
				for( Ordinal i=0; i<space()->end(U,3) - std::max(space()->begin(U,3),1) + 1; ++i )
					fields_[i] = Teuchos::rcp( new ModeField<IFT>( vF.getConstField(i+std::max(space()->begin(U,3),1)), copyType ) );
			}

			allocate();
			switch( copyType ) {
				case ECopy::Shallow:
					initField();
					break;
				case ECopy::Deep:
					for( int i=0; i<getStorageSize(); ++i )
						s_[i] = vF.s_[i];
					break;
			}
	};


	~MultiHarmonicField() { delete[] s_; }

  Teuchos::RCP<FieldT> clone( ECopy ctype=ECopy::Deep ) const {

    return( Teuchos::rcp( new FieldT(*this,ctype) ) );
  }


  /// \name Attribute methods
  /// \{

protected:

	constexpr Ordinal index( const Ordinal& i ) const {
		return( i - 1 + 
				(( global_||0==space()->begin(U,3) )?
					0:
					(-space()->begin(U,3)+1) )
				);
	};

public:

	constexpr const bool& global() const { return( global_ ); }

  constexpr IFT&                    get0Field()               { return( *field0_ ); }
  constexpr const IFT&              getConst0Field()    const { return( *field0_ ); }
  constexpr Teuchos::RCP<      IFT> get0FieldPtr()            { return(  field0_ ); }
  constexpr Teuchos::RCP<const IFT> getConst0FieldPtr() const { return(  field0_ ); }

  constexpr ModeField<IFT>&                     getField        ( const Ordinal& i )       { return( *fields_[index(i)] ); }
  constexpr const ModeField<IFT>&               getConstField   ( const Ordinal& i ) const { return( *fields_[index(i)] ); }
  constexpr Teuchos::RCP<      ModeField<IFT> > getFieldPtr     ( const Ordinal& i )       { return( fields_[index(i)]  ); }
  constexpr Teuchos::RCP<const ModeField<IFT> > getConstFieldPtr( const Ordinal& i ) const { return( fields_[index(i)]  ); }

  constexpr IFT&                    getCField        ( const Ordinal& i )       { return( fields_[index(i)]->getCField()      ); }
  constexpr const IFT&              getConstCField   ( const Ordinal& i ) const { return( fields_[index(i)]->getConstCField() ); }
  constexpr Teuchos::RCP<      IFT> getCFieldPtr     ( const Ordinal& i )       { return( fields_[index(i)]->getCFieldPtr()   ); }
  constexpr Teuchos::RCP<const IFT> getConstCFieldPtr( const Ordinal& i ) const { return( fields_[index(i)]->getConstCFieldPtr() ); }

  constexpr IFT&                    getSField        ( const Ordinal& i )       { return( fields_[index(i)]->getSField()         ); }
  constexpr const IFT&              getConstSField   ( const Ordinal& i ) const { return( fields_[index(i)]->getConstSField()    ); }
  constexpr Teuchos::RCP<      IFT> getSFieldPtr     ( const Ordinal& i )       { return( fields_[index(i)]->getSFieldPtr()      ); }
  constexpr Teuchos::RCP<const IFT> getConstSFieldPtr( const Ordinal& i ) const { return( fields_[index(i)]->getConstSFieldPtr() ); }


  constexpr const Teuchos::RCP<const SpaceT>& space() const { return( AF::space_ ); }


  constexpr const MPI_Comm& comm() const { return( space()->getProcGrid()->getCommWorld() ); }


  /// \brief returns the length of Field.
  ///
  /// the vector length is with regard to the inner points
  constexpr Ordinal getLength() const {

		Ordinal len = 0;

		len += getConst0FieldPtr()->getLength();
		//len += space()->nGlo(3)*getConstFieldPtr( std::max(space()->begin(U,3),1) )->getLength();
		len += 2*space()->nGlo(3)*getConst0FieldPtr()->getLength();

    return( len );
  }


  /// \brief get number of stored Field's
  constexpr int getNumberVecs() const { return( 1 ); }


  /// \}
  /// \name Update methods
  /// \{

  /// \brief Replace \c this with \f$\alpha A + \beta B\f$.
	/// \todo add test for consistent VectorSpaces in debug mode
  void add( const Scalar& alpha, const FieldT& A, const Scalar& beta, const FieldT& B ) {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->add(alpha, A.getConst0Field(), beta, B.getConst0Field() );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getFieldPtr(i)->add( alpha, A.getConstField(i), beta, B.getConstField(i) );

		changed();
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void abs( const FieldT& y) {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->abs( y.getConst0Field() );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getFieldPtr(i)->abs( y.getConstField(i) );

		changed();
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal( const FieldT& y){

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->reciprocal( y.getConst0Field() );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getFieldPtr(i)->reciprocal( y.getConstField(i) );

		changed();
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale( const Scalar& alpha ) {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->scale( alpha );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getFieldPtr(i)->scale( alpha );

		changed();

  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void scale( const FieldT& a) {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->scale( a.getConst0Field() );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getFieldPtr(i)->scale( a.getConstField(i) );

		changed();
  }




  /// \}
  /// \name Norm method and SP
  /// \{

  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  constexpr Scalar dotLoc( const FieldT& a ) const {

    Scalar b = 0.;

		if( 0==space()->begin(U,3) )
			b += getConst0FieldPtr()->dotLoc( a.getConst0Field() );
		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			b += getConstFieldPtr(i)->dotLoc( a.getConstField(i) );

    return( b );
  }

	/// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
	constexpr Scalar dot( const FieldT& y ) const {

		return( this->reduce( comm(), dotLoc( y ) ) );
	}

  /// \brief Compute the norm of Field.
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  constexpr Scalar normLoc( Belos::NormType type = Belos::TwoNorm ) const {

    Scalar normvec = 0.;

		if( 0==space()->begin(U,3) )
			normvec += getConst0FieldPtr()->normLoc(type);

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			normvec =
				(Belos::InfNorm==type)?
				std::max( getConstFieldPtr(i)->normLoc(type), normvec ):
				(normvec+getConstFieldPtr(i)->normLoc(type));

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

    Scalar normvec= Teuchos::ScalarTraits<Scalar>::zero();

		if( 0==space()->begin(U,3) )
			normvec += getConst0FieldPtr()->normLoc( weights.getConst0Field() );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			normvec += getConstFieldPtr(i)->normLoc(weights.getConstField(i));

    return( normvec );
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

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->assign( a.getConst0Field() );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getFieldPtr(i)->assign( a.getConstField(i) );

		changed();
  }


  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed = false, int seed = 1) {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->random();

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getFieldPtr(i)->random();

		changed();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->init( alpha );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getFieldPtr(i)->init( alpha );

		changed();
  }


	void initField() {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->initField();

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getFieldPtr(i)->initField();

		changed();
	}


	///  \brief initializes including boundaries to zero 
	void initField( Teuchos::ParameterList& para ) {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->initField( para.sublist("0 mode") );

		if( space()->begin(U,3)<=1 && 1<=space()->end(U,3) ) {
			getCFieldPtr(1)->initField( para.sublist("cos mode") );
			getSFieldPtr(1)->initField( para.sublist("sin mode") );
		}
		changed();
	}

	void extrapolateBC( const Belos::ETrans& trans=Belos::NOTRANS ) {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->extrapolateBC( trans );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getConstFieldPtr(i)->extrapolateBC( trans );

		changed();
  }

  void level()  {

		if( 0==space()->begin(U,3) )
			get0FieldPtr()->level();

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getConstFieldPtr(i)->level();

		changed();
  }

  /// \}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os=std::cout ) const {

		if( 0==space()->begin(U,3) )
			getConst0FieldPtr()->print( os );

		for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
			getConstFieldPtr(i)->print( os );
  }


  void write( int count=0, bool time_evol_yes=false ) const {

		if( time_evol_yes ) {
			exchange();

			if( space()->getProcGrid()->getIB(3)==1 ) {
				Scalar pi = 4.*std::atan(1.);
				Ordinal nf = space()->nGlo(3);
				Ordinal nt = 4*nf;
				Teuchos::RCP<IFT> temp = getConst0FieldPtr()->clone( Pimpact::ECopy::Shallow );
				for( Ordinal i=0; i<nt;  ++i ) {
					temp->assign( getConst0Field() );
					//				temp->initField(  );
					for( Ordinal j=1; j<=nf; ++j ) {
						temp->add(
								1., *temp,
								std::sin( 2.*pi*i*((Scalar)j)/nt ), getConstSField(j) );
						temp->add(
								std::cos( 2.*pi*i*((Scalar)j)/nt ), getConstCField(j),
								1., *temp );
						temp->write( count+i );
					}
				}
			}
		}
		else{

			if( 0==space()->begin(U,3) )
				getConst0FieldPtr()->write(count);

			for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i ) {
				getConstCFieldPtr(i)->write( count+2*i-1 );
				getConstSFieldPtr(i)->write( count+2*i   );
			}
		}
  }


	/// \name comunication methods.
	/// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
	///
	/// \{

	void changed() const { exchangedState_ = false; }


	void exchange() const {

#ifndef NDEBUG
		TEUCHOS_TEST_FOR_EXCEPT( !global_ );
#endif

		// check if exchange is necessary
		if( exchangedState_==false ) {

			// exchange spatial
			if( 0==space()->begin(U,3) )
				get0FieldPtr()->exchange();

			for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
				getConstFieldPtr(i)->exchange();

			// mpi stuff
			int nx = getConst0FieldPtr()->getStorageSize();

			// --- sendcount ---
			int sendcount = 0;
			if( 0==space()->begin(U,3) )
				sendcount += nx;
			for( Ordinal i=std::max(space()->begin(U,3),1); i<=space()->end(U,3); ++i )
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
			if( 0!=space()->begin(U,3) )
				getConst0FieldPtr()->setExchanged();

			for( Ordinal i=1; i<space()->begin(U,3); ++i )
				getConstFieldPtr(i)->setExchanged();

			for( Ordinal i=space()->end(U,3)+1; i<=space()->nGlo(3); ++i )
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



/// \brief creates a multi-harmonic scalar field.
///
/// \relates MultiHarmonicField
/// \param space scalar Vector Space to which returned vector belongs
/// \param global 
template<class FieldT>
Teuchos::RCP< MultiHarmonicField< FieldT > > createMultiHarmonic(
    const Teuchos::RCP<const typename FieldT::SpaceT >& space, bool global ) {

	return( Teuchos::rcp( new  MultiHarmonicField<FieldT>( space, global ) ) );
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
