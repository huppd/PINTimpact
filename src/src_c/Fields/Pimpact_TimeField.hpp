#pragma once
#ifndef PIMPACT_TIMEFIELD_HPP
#define PIMPACT_TIMEFIELD_HPP


#include <cmath>
#include <vector>

#include "mpi.h"

#include <Teuchos_Array.hpp>
#include "Teuchos_RCP.hpp"
#include <Teuchos_Range1D.hpp>

#include <BelosTypes.hpp>

#include "Pimpact_AbstractField.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"

#include "Pimpact_Types.hpp"



namespace Pimpact {



/// \brief templated class which is the interface to \c Belos and \c NOX
///
/// has multiple \c Field's, where \c Field can be a \c Pimpact:ScalarField, \c
/// Pimpact:VectorField or some combination with \c Pimpact::ModeField or \c
/// Pimpact::CompoundField \note if this is heavily used for many Field's, then
/// the implementation should be improved such that communication is done such
/// that only done once per MV not per Field
///
/// \todo decide if it is better to modifie ghostlayer or exchange eventuell more
/// \todo maybe move functionality in Scalar/VectorField<S,O,4>
/// \ingroup Field
template<class Field>
class TimeField : private AbstractField<typename Field::SpaceT> {

	template<class Op,bool CNY>
	friend class TimeOpWrap;
	template<class SpaceTT>
	friend class DtTimeOp;
	template<class SpaceTT, bool CNY>
	friend class TimeNonlinearJacobian;

public:

	typedef typename Field::Scalar Scalar;
	typedef typename Field::Ordinal Ordinal;

	static const int dimension = Field::dimension;

	typedef typename Field::SpaceT SpaceT;

	typedef Pimpact::TimeField<Field> MV;

	typedef Scalar* ScalarArray;

	typedef AbstractField<SpaceT> AF;


	Teuchos::Array< Teuchos::RCP<Field> > mfs_;

public:

	typedef Teuchos::Array< Teuchos::RCP<Field> > FieldArray;
	typedef typename FieldArray::iterator Iter;

protected:

	ScalarArray array_;

	mutable bool exchangedState_;

public:

	TimeField( Teuchos::RCP<const SpaceT> space, EField dummy=EField::S ):
		AF( space ), exchangedState_(true) {

			Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

			mfs_ = Teuchos::Array< Teuchos::RCP<Field> >( nt );

			for( int i=0; i<nt; ++i )
				mfs_[i] = Teuchos::rcp( new Field( space, false ) );

			Ordinal nx = mfs_[0]->getStorageSize();

			array_ = new Scalar[nx*nt];

			for( int i=0; i<nx*nt; ++i )
				array_[i] = 0.;

			for( int i=0; i<nt; ++i )
				mfs_[i]->setStoragePtr( array_+i*nx );

		}



  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param field 
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	TimeField( const TimeField& field, ECopyType copyType=DeepCopy ):
		AF( field.space() ),
		exchangedState_( field.exchangedState_ ) {

			Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

			mfs_ = Teuchos::Array< Teuchos::RCP<Field> >(nt);

			for( int i=0; i<nt; ++i )
				// mfs_[i] = Teuchos::rcp( new Field( space(), false ) );
				mfs_[i] = field.mfs_[i]->clone(copyType);

			Ordinal nx = mfs_[0]->getStorageSize();

			array_ = new Scalar[nx*nt];

			for( int i=0; i<nt; ++i )
				mfs_[i]->setStoragePtr( array_+i*nx );

			if( DeepCopy==copyType )
				for( int i=0; i<nt; ++i ) {
					mfs_[i]->assign( *(field.mfs_[i]) );
				}
			else {
				for( int i=0; i<nt*nx; ++i )
					array_[i] = 0.;
				exchangedState_ = true;
			}

		}



	~TimeField() { delete[] array_; }


	/// \brief Create a new \c TimeField with
	Teuchos::RCP< MV > clone( ECopyType ctype = DeepCopy ) const {
		auto mv_ = Teuchos::rcp( new MV(*this,ctype) );
		return( mv_ );
	}


	/// \brief returns the length of Field.
	///
	/// \param noxVec if \c TimeField is used for NOX the Vector length is
	/// considered for all Fields
	Ordinal getLength( bool noxVec=true ) const {
		return( space()->nGlo()[3]*mfs_[0]->getLength(noxVec) );
	}


public:

	/// \brief is true
	bool HasConstantStride() const { return( true ); }
	int getNumberVecs() const {  return( mfs_.size() ); }
	/// \}
	/// \name Update methods
	/// \{


	/// \brief <tt>mv := alpha*A + beta*B</tt>
	///
	///	The Tpetra specialization of this method ignores and completely
	///	overwrites any NaN or Inf entries in A.  Thus, it does <i>not</i> mean
	///	the same thing as <tt>mv := 0*mv + alpha*A + beta*B</tt> in IEEE 754
	///	floating-point arithmetic. (Remember that NaN*0 = NaN.)
	void add( Scalar alpha, const MV& A, Scalar beta, const MV& B ) {

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->add( alpha, *A.mfs_[i], beta, *B.mfs_[i] );
		changed();

	}


	/// \brief Put element-wise absolute values of source vector \c y into this
	/// vector.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
	/// \return Reference to this object
	void abs(const MV& y) {

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->abs( *y.mfs_[i] );
		changed();

	}


	/// \brief Put element-wise reciprocal of source vector \c y into this vector.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
	/// \return Reference to this object
	void reciprocal(const MV& y){
		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->reciprocal( *y.mfs_[i] );
		changed();
	}


	/// \brief Scale each element of every \c Field by \c gamma.
	///
	/// Here x represents on \c Field, and we update it as
	/// \f[ x_i = \alpha x_i \quad \mbox{for } i=1,\dots,n \f]
	void scale( const Scalar& alpha ) {
		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->scale(alpha);
		changed();
	}


	/// \brief Scale this vector <em>element-by-element</em> by the vector a.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i = x_i \cdot y_i \quad \mbox{for } i=1,\dots,n \f]
	/// \return Reference to this object
	void scale(const MV& y) {
		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->scale( *y.mfs_[i] );
		changed();
	}

	/// \}



	/// \brief Compute the inner product for the \c TimeField considering it as one Vector.
	Scalar dot( const MV& A, bool global=true ) const {

		Scalar b = 0.;

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			b+= mfs_[i]->dot( *A.mfs_[i], false );

		if( global ) this->reduceNorm( comm(), b );

		return( b );

	}



	/// \brief Compute the norm for the \c TimeField as it is considered as one Vector .
	Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

		Scalar normvec = 0.;

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i ) {
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

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			nor+= mfs_[i]->norm( *weights.mfs_[i] );

		if( global ) this->reduceNorm( comm(), nor, Belos::TwoNorm );

		return( nor );
	}


	/// \brief mv := A.
	///
	/// assign (deep copy) A into mv.
	void assign( const MV& A ) {

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->assign( *A.mfs_[i] );
		changed();

	}


	/// \brief Replace the vectors with a random vectors.
	void random(bool useSeed = false, int seed = 1) {

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->random();
		changed();

	}


	/// \brief \f[ *this = \alpha \f]
	void init( Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->init(alpha);
		changed();

	}

	void initField() {
		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->initField();
		changed();
	}

	void setCornersZero() const {
		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->setCornersZero();
		changed();
	}

	void level() const {

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->level();
		changed();

	}


	/// \param os
	void print( std::ostream& os=std::cout ) const {
		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->print( os );
	}



	void write( int count=0 ) const {
		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i )
			mfs_[i]->write(count++ + space()->getShift(3) );
	}


	const MPI_Comm& comm() const { return( space()->getProcGrid()->getCommWorld() ); }

	Teuchos::RCP<const SpaceT> space() const { return( AF::space_ ); }

	public:

	void changed() const {
		exchangedState_ = false;
	}

	/// \note shoud be constant but weirdly then Iter becomes const iter and can't be converted to int
	void exchange() const {

		if( !exchangedState_ ) {
			if( space()->np(3)>1 ) {

				int transL = std::abs( space()->bl(3) );
				int transU = std::abs( space()->bu(3) );

				// std::cout << "transL: " <<  transL << "\n";
				// std::cout << "transU: " <<  transU << "\n";

				int rankU = space()->getProcGrid()->getRankU(3);
				int rankL = space()->getProcGrid()->getRankL(3);

				MPI_Request reqL;
				MPI_Request reqU;

				MPI_Status statusL;
				MPI_Status statusU;

				Ordinal lengthL = transL * mfs_[0]->getStorageSize();
				Ordinal lengthU = transU * mfs_[0]->getStorageSize();

				Scalar* ghostUR = mfs_[0]->getRawPtr();
				Scalar* ghostLR = mfs_[space()->eInd(S,3)]->getRawPtr();

				Scalar* ghostUS = mfs_[space()->eInd(S,3)-transL]->getRawPtr();
				Scalar* ghostLS = mfs_[space()->sInd(S,3)       ]->getRawPtr();

				if( transL>0 ) MPI_Irecv( ghostUR, lengthL, MPI_REAL8, rankL, 1, comm(), &reqL);
				if( transU>0 ) MPI_Irecv( ghostLR, lengthU, MPI_REAL8, rankU, 2, comm(), &reqU);

				if( transL>0 ) MPI_Send ( ghostUS, lengthL, MPI_REAL8, rankU, 1, comm() );
				if( transU>0 ) MPI_Send ( ghostLS, lengthU, MPI_REAL8, rankL, 2, comm() );

				if( transL>0 ) MPI_Wait( &reqL, &statusL );
				if( transU>0 ) MPI_Wait( &reqU, &statusU );

				// depends on if field from sender was exchanged, so to be sure
				mfs_[ 0                  ]->changed();
				mfs_[ space()->eInd(S,3) ]->changed();

			}
			else {
				if( std::abs( space()->bl(3) )>0 ) {
					mfs_[space()->sInd(S,3)-1]->assign( *mfs_[space()->eInd(S,3)-1] );
					mfs_[space()->sInd(S,3)-1]->changed();
				}
				if( std::abs( space()->bu(3) )>0 ) {
					mfs_[space()->eInd(S,3)]->assign( *mfs_[space()->sInd(S,3)] );
					mfs_[space()->eInd(S,3)]->changed();
				}
			}
		}

		exchangedState_ = true;
	}


	Teuchos::RCP<Field> getFieldPtr( int i ) { return(  mfs_[i] ); }
	Field& getField   ( int i ) { return( *mfs_[i] ); }


	Teuchos::RCP<const Field> getConstFieldPtr( int i ) const { return(  mfs_[i] ); }
	const Field&  getConstField   ( int i ) const { return( *mfs_[i] ); }


	ScalarArray getRawPtr() { return( array_ ); }


	const Scalar* getConstRawPtr() const { return( array_ ); }


}; // end of class TimeField



/// \brief factory for \c TimeField
/// \relates TimeField
/// \todo substitue
/// \deprecated
template<class FieldT, class SpaceT>
Teuchos::RCP< TimeField<FieldT> >
createTimeField( const Teuchos::RCP<const SpaceT>& space ) {

	return( Teuchos::rcp( new TimeField<FieldT>( space ) ) );

}



/// \todo move to initField
template<class SpaceT>
void
initVectorTimeField(
		Teuchos::RCP<TimeField<VectorField<SpaceT> > > field,
		EFlowType flowType=Zero2DFlow,
		typename SpaceT::Scalar xm=0.5,
		typename SpaceT::Scalar ym=0.5,
		typename SpaceT::Scalar rad=0.1,
		typename SpaceT::Scalar amp=0.25 ) {

	typedef typename SpaceT::Scalar S;
	typedef typename SpaceT::Ordinal O;

	auto space = field->space();

	S pi = 4.*std::atan(1.);

	S nt = space->nGlo(3);
	S offset = space->getShift(3) - space->sInd(EField::S,3);

	for( O i=space->sInd(EField::S,3); i<space->eInd(EField::S,3); ++i )
		switch( flowType ) {
			case Zero2DFlow:
				field->getFieldPtr(i)->initField( ZeroFlow );
				break;
			case Const2DFlow:
				field->getFieldPtr(i)->initField( ConstFlow, xm, ym, rad );
				break;
			case Poiseuille_inX:
				field->getFieldPtr(i)->initField( PoiseuilleFlow2D_inX );
				break;
			case Poiseuille_inY:
				field->getFieldPtr(i)->initField( PoiseuilleFlow2D_inY );
				break;
			case Streaming2DFlow: {
															S ampt = std::sin( 2.*pi*((S)i+offset)/nt );
															field->getFieldPtr(i)->initField( Streaming2D, ampt );
															break;
														}
			case OscilatingDisc2D: {
															 //			std::cout << "\ti: " << i << "\tt: " << 2.*pi*((S)i+offset)/nt << "\tt: " << space->getCoordinatesLocal()->getX( ECoord::T, EField::S )[i] << "\n";
															 S ymt = ym+amp*std::sin( space->getCoordinatesLocal()->getX( ECoord::T, EField::S )[i] );
															 S xmt = xm;
															 field->getFieldPtr(i)->initField( Disc2D, xmt, ymt, rad );
															 break;
														 }
			case OscilatingDisc2DVel: {
																	S yvelt = amp*std::cos( 2.*pi*((S)i+offset)/nt );
																	S xvelt = 0;
																	field->getFieldPtr(i)->init( Teuchos::tuple( xvelt, yvelt, 0.) );
																	break;
																}
			case ConstVel_inX:{
													field->getFieldPtr(i)->init( Teuchos::tuple( -2*xm*std::cos(space->getCoordinatesLocal()->getX( ECoord::T, EField::S )[i]), 0., 0.) ); // here xm = p
													break;
												}
			case Pulsatile_inX: {
														//field->getFieldPtr(i)->initField( Pulsatile2D_inX, xm, space->getCoordinatesLocal()->getX( ECoord::T, EField::S )[i], ym, rad); // the arguments are (xmt,i,ymt,rad) --> (re,t,px,alpha)
														break;
													}
			default:
													field->getFieldPtr(i)->initField( ZeroFlow );
													break;
		}

	field->changed();

}


} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > >;
extern template class Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_TIMEFIELD_HPP
