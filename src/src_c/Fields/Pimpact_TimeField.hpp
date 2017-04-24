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
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \brief templated class which is the interface to \c Belos and \c NOX
///
/// has multiple \c Field's, where \c Field can be a \c Pimpact:ScalarField, \c
/// Pimpact:VectorField or some combination with \c Pimpact::ModeField or \c
/// Pimpact::CompoundField \note if this is heavily used for many Field's, then
/// the implementation should be improved such that communication is done such
/// that only done once per FieldT not per Field
///
/// \todo rm RCP from vector
/// \ingroup Field
template<class Field>
class TimeField : private AbstractField<typename Field::SpaceT> {

	template<class Op,bool CNY>
	friend class TimeOpWrap;
	template<class SpaceTT>
	friend class DtTimeOp;

public:

	using SpaceT = typename Field::SpaceT;

	using Scalar = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;

	using FieldT = Pimpact::TimeField<Field>;

	using ScalarArray = Scalar*;

	using AF = AbstractField<SpaceT>;

protected:

	Teuchos::Array< Teuchos::RCP<Field> > mfs_; 

	ScalarArray array_;

	mutable bool exchangedState_;

public:

	TimeField( Teuchos::RCP<const SpaceT> space, F dummy=F::S ):
		AF( space ), exchangedState_(true) {

			Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

			mfs_ = Teuchos::Array< Teuchos::RCP<Field> >( nt );

			for( int i=0; i<nt; ++i )
				mfs_[i] = Teuchos::rcp( new Field( space, false ) );

			Ordinal nx = at(0).getStorageSize();

			array_ = new Scalar[nx*nt];

			for( int i=0; i<nx*nt; ++i )
				array_[i] = 0.;

			for( int i=0; i<nt; ++i )
				at(i).setStoragePtr( array_+i*nx );
	}



  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param field 
  /// \param copyType by default a ECopy::Shallow is done but allows also to deepcopy the field
	TimeField( const TimeField& field, ECopy copyType=ECopy::Deep ):
		AF( field.space() ),
		exchangedState_( field.exchangedState_ ) {

			Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

			mfs_ = Teuchos::Array< Teuchos::RCP<Field> >(nt);

			for( int i=0; i<nt; ++i )
				 mfs_[i] = Teuchos::rcp( new Field( space(), false ) );

			Ordinal nx = mfs_[0]->getStorageSize();

			array_ = new Scalar[nx*nt];

			for( int i=0; i<nt; ++i )
				at(i).setStoragePtr( array_+i*nx );

			switch( copyType ) {
				case( ECopy::Deep ) : {
					for( int i=0; i<nt; ++i )
						at(i) = field(i);
					break;
				}
				case( ECopy::Shallow ) : {
					for( int i=0; i<nt*nx; ++i )
						array_[i] = 0.;
					exchangedState_ = true;
					break;
				}
			}
	}



	~TimeField() { delete[] array_; }


	/// \brief Create a new \c TimeField with
	Teuchos::RCP< FieldT > clone( ECopy ctype = ECopy::Deep ) const {
		Teuchos::RCP< FieldT > mv_ = Teuchos::rcp( new FieldT(*this,ctype) );
		return( mv_ );
	}


	/// \brief returns the length of Field.
	constexpr Ordinal getLength() {
		return( space()->nGlo(3)*at(0).getLength() );
	}

public:

	/// \}
	/// \name Update methods
	/// \{


	/// \brief <tt>mv := alpha*a + beta*b</tt>
	void add( Scalar alpha, const FieldT& a, Scalar beta, const FieldT& b, const B& wb=B::Y ) {

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).add( alpha, a(i), beta, b(i), wb );
		changed();
	}


	/// \brief Put element-wise absolute values of source vector \c y into this
	/// vector.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
	/// \return Reference to this object
	void abs( const FieldT& y) {

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).abs( y(i) );
		changed();
	}


	/// \brief Put element-wise reciprocal of source vector \c y into this vector.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
	/// \return Reference to this object
	void reciprocal( const FieldT& y){

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).reciprocal( y(i) );
		changed();
	}


	/// \brief Scale each element of every \c Field by \c gamma.
	///
	/// Here x represents on \c Field, and we update it as
	/// \f[ x_i = \alpha x_i \quad \mbox{for } i=1,\dots,n \f]
	void scale( const Scalar& alpha, const B& wB=B::Y ) {

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).scale( alpha, wB );
		changed();
	}


	/// \brief Scale this vector <em>element-by-element</em> by the vector a.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i = x_i \cdot y_i \quad \mbox{for } i=1,\dots,n \f]
	/// \return Reference to this object
	void scale( const FieldT& y) {
		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).scale( y(i) );
		changed();
	}

	/// \}


	/// \brief Compute the inner product for the \c TimeField considering it as one Vector.
	constexpr Scalar dotLoc( const FieldT& A ) const {

		Scalar b = 0.;

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			b+= at(i).dotLoc( A(i) );

		return( b );
	}

	/// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
	constexpr Scalar dot( const FieldT& y ) const {

		return( this->reduce( comm(), dotLoc( y ) ) );
	}


	/// \brief Compute the norm for the \c TimeField as it is considered as one Vector .
	constexpr Scalar normLoc( Belos::NormType type = Belos::TwoNorm, const B& bcYes=B::Y  ) const {

		Scalar normvec = 0.;

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			normvec = 
				( (Belos::InfNorm==type)?
				std::max( at(i).normLoc(type, bcYes), normvec ) :
				(normvec + at(i).normLoc(type, bcYes) ) );

		return( normvec );
	}

 /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm/
  constexpr Scalar norm( Belos::NormType type = Belos::TwoNorm, const B& bcYes=B::Y  ) const {

		Scalar normvec = this->reduce(
				comm(),
				normLoc( type, bcYes ),
				(Belos::InfNorm==type)?MPI_MAX:MPI_SUM );

		normvec = (
			(Belos::TwoNorm==type) ?
				std::sqrt(normvec) :
				normvec );

    return( normvec );
  }


	/// \brief Weighted 2-Norm.
	///
	/// Here x represents this vector, and we compute its weighted norm as follows:
	/// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
	/// \return \f$ \|x\|_w \f$
	constexpr Scalar normLoc( const FieldT& weights ) const {

		Scalar nor = Teuchos::ScalarTraits<Scalar>::zero();

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			nor+= at(i).norm( weights(i) );

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


	/// \brief *this := a.
	///
	/// assign (deep copy) A into mv.
	TimeField& operator=( const TimeField& a ) {

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i) = a(i);
		changed();

		return *this;
	}


	/// \brief Replace the vectors with a random vectors.
	void random(bool useSeed = false, int seed = 1) {

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).random();
		changed();
	}


	/// \brief \f[ *this = \alpha \f]
	void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero(), const B& wB=B::Y ) {

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).init(alpha,wB);
		changed();
	}

	void extrapolateBC( const Belos::ETrans& trans=Belos::NOTRANS ) {

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).extrapolateBC( trans );
		changed();
	}

	void level() const {

		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).level();
		changed();
	}


	/// \param os
	void print( std::ostream& os=std::cout ) const {
		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).print( os );
	}



	void write( int count=0 ) const {
		for( Ordinal i=space()->si(F::S,3); i<=space()->ei(F::S,3); ++i )
			at(i).write(count++ + space()->getShift(3) );
	}


	constexpr const MPI_Comm& comm() const { return( space()->getProcGrid()->getCommWorld() ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return( AF::space_ ); }

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

				Ordinal lengthL = transL * at(0).getStorageSize();
				Ordinal lengthU = transU * at(0).getStorageSize();

				Scalar* ghostUR = at( space()->si(F::S,3)-transU).getRawPtr();
				Scalar* ghostLR = at( space()->ei(F::S,3)  +transL).getRawPtr();

				Scalar* ghostUS = at(space()->ei(F::S,3)  ).getRawPtr();
				Scalar* ghostLS = at(space()->si(F::S,3)).getRawPtr();

				if( transL>0 ) MPI_Irecv( ghostUR, lengthL, MPI_REAL8, rankL, 1, comm(), &reqL);
				if( transU>0 ) MPI_Irecv( ghostLR, lengthU, MPI_REAL8, rankU, 2, comm(), &reqU);

				if( transL>0 ) MPI_Send ( ghostUS, lengthL, MPI_REAL8, rankU, 1, comm() );
				if( transU>0 ) MPI_Send ( ghostLS, lengthU, MPI_REAL8, rankL, 2, comm() );

				if( transL>0 ) MPI_Wait( &reqL, &statusL );
				if( transU>0 ) MPI_Wait( &reqU, &statusU );

				// depends on if field from sender was exchanged, so to be sure
				at( 0                    ).changed();
				at( space()->ei(F::S,3) ).changed();

			}
			else {
				if( std::abs( space()->bl(3) )>0 ) {
					*mfs_[ space()->si(F::S,3)-1 ] = at( space()->ei(F::S,3) );
					at( space()->si(F::S,3)-1 ).changed();
				}
				if( std::abs( space()->bu(3) )>0 ) {
					*mfs_[ space()->ei(F::S,3)+1 ] = at( space()->si(F::S,3) );
					at( space()->ei(F::S,3)+1 ).changed();
				}
			}
		}
		exchangedState_ = true;
	}


protected:
	                Field& at( const int& i ) { return( *mfs_[i] ); }
	constexpr const Field& at( const int& i ) { return( *mfs_[i] ); }

public:

	                Field& operator()( const int& i ) { return( at(i) ); }
	constexpr const Field& operator()( const int& i ) { return( at(i) ); }

	constexpr ScalarArray getRawPtr() { return( array_ ); }

	constexpr const Scalar* getConstRawPtr() const { return( array_ ); }

}; // end of class TimeField



/// \brief factory for \c TimeField
/// \relates TimeField
/// \deprecated
template<class FieldT, class SpaceT>
Teuchos::RCP< TimeField<FieldT> >
createTimeField( const Teuchos::RCP<const SpaceT>& space ) {

	return( Teuchos::rcp( new TimeField<FieldT>( space ) ) );
}



/// \deprecated 
template<class SpaceT>
void
initVectorTimeField(
		TimeField<VectorField<SpaceT> >& field,
		EFlowType flowType=Zero2DFlow,
		typename SpaceT::Scalar xm=0.5,
		typename SpaceT::Scalar ym=0.5,
		typename SpaceT::Scalar rad=0.1,
		typename SpaceT::Scalar amp=0.25 ) {

	using S = typename SpaceT::Scalar;
	//using O = typename SpaceT::Ordinal;

	Teuchos::RCP<const SpaceT> space = field.space();

	S pi = 4.*std::atan(1.);

	S nt = space->nGlo(3);
	S offset = space->getShift(3) - space->si(F::S,3);

	bool notImplemented = true;
	TEUCHOS_TEST_FOR_EXCEPT( notImplemented );
	//for( O i=space->si(F::S,3); i<=space->ei(F::S,3); ++i )
		//switch( flowType ) {
			//case Zero2DFlow:
				//field(i)->initField( ZeroFlow );
				//break;
			//case Const2DFlow:
				//field(i)->initField( ConstFlow, xm, ym, rad );
				//break;
			//case Poiseuille_inX:
				//field(i)->initField( PoiseuilleFlow2D_inX );
				//break;
			//case Poiseuille_inY:
				//field(i)->initField( PoiseuilleFlow2D_inY );
				//break;
			//case Streaming2DFlow: {
				//S ampt = std::sin( 2.*pi*((F::S)i+offset)/nt );
				//field(i)->initField( Streaming2D, ampt );
				//break;
			//}
			//case OscilatingDisc2D: {
				////			std::cout << "\ti: " << i << "\tt: " << 2.*pi*((F::S)i+offset)/nt << "\tt: " << space->getCoordinatesLocal()->getX(  F::S, ECoord::T )[i] << "\n";
				//S ymt = ym+amp*std::sin( space->getCoordinatesLocal()->getX( F::S, ECoord::T )[i] );
				//S xmt = xm;
				//field(i)->initField( Disc2D, xmt, ymt, rad );
				//break;
			//}
			//case OscilatingDisc2DVel: {
				//S yvelt = amp*std::cos( 2.*pi*((F::S)i+offset)/nt );
				//S xvelt = 0;
				//field(i)->init( Teuchos::tuple( xvelt, yvelt, 0.) );
				//break;
			//}
			//case ConstVel_inX:{
				//field(i)->init( Teuchos::tuple( -2*xm*std::cos(space->getCoordinatesLocal()->getX( F::S, ECoord::T )[i]), 0., 0.) ); // here xm = p
				//break;
			//}
			//case Pulsatile_inX: {
				////field(i)->initField( Pulsatile2D_inX, xm, space->getCoordinatesLocal()->getX( F::S, ECoord::T )[i], ym, rad); // the arguments are (xmt,i,ymt,rad) --> (re,t,px,alpha)
				//break;
			//}
			//default:
				//field(i)->initField( ZeroFlow );
				//break;
		//}

	//field.changed();
}


} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > >;
extern template class Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_TIMEFIELD_HPP
