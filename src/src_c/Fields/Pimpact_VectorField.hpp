#pragma once
#ifndef PIMPACT_VECTORFIELD_HPP
#define PIMPACT_VECTORFIELD_HPP


#include <iostream>
#include <vector>

#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_AbstractField.hpp"
#include "Pimpact_extern_ScalarField.hpp"
#include "Pimpact_extern_VectorField.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact {






/// \brief important basic Vector class  it wraps three ScalarFields.
/// \ingroup Field
/// \relates ScalarField
template<class SpaceType>
class VectorField : private AbstractField<SpaceType> {

public:

	using SpaceT = SpaceType;

protected:

	using Scalar = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;

	using ScalarArray = Scalar*;

	using SF = ScalarField<SpaceT>;

	ScalarArray s_;

	const bool owning_; /// < not template parameter

	SF sFields_[3];

	void allocate() {
		Ordinal n = getStorageSize()/3;
		s_ = new Scalar[3*n];
		for( int i=0; i<3; ++i )
			sFields_[i].setStoragePtr( s_+i*n );
	}

public:

	VectorField( const Teuchos::RCP< const SpaceT >& space, bool owning=true ):
		AbstractField<SpaceT>( space ),
		owning_(owning),
		sFields_{ {space,false,U}, {space,false,V}, {space,false,W} }
	{

			if( owning_ ) {
				allocate();
				initField();
			}
	};


	/// \brief copy constructor.
	///
	/// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
	/// \param vF
	/// \param copyType by default a ECopy::Shallow is done but allows also to deepcopy the field
	VectorField( const VectorField& vF, ECopy copyType=ECopy::Deep ):
		AbstractField<SpaceT>( vF.space() ),
		owning_(vF.owning_),
		sFields_{ {vF(U),copyType}, {vF(V),copyType}, {vF(W),copyType} }
	{

			//for( int i=0; i<3; ++i )
				//sFields_[i] = Teuchos::rcp( new SF( vF(i), copyType ) ); // copytype doesnot matter here, because it's not owning

			if( owning_ ) {

				allocate();

				switch( copyType ) {
					case ECopy::Shallow:
						initField();
						break;
					case ECopy::Deep:
						*this = vF;
						break;
				}
			}
		};


	~VectorField() { if( owning_ ) delete[] s_; }

	Teuchos::RCP<VectorField> clone( ECopy copyType=ECopy::Deep ) const {

		Teuchos::RCP<VectorField> vf = Teuchos::rcp( new VectorField( space() ) );

		switch( copyType ) {
			case ECopy::Shallow:
				break;
			case ECopy::Deep:
				*vf = *this;
				break;
		}

		return( vf );
	}

  /// \name Attribute methods
  /// @{


	/// \brief returns the length of Field.
	///
	/// the vector length is with regard to the inner points such that
	/// \f[ N_u = (N_x-1)(N_y-2)(N_z-2) \f]
	/// \f[ N_v = (N_x-2)(N_y-1)(N_z-2) \f]
	/// \f[ N_w = (N_x-2)(N_y-2)(N_z-1) \f]
	/// \return vect length \f[= N_u+N_v+N_w\f]
	constexpr Ordinal getLength() const {
		Ordinal n = 0;
		for( int i=0; i<SpaceT::sdim; ++i )
			n += sFields_[i].getLength();

		return( n );
	}


  /// @}
  /// \name Update methods
  /// @{

  /// \brief Replace \c this with \f$\alpha a + \beta b\f$.
  ///
  /// only inner points
	void add( const Scalar& alpha, const VectorField& a, const Scalar& beta, const
			VectorField& b, const B& wB=B::Y ) {

		// add test for consistent VectorSpaces in debug mode
		for( int i=0; i<SpaceT::sdim; ++i ) {
			sFields_[i].add( alpha, a(i), beta, b(i), wB );
		}
		changed();
	}


	/// \brief Put element-wise absolute values of source vector \c y into this
	/// vector.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
	/// \return Reference to this object
	void abs( const VectorField& y, const B& bcYes=B::Y ) {
		for( int i=0; i<SpaceT::sdim; ++i )
			sFields_[i].abs( y(i), bcYes );
		changed();
	}


	/// \brief Put element-wise reciprocal of source vector \c y into this vector.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
	/// \return Reference to this object
	void reciprocal( const VectorField& y, const B& bcYes=B::Y ) {
		// add test for consistent VectorSpaces in debug mode
		for( int i=0; i<SpaceT::sdim; ++i)
			sFields_[i].reciprocal( y(i), bcYes );
		changed();
	}


	/// \brief Scale each element of the vectors in \c this with \c alpha.
	void scale( const Scalar& alpha, const B& bcYes=B::Y ) {
		for(int i=0; i<SpaceT::sdim; ++i)
			sFields_[i].scale( alpha, bcYes );
		changed();
	}


	/// \brief Scale this vector <em>element-by-element</em> by the vector a.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
	/// \return Reference to this object
	void scale( const VectorField& a, const B& bcYes=B::Y ) {
		// add test for consistent VectorSpaces in debug mode
		for(int i=0; i<SpaceT::sdim; ++i)
			sFields_[i].scale( a(i), bcYes );
		changed();
	}


	/// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
	constexpr Scalar dotLoc ( const VectorField& a, const B& bcYes=B::Y ) const {
		Scalar b = 0.;

		for( int i=0; i<SpaceT::sdim; ++i )
			b += sFields_[i].dotLoc( a(i), bcYes );

		return( b );
	}


	/// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
	constexpr Scalar dot( const VectorField& y, const B& bcYes=B::Y ) const {

		return( this->reduce( comm(), dotLoc( y, bcYes ) ) );
	}


	/// @}
	/// \name Norm method
	/// @{

	constexpr Scalar normLoc( Belos::NormType type = Belos::TwoNorm, const B& bcYes=B::Y ) const {

		Scalar normvec = 0.;

		for( int i=0; i<SpaceT::sdim; ++i )
			normvec =
				(type==Belos::InfNorm)?
					std::max( sFields_[i].normLoc(type,bcYes), normvec ):
					( normvec+sFields_[i].normLoc(type,bcYes) );

		return( normvec );
	}


 /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  constexpr Scalar norm( Belos::NormType type = Belos::TwoNorm, const B& bcYes=B::Y ) const {

		Scalar normvec = this->reduce(
				comm(),
				normLoc( type, bcYes ),
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
	constexpr Scalar normLoc( const VectorField& weights, const B& bcYes=B::Y ) const {
		Scalar normvec = 0.;

		for( int i=0; i<SpaceT::sdim; ++i )
			normvec += sFields_[i].normLoc( weights(i), bcYes );

		return( normvec );
	}


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr Scalar norm( const VectorField& weights, const B& bcYes=B::Y ) const {
		return( std::sqrt( this->reduce( comm(), normLoc( weights, bcYes ) ) ) );
	}


	/// @}
	/// \name Initialization methods
  /// @{


	/// \brief *this := a
	///
	/// Assign (deep copy) a into mv.
	VectorField& operator=( const VectorField& a ) {

		for( int i=0; i<SpaceT::sdim; ++i)
			sFields_[i] = a(i);

		return *this;
	}


	/// \brief Replace the vectors with a random vectors.
	///
	/// depending on Fortrans \c Random_number implementation, with always same
	/// seed => not save, if good randomness is required
	void random( bool useSeed=false, const B& bcYes=B::Y, int seed=1 ) {

		for( int i=0; i<SpaceT::sdim; ++i )
			sFields_[i].random( useSeed, bcYes, seed );

		changed();
	}

	/// \brief Replace each element of the vector  with \c alpha.
	void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero(), const B& bcYes=B::Y ) {
		for( int i=0; i<SpaceT::sdim; ++i )
			sFields_[i].init( alpha, bcYes );
		changed();
	}


	/// \brief Replace each element of the vector with \c alpha[i].
	void init( const Teuchos::Tuple<Scalar,3>& alpha, const B& bcYes=B::Y ) {
		for( int i=0; i<SpaceT::sdim; ++i )
			sFields_[i].init( alpha[i], bcYes );
		changed();
	}


	///  \brief initializes VectorField including boundaries to zero 
	void initField() {
		for( int i=0; i<SpaceT::sdim; ++i )
			sFields_[i].initField();
	}

private:

	/// \brief kind of VectorField profile
	/// \relates VectorField::initField
	enum EVectorField {
		ZeroFlow,
		PoiseuilleFlow2D_inX,
		PoiseuilleFlow2D_inY,
		PoiseuilleFlow2D_inZ,
		Pulsatile2D_inXC,
		Pulsatile2D_inXS,
		Pulsatile2D_inYC,
		Pulsatile2D_inYS,
		Streaming2D,
		Circle2D,
		Circle2D_inXZ,
		RankineVortex2D,
		GaussianForcing1D,
		BoundaryFilter1D,
		GaussianForcing2D,
		BoundaryFilter2D,
		Streaming2DC,
		Streaming2DS,
		VPoint2D,
		Disc2D,
		RotationDisc2D,
		ConstFlow,
		SweptHiemenzFlow,
		Disturbance,
		ScalarFields,
		Couette,
		Cavity
	};

	/// \brief helper function getting \c EVectorField for switch statement from name 
	///
	/// \param name input name
	/// \return according int number
	EVectorField string2enum( const std::string& name ) {

		std::string lcName = name;
		std::transform(lcName.begin(), lcName.end(), lcName.begin(), ::tolower);

		if( "zero" == lcName ) return( ZeroFlow );
		else if( "constant" == lcName) return( ConstFlow );
		else if( "poiseuille" == lcName ) return( PoiseuilleFlow2D_inX );
		else if( "poiseuille in x" == lcName ) return( PoiseuilleFlow2D_inX );
		else if( "poiseuille in y" == lcName ) return( PoiseuilleFlow2D_inY );
		else if( "poiseuille in z" == lcName ) return( PoiseuilleFlow2D_inZ );
		else if( "pulsatile in x cos" == lcName ) return( Pulsatile2D_inXC );
		else if( "pulsatile in x sin" == lcName ) return( Pulsatile2D_inXS );
		else if( "pulsatile in y cos" == lcName ) return( Pulsatile2D_inYC );
		else if( "pulsatile in y sin" == lcName ) return( Pulsatile2D_inYS );
		else if( "streaming" == lcName ) return( Streaming2D );
		else if( "circle" == lcName ) return( Circle2D );
		else if( "circle xy" == lcName ) return( Circle2D );
		else if( "circle xz" == lcName ) return( Circle2D_inXZ );
		else if( "rankine vortex" == lcName ) return( RankineVortex2D );
		else if( "gaussian forcing 1d" == lcName ) return( GaussianForcing1D );
		else if( "boundary filter 1d" == lcName ) return( BoundaryFilter1D );
		else if( "gaussian forcing 2d" == lcName ) return( GaussianForcing2D );
		else if( "boundary filter 2d" == lcName ) return( BoundaryFilter2D );
		else if( "streaming 2d cos" == lcName ) return( Streaming2DC );
		else if( "streaming 2d sin" == lcName ) return( Streaming2DS );
		else if( "v point 2d" == lcName ) return( VPoint2D );
		else if( "disc 2d" == lcName ) return( Disc2D );
		else if( "rotation disc 2d" == lcName ) return( RotationDisc2D );
		else if( "shbl" == lcName ) return( SweptHiemenzFlow );
		else if( "swept hiemenz flow" == lcName ) return( SweptHiemenzFlow );
		else if( "disturbance" == lcName ) return( Disturbance );
		else if( "scalar" == lcName ) return( ScalarFields );
		else if( "couette" == lcName ) return( Couette );
		else if( "cavity" == lcName ) return( Cavity );
		else {
			const bool& Flow_Type_not_known = true; 
			TEUCHOS_TEST_FOR_EXCEPT( Flow_Type_not_known );
		}
		return( ZeroFlow ); // just to please the compiler
	}

public:

	///  \brief initializes including boundaries to zero 
	void initField( Teuchos::ParameterList& para ) {

		EVectorField type =
			string2enum( para.get<std::string>( "Type", "zero" ) );

		switch( type ) {
			case ZeroFlow : {
				for( int i=0; i<SpaceT::sdim; ++i )
					sFields_[i].initField( ConstField );
				break;
			}
			case ConstFlow : {
				sFields_[U].initField( ConstField, para.get<Scalar>( "U", 1.) );
				sFields_[V].initField( ConstField, para.get<Scalar>( "V", 1.) );
				sFields_[W].initField( ConstField, para.get<Scalar>( "W", 1.) );
				break;
			}
			case PoiseuilleFlow2D_inX : {
				for( int i=0; i<SpaceT::sdim; ++i )
					if( U==i )
						sFields_[i].initField( Poiseuille2D_inY );
					else
						sFields_[i].initField( ConstField );
				break;
			}
			case PoiseuilleFlow2D_inY : {
				for( int i=0; i<SpaceT::sdim; ++i )
					if( V==i )
						sFields_[i].initField( Poiseuille2D_inX );
					else
						sFields_[i].initField( ConstField );
				break;
			}
			case PoiseuilleFlow2D_inZ : {
				for( int i=0; i<SpaceT::sdim; ++i )
					if( W==i )
						sFields_[i].initField( Poiseuille2D_inX );
					else
						sFields_[i].initField( ConstField );
				break;
			}
			case Pulsatile2D_inXC :
				VF_init_2DPulsatileXC(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(1),
						space()->getCoordinatesLocal()->getX(Y,EField::S),
						space()->getDomainSize()->getRe(),     // TODO: verify
						space()->getDomainSize()->getAlpha2(), // TODO: verify
						para.get<Scalar>( "px", 1. ),          // TODO: verify
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case Pulsatile2D_inYC :
				VF_init_2DPulsatileYC(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(0),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getDomainSize()->getRe(),     // TODO: verify
						space()->getDomainSize()->getAlpha2(), // TODO: verify
						para.get<Scalar>( "px", 1. ),          // TODO: verify
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case Pulsatile2D_inXS :
				VF_init_2DPulsatileXS(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(1),
						space()->getCoordinatesLocal()->getX(Y,EField::S),
						space()->getDomainSize()->getRe(),     // TODO: verify
						space()->getDomainSize()->getAlpha2(), // TODO: verify
						para.get<Scalar>( "px", 1. ),          // TODO: verify
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case Pulsatile2D_inYS :
				VF_init_2DPulsatileYS(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(0),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getDomainSize()->getRe(),     // TODO: verify
						space()->getDomainSize()->getAlpha2(), // TODO: verify
						para.get<Scalar>( "px", 1. ),          // TODO: verify
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case Streaming2D :
				VF_init_StreamingS(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(0),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getDomainSize()->getRe(),     // TODO: verify
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case Streaming2DC :
				VF_init_StreamingC(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(0),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getDomainSize()->getRe(),     // TODO: verify
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case Streaming2DS:
				VF_init_StreamingS(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(0),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getDomainSize()->getRe(),     // TODO: verify
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case Circle2D : {
				sFields_[U].initField( Grad2D_inY, -1. );
				sFields_[V].initField( Grad2D_inX,  1. );
				sFields_[W].initField();
				break;
			}
			case Circle2D_inXZ : {
				sFields_[U].initField( Grad2D_inY, -1. );
				sFields_[V].initField();
				sFields_[W].initField( Grad2D_inX,  1. );
				break;
			}
			case RankineVortex2D :
				VF_init_RankineVortex(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getCoordinatesLocal()->getX(Y,EField::S),
						space()->getCoordinatesLocal()->getX(X,EField::U),
						space()->getCoordinatesLocal()->getX(Y,EField::V),
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case GaussianForcing1D :
				VF_init_GaussianForcing1D(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(0),
						space()->getCoordinatesLocal()->getX(X,EField::U),
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case BoundaryFilter1D :
				VF_init_BoundaryFilter1D(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(0),
						space()->getCoordinatesLocal()->getX(X,EField::U),
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case GaussianForcing2D :
				VF_init_GaussianForcing2D(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getCoordinatesLocal()->getX(Y,EField::S),
						space()->getCoordinatesLocal()->getX(X,EField::U),
						space()->getCoordinatesLocal()->getX(Y,EField::V),
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case BoundaryFilter2D :
				VF_init_BoundaryFilter2D(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getCoordinatesLocal()->getX(Y,EField::S),
						space()->getCoordinatesLocal()->getX(X,EField::U),
						space()->getCoordinatesLocal()->getX(Y,EField::V),
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case VPoint2D :
				VF_init_Vpoint(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getDomainSize()->getSize(),
						space()->getCoordinatesLocal()->getX(X,EField::U),
						space()->getCoordinatesLocal()->getX(Y,EField::S),
						space()->getDomainSize()->getRe(),     // TODO: verify
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case Disc2D :
				VF_init_Disc(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getCoordinatesLocal()->getX(Y,EField::S),
						space()->getCoordinatesLocal()->getX(Z,EField::S),
						space()->getCoordinatesLocal()->getX(X,EField::U),
						space()->getCoordinatesLocal()->getX(Y,EField::V),
						//re, om, px,sca,
						para.get<Scalar>( "center x", 1. ),
						para.get<Scalar>( "center y", 1. ),
						para.get<Scalar>( "radius", 1. ),
						para.get<Scalar>( "sca", 0.1 ),
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case RotationDisc2D :
				VF_init_RotatingDisc(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getCoordinatesLocal()->getX(X,EField::S),
						space()->getCoordinatesLocal()->getX(Y,EField::S),
						para.get<Scalar>( "center x", 1. ),
						para.get<Scalar>( "center y", 1. ),
						para.get<Scalar>( "omega", 1. ),
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			case SweptHiemenzFlow : {
				Scalar pi = 4.*std::atan(1.);

				std::cout << "hello\n";

				Ordinal nTemp = space()->gu(X) - space()->gl(X) + 1;
				//for( Ordinal i=0; i<=space()->nLoc(X); ++i ) {
				//std::cout << "i: " << i<< " (\t";
				//for( Ordinal ii=0; ii<nTemp; ++ii )
				//std::cout << space()->getInterpolateV2S()->getC(X)[ i*nTemp + ii ] << ",\t";
				//std::cout << ")\n";
				//}
				//Scalar c[6] = { 
				//space()->getInterpolateV2S()->getC(X)[ nTemp + 0 ],
				//space()->getInterpolateV2S()->getC(X)[ nTemp + 1 ],
				//space()->getInterpolateV2S()->getC(X)[ nTemp + 2 ],
				//space()->getInterpolateV2S()->getC(X)[ nTemp + 3 ],
				//space()->getInterpolateV2S()->getC(X)[ nTemp + 4 ],
				//space()->getInterpolateV2S()->getC(X)[ nTemp + 5 ] };
				//std::cout << "c\n";
				//for( Ordinal ii=0; ii<nTemp; ++ii )
				//std::cout << c[  ii ] << ",\t";
				//std::cout << "c\n";

				VF_init_SHBF( 
						space()->rankST(),      
						space()->getShift(0),
						space()->getProcGrid()->getIB(0),
						space()->nGlo(),
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->dl(),
						space()->du(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getCoordinatesGlobal()->getX( X, EField::S ),
						space()->getCoordinatesGlobal()->getX( X, EField::U ),
						space()->getCoordinatesLocal()->getX( Z, EField::W ),
						space()->getInterpolateV2S()->getC(X)+2*nTemp-1, /// \todo rm dirty hack
						space()->getDomainSize()->getRe(),
						para.get<int>( "nonDim", 0 ),
						para.get<Scalar>( "kappa", 0. ),
						para.get<Scalar>( "seep angle", 0. ),
						para.get<Scalar>( "seep angle", 0. )*pi/180.,
						para.get<Scalar>( "attack angle", 0. ),
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );

				std::cout << "hello\n"
					<< space()->getInterpolateV2S()->getC(X)[9] << "\n";

				break;
			}
			case Disturbance : {
				VF_init_Dist(
						space()->rankST(),      
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(U),
						space()->eIndB(U),
						space()->sIndB(V),
						space()->eIndB(V),
						space()->sIndB(W),
						space()->eIndB(W),
						space()->getBCGlobal()->getBCL( Z ),
						space()->getCoordinatesLocal()->getX( X, EField::U ),
						space()->getCoordinatesLocal()->getX( X, EField::S ),
						space()->getCoordinatesLocal()->getX( Y, EField::S ),
						space()->getCoordinatesLocal()->getX( Z, EField::W ),
						space()->getCoordinatesLocal()->getX( Z, EField::S ),
						3, // dist_type,          
						0.15, // vortex_ampli_prim,  
						3., // vortex_x1pos,       
						3., // vortex_x3pos,       
						3., // vortex_radius,      
						10, // vortex_band,        
						sFields_[U].getRawPtr(),
						sFields_[V].getRawPtr(),
						sFields_[W].getRawPtr() );
				break;
			}
			case ScalarFields : {
				sFields_[U].initField( para.sublist( "U" ) );
				sFields_[V].initField( para.sublist( "V" ) );
				sFields_[W].initField( para.sublist( "W" ) );
				break;
			}
			case Couette : {
				sFields_[U].initFromFunction( 
						[]( Scalar x, Scalar y, Scalar z ) -> Scalar {
						return( y );
						} );
				break;
			}
			case Cavity : {
				sFields_[U].initFromFunction( 
						[]( Scalar x, Scalar y, Scalar z ) -> Scalar {
						if( std::fabs(x-1.)<0.5 )
						return( 1. );
						else
						return( 0. );
						} );
				break;
			}
		}
		changed();
	}


	/// \brief extrapolates on the boundaries such that it is zero
	/// \note dirty hack(necessary for TripleCompostion)
	void extrapolateBC( const Belos::ETrans& trans=Belos::NOTRANS ) {
		for( int i=0; i<SpaceT::sdim; ++i )
			sFields_[i].extrapolateBC( trans );
	}


	/// dirty hack(necessary for MG)
  void level() const {}

  /// @}


  /// \brief Print the vector.  To be used for debugging only.
  void print( std::ostream& out=std::cout ) const {
    for(int i=0; i<SpaceT::sdim; ++i) {
      sFields_[i].print( out );
    }
  }


	/// \brief writes vector field to hdf5 file
	///
	/// \param count
	/// \param restart decides if velocity is interpolated to pressure points
	/// \todo implement/test restart and read
	void write( int count=0, bool restart=false ) const {

		if( 0==space()->rankS() ) {
			Teuchos::Tuple<Ordinal,3> N;
			for( int i=0; i<3; ++i ) {
				N[i] = space()->nGlo(i);
				if( space()->getBCGlobal()->getBCL(i)==Pimpact::PeriodicBC )
					N[i] = N[i]-1;
			}
			std::ofstream xfile;
			std::ostringstream ss;
			ss << std::setw( 5 ) << std::setfill( '0' ) << count;
			//        std::string fname = "v_"+ss.str();
			xfile.open( "vel_"+ ss.str() +".xmf", std::ofstream::out );
			xfile<< "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.1\">\n";
			xfile << "\t<Domain>\n";
			xfile << "\t\t<Grid Name=\"3DRectMesh\" GridType=\"Uniform\">\n";
			xfile << "\t\t\t<Topology TopologyType=\"3DRectMesh\" Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\"/>\n";
			xfile << "\t\t\t<Geometry GeometryType=\"VXVYVZ\">\n";
			xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
			xfile << "\t\t\t\t\tDimensions=\""<< N[0] << "\"\n";
			xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
			xfile << "\t\t\t\t\tPrecision=\"8\"\n";
			xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
			xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/VectorX\n";
			xfile << "\t\t\t\t</DataItem>\n";
			xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
			xfile << "\t\t\t\t\tDimensions=\""<< N[1] << "\"\n";
			xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
			xfile << "\t\t\t\t\tPrecision=\"8\"\n";
			xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
			xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/VectorY\n";
			xfile << "\t\t\t\t</DataItem>\n";
			xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
				xfile << "\t\t\t\t\tDimensions=\""<< N[2] << "\"\n";
			xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
			xfile << "\t\t\t\t\tPrecision=\"8\"\n";
			xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
			xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/VectorZ\n";
			xfile << "\t\t\t\t</DataItem>\n";
			xfile << "\t\t\t</Geometry>\n";
			xfile << "\t\t\t<Attribute Name=\"VelX\" AttributeType=\"Scalar\" Center=\"Node\">\n";
			xfile << "\t\t\t\t<DataItem Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
			xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/velX\n";
			xfile << "\t\t\t\t</DataItem>\n";
			xfile << "\t\t\t</Attribute>\n";
			xfile << "\t\t\t<Attribute Name=\"VelY\" AttributeType=\"Scalar\" Center=\"Node\">\n";
			xfile << "\t\t\t\t<DataItem Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
			xfile << "\t\t\t\t\tvelY_"<< ss.str() << ".h5:/velY\n";
			xfile << "\t\t\t\t</DataItem>\n";
			xfile << "\t\t\t</Attribute>\n";
			xfile << "\t\t\t<Attribute Name=\"VelZ\" AttributeType=\"Scalar\" Center=\"Node\">\n";
			xfile << "\t\t\t\t<DataItem Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
			xfile << "\t\t\t\t\tvelZ_"<< ss.str() << ".h5:/velZ\n";
			xfile << "\t\t\t\t</DataItem>\n";
			xfile << "\t\t\t</Attribute>\n";
			xfile << "\t\t</Grid>\n";
			xfile << "\t</Domain>\n";
			xfile << "</Xdmf>\n";
			xfile.close();
		}
		for( int i=0; i<SpaceT::sdim; ++i )
			sFields_[i].write( count, restart );
	}


public:

	/// \name storage methods.
	/// \brief highly dependent on underlying storage should only be used by
	/// Operator or on top field implementer.
	/// @{ 

	constexpr Ordinal getStorageSize() const { return( sFields_[0].getStorageSize()*3 ); }

  void setStoragePtr( Scalar*  array ) {
    s_ = array;
    Ordinal n = sFields_[0].getStorageSize();
    for( int i=0; i<3; ++i )
      sFields_[i].setStoragePtr( s_+i*n );
  }

	constexpr Scalar* getRawPtr() { return( s_ ); }

	constexpr const Scalar* getConstRawPtr() const { return( s_ ); }

	/// \deprecated
  Scalar* getRawPtr ( int i )       { return( sFields_[i].getRawPtr() ); }

	/// \deprecated
	constexpr const Scalar* getConstRawPtr ( int i )  const  { return( sFields_[i].getConstRawPtr() ); }

	///  @} 

//protected:

  //SF& getField( int i ) { return( sFields_[i]); }
	//constexpr const SF&  getConstField( int i ) const { return( sFields_[i] ); }

public:

  SF& operator()( int i ) {
		assert( 0<=i );
		assert( i< 3 );
		return( sFields_[i] );
	}

  constexpr const SF& operator()( int i ) const {
		assert( 0<=i );
		assert( i< 3 );
		return( sFields_[i] );
	}


  constexpr const Teuchos::RCP<const SpaceT>& space() const { return( AbstractField<SpaceT>::space_ ); }

  constexpr const MPI_Comm& comm() const { return(space()->comm()); }


  /// \name comunication methods.
  /// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
  ///
  /// @{
	
  void changed( const int& vel_dir, const int& dir ) const {
    sFields_[vel_dir].changed( dir );
  }

  void changed() const {
    for( int vel_dir=0; vel_dir<SpaceT::sdim; ++vel_dir )
      for( int dir=0; dir<SpaceT::sdim; ++dir ) {
        changed( vel_dir, dir );
      }
  }

  void exchange( const int& vel_dir, const int& dir ) const {
    sFields_[vel_dir].exchange(dir);
  }

  void exchange() const {
    for( int vel_dir=0; vel_dir<SpaceT::sdim; ++vel_dir )
      for( int dir=0; dir<SpaceT::sdim; ++dir )
        exchange( vel_dir, dir );
  }

  void setExchanged( const int& vel_dir, const int& dir ) const {
    sFields_[vel_dir].setExchanged( dir );
  }

	void setExchanged() const {
		for( int vel_dir=0; vel_dir<SpaceT::sdim; ++vel_dir )
			for( int dir=0; dir<SpaceT::sdim; ++dir ) {
				setExchanged( vel_dir, dir );
			}
	}

  /// @}


}; // end of class VectorField




} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::VectorField< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::VectorField< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::VectorField< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::VectorField< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_VECTORFIELD_HPP
