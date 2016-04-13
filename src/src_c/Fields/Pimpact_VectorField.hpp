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



/// \brief kind of VectorField profile
/// \relates VectorField::initField
enum EVectorField {
  ZeroFlow=0,
  PoiseuilleFlow2D_inX=1,
	PoiseuilleFlow2D_inY=2,
	PoiseuilleFlow2D_inZ=20,
  Pulsatile2D_inXC=3,
	Pulsatile2D_inXS=5,
  Pulsatile2D_inYC=4,
	Pulsatile2D_inYS=6,
  Streaming2D=7,
  Circle2D=8,
  Circle2D_inXZ=21,
  RankineVortex2D=9,
  GaussianForcing1D=10,
  BoundaryFilter1D=11,
  GaussianForcing2D=12,
  BoundaryFilter2D=13,
  Streaming2DC=14,
  Streaming2DS=15,
  VPoint2D=16,
  Disc2D=17,
  RotationDisc2D=18,
  ConstFlow=19,
	SweptHiemenzFlow=22,
	Disturbance=23,
	ScalarFields=24
};



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

	using FieldT = VectorField<SpaceT>;
	using SF = ScalarField<SpaceT>;

	ScalarArray s_;

	const bool owning_;

	Teuchos::Tuple< Teuchos::RCP<SF>, 3 > sFields_;

	void allocate() {
		Ordinal n = getStorageSize()/3;
		s_ = new Scalar[3*n];
		for( int i=0; i<3; ++i )
			sFields_[i]->setStoragePtr( s_+i*n );
	}

public:

	VectorField( const Teuchos::RCP< const SpaceT >& space, bool owning=true ):
		AbstractField<SpaceT>( space ),
		owning_(owning) {

			for( int i=0; i<3; ++i )
				sFields_[i] = Teuchos::rcp( new SF( space, false, static_cast<EField>(i) ) );

			if( owning_ ) {
				allocate();
				initField();
			}

	};


	/// \brief copy constructor.
	///
	/// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
	/// \param vF
	/// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
	VectorField( const VectorField& vF, ECopyType copyType=DeepCopy ):
		AbstractField<SpaceT>( vF.space() ),
		owning_(vF.owning_) {

			for( int i=0; i<3; ++i )
				sFields_[i] = Teuchos::rcp( new SF( vF.getConstField(i), copyType ) ); // copytype doesnot matter here, because it's not owning

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


	~VectorField() { if( owning_ ) delete[] s_; }

	Teuchos::RCP<FieldT> clone( ECopyType copyType=DeepCopy ) const {

		Teuchos::RCP<FieldT> vf = Teuchos::rcp( new FieldT( space() ) );

		switch( copyType ) {
			case ShallowCopy:
				break;
			case DeepCopy:
				for( int i=0; i<getStorageSize(); ++i )
					vf->getRawPtr()[i] = s_[i];
				break;
		}

		return( vf );
	}

  /// \name Attribute methods
  /// @{


	/// \brief returns the length of Field.
	///
	/// the vector length is withregard to the inner points such that
	/// \f[ N_u = (N_x-1)(N_y-2)(N_z-2) \f]
	/// \f[ N_v = (N_x-2)(N_y-1)(N_z-2) \f]
	/// \f[ N_w = (N_x-2)(N_y-2)(N_z-1) \f]
	/// \return vect length \f[= N_u+N_v+N_w\f]
	constexpr Ordinal getLength() const {
		Ordinal n = 0;
		for( int i=0; i<space()->dim(); ++i )
			n += sFields_[i]->getLength();

		return( n );
	}


	/// \brief get number of stored Field's
	constexpr int getNumberVecs() const { return( 1 ); }


  /// @}
  /// \name Update methods
  /// @{

  /// \brief Replace \c this with \f$\alpha A + \beta B\f$.
  ///
  /// only inner points
	void add( const Scalar& alpha, const FieldT& A, const Scalar& beta, const FieldT& B ) {
		// add test for consistent VectorSpaces in debug mode
		for( int i=0; i<space()->dim(); ++i ) {
			sFields_[i]->add( alpha, *A.sFields_[i], beta, *B.sFields_[i] );
		}
		changed();
	}


	/// \brief Put element-wise absolute values of source vector \c y into this
	/// vector.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
	/// \return Reference to this object
	void abs(const FieldT& y) {
		for( int i=0; i<space()->dim(); ++i )
			sFields_[i]->abs( *y.sFields_[i] );
		changed();
	}


	/// \brief Put element-wise reciprocal of source vector \c y into this vector.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
	/// \return Reference to this object
	void reciprocal(const FieldT& y){
		// add test for consistent VectorSpaces in debug mode
		for( int i=0; i<space()->dim(); ++i)
			sFields_[i]->reciprocal( *y.sFields_[i] );
		changed();
	}


	/// \brief Scale each element of the vectors in \c this with \c alpha.
	void scale( const Scalar& alpha ) {
		for(int i=0; i<space()->dim(); ++i)
			sFields_[i]->scale( alpha );
		changed();
	}


	/// \brief Scale this vector <em>element-by-element</em> by the vector a.
	///
	/// Here x represents this vector, and we update it as
	/// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
	/// \return Reference to this object
	void scale(const FieldT& a) {
		// add test for consistent VectorSpaces in debug mode
		for(int i=0; i<space()->dim(); ++i)
			sFields_[i]->scale( *a.sFields_[i] );
		changed();
	}


	/// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
	constexpr Scalar dotLoc ( const FieldT& a ) const {
		Scalar b = 0.;

		for( int i=0; i<space()->dim(); ++i )
			b += sFields_[i]->dotLoc( *a.sFields_[i] );

		return( b );

	}

	/// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
	constexpr Scalar dot( const FieldT& y ) const {

		return( this->reduce( comm(), dotLoc( y ) ) );

	}


	/// @}
	/// \name Norm method
	/// @{

	constexpr Scalar normLoc(  Belos::NormType type = Belos::TwoNorm ) const {

		Scalar normvec = 0.;

		for( int i=0; i<space()->dim(); ++i )
			normvec =
				(type==Belos::InfNorm)?
					std::max( sFields_[i]->normLoc(type), normvec ):
					( normvec+sFields_[i]->normLoc(type) );

		return( normvec );

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
		Scalar normvec = 0.;

		for( int i=0; i<space()->dim(); ++i )
			normvec += sFields_[i]->normLoc( *weights.sFields_[i] );

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


	/// @}
	/// \name Initialization methods
  /// @{


	/// \brief mv := A
	///
	/// Assign (deep copy) A into mv.
	void assign( const FieldT& a ) {

		for( int i=0; i<space()->dim(); ++i)
			sFields_[i]->assign( *a.sFields_[i] );
		changed();

	}


	/// \brief Replace the vectors with a random vectors.
	///
	/// depending on Fortrans \c Random_number implementation, with always same
	/// seed => not save, if good randomness is required
	void random(bool useSeed = false, int seed = 1) {

		for( int i=0; i<space()->dim(); ++i )
			sFields_[i]->random( useSeed, seed );

		changed();
	}

	/// \brief Replace each element of the vector  with \c alpha.
	void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
		for( int i=0; i<space()->dim(); ++i )
			sFields_[i]->init( alpha );
		changed();
	}


	/// \brief Replace each element of the vector \c getRawPtr[i] with \c alpha[i].
	void init( const Teuchos::Tuple<Scalar,3>& alpha ) {
		for( int i=0; i<space()->dim(); ++i )
			sFields_[i]->init( alpha[i] );
		changed();
	}


	///  \brief initializes VectorField including boundaries to zero 
	void initField() {
		for( int i=0; i<space()->dim(); ++i )
			sFields_[i]->initField();
	}

private:

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
		else {
			const bool& Flow_Type_not_known = true; 
			TEUCHOS_TEST_FOR_EXCEPT( Flow_Type_not_known );
//		case SweptHiemenzFlow: 
//			Scalar kappa = 0.;
//			Scalar sweep_angle = 60.;
//			Scalar pi = 4.*std::atan(1.);
//			VF_init_SHBF( 
//					space()->rankST(),      
//					space()->getProcGrid()->getShift(0),
//					space()->getProcGrid()->getIB(0),
//					space()->nGlo(),
//					space()->nLoc(),
//					space()->bl(),
//					space()->bu(),
//					space()->dl(),
//					space()->du(),
//				 space()->sIndB(U),
//				 space()->eIndB(U),
//				 space()->sIndB(V),
//				 space()->eIndB(V),
//				 space()->sIndB(W),
//				 space()->eIndB(W),
//				 space()->getCoordinatesGlobal()->getX(X,EField::S),
//				 space()->getCoordinatesGlobal()->getX(X,EField::U),
//				 space()->getCoordinatesLocal()->getX(Z,EField::W),
//					space()->getInterpolateV2S()->getC( X ),
//					space()->getDomain()->getDomainSize()->getRe(),
//					0,  // nonDim  
//					kappa, // kappa
//					sweep_angle, // sweep_angle_degrees,  
//					sweep_angle*pi/180., // sweep_angle,          
//					0., // angle_attack,         
//				 sFields_[U]->getRawPtr(),
//				 sFields_[V]->getRawPtr(),
//				 sFields_[W]->getRawPtr() );
//			break;
//		case Disturbance: 
//			VF_init_Dist(
//					space()->rankST(),      
//					space()->nLoc(),
//					space()->bl(),
//					space()->bu(),
//					space()->sIndB(U),
//					space()->eIndB(U),
//					space()->sIndB(V),
//					space()->eIndB(V),
//					space()->sIndB(W),
//					space()->eIndB(W),
//					space()->getDomain()->getBCGlobal()->getBCL( Z ),
//					space()->getCoordinatesLocal()->getX( X, EField::U ),
//					space()->getCoordinatesLocal()->getX( X, EField::S ),
//					space()->getCoordinatesLocal()->getX( Y, EField::S ),
//					space()->getCoordinatesLocal()->getX( Z, EField::W ),
//					space()->getCoordinatesLocal()->getX( Z, EField::S ),
//					3, // dist_type,          
//					0.15, // vortex_ampli_prim,  
//					3., // vortex_x1pos,       
//					3., // vortex_x3pos,       
//					3., // vortex_radius,      
//					10, // vortex_band,        
//					sFields_[U]->getRawPtr(),
//					sFields_[V]->getRawPtr(),
//					sFields_[W]->getRawPtr() );
//			break;
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
				for( int i=0; i<space()->dim(); ++i )
					sFields_[i]->initField( ConstField );
				break;
			}
			case ConstFlow : {
				sFields_[U]->initField( ConstField, para.get<Scalar>( "U", 1.) );
				sFields_[V]->initField( ConstField, para.get<Scalar>( "V", 1.) );
				sFields_[W]->initField( ConstField, para.get<Scalar>( "W", 1.) );
				break;
			}
			case PoiseuilleFlow2D_inX : {
				for( int i=0; i<space()->dim(); ++i )
					if( U==i )
						sFields_[i]->initField( Poiseuille2D_inY );
					else
						sFields_[i]->initField( ConstField );
				break;
			}
			case PoiseuilleFlow2D_inY : {
				for( int i=0; i<space()->dim(); ++i )
					if( V==i )
						sFields_[i]->initField( Poiseuille2D_inX );
					else
						sFields_[i]->initField( ConstField );
				break;
			}
			case PoiseuilleFlow2D_inZ : {
				for( int i=0; i<space()->dim(); ++i )
					if( W==i )
						sFields_[i]->initField( Poiseuille2D_inX );
					else
						sFields_[i]->initField( ConstField );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
				break;
			case Circle2D : {
				sFields_[U]->initField( Grad2D_inY, -1. );
				sFields_[V]->initField( Grad2D_inX,  1. );
				sFields_[W]->initField();
				break;
			}
			case Circle2D_inXZ : {
				sFields_[U]->initField( Grad2D_inY, -1. );
				sFields_[V]->initField();
				sFields_[W]->initField( Grad2D_inX,  1. );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
				break;
			case SweptHiemenzFlow : {
				Scalar kappa       = 0.;
				Scalar sweep_angle = 0.;
				Scalar pi = 4.*std::atan(1.);
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
						space()->getInterpolateV2S()->getC( X ),
						space()->getDomainSize()->getRe(),
						para.get<int>( "nonDim", 0 ),
						para.get<Scalar>( "kappa", 0. ),
						para.get<Scalar>( "seep angle", 0. ),
						para.get<Scalar>( "seep angle", 0. )*pi/180.,
						para.get<Scalar>( "attack angle", 0. ),
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
				break;
			}
			case ScalarFields : {
				sFields_[U]->initField( para.sublist( "U" ) );
				sFields_[V]->initField( para.sublist( "V" ) );
				sFields_[W]->initField( para.sublist( "W" ) );
				break;
			}
		}
		changed();
	}


	/// \brief initializes VectorField with the initial field defined in Fortran
	/// \deprecated 
	void initField( EVectorField flowType,
			Scalar re=0.,
			Scalar om=0.,
			Scalar px = 0.,
			Scalar sca=1. ) {

		switch( flowType ) {
			case ZeroFlow : {
				for( int i=0; i<space()->dim(); ++i )
					sFields_[i]->initField( ConstField );
				break;
			}
			case ConstFlow : {
				sFields_[U]->initField( ConstField, re );
				sFields_[V]->initField( ConstField, om );
				sFields_[W]->initField( ConstField, px );
				break;
			}
			case PoiseuilleFlow2D_inX : {
				for( int i=0; i<space()->dim(); ++i )
					if( U==i )
						sFields_[i]->initField( Poiseuille2D_inY );
					else
						sFields_[i]->initField( ConstField );
				break;
			}
			case PoiseuilleFlow2D_inY : {
				for( int i=0; i<space()->dim(); ++i )
					if( V==i )
						sFields_[i]->initField( Poiseuille2D_inX );
					else
						sFields_[i]->initField( ConstField );
				break;
			}
			case PoiseuilleFlow2D_inZ : {
				for( int i=0; i<space()->dim(); ++i )
					if( W==i )
						sFields_[i]->initField( Poiseuille2D_inX );
					else
						sFields_[i]->initField( ConstField );
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
						re, om, px,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						re, om, px,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						re, om, px,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						re, om, px,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						re,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						re,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						re,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
				break;
			case Circle2D : {
				sFields_[U]->initField( Grad2D_inY, -1. );
				sFields_[V]->initField( Grad2D_inX,  1. );
				sFields_[W]->initField();
				break;
			}
			case Circle2D_inXZ : {
				sFields_[U]->initField( Grad2D_inY, -1. );
				sFields_[V]->initField();
				sFields_[W]->initField( Grad2D_inX,  1. );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						re,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						re, om, px,sca,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						re, om, px,
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
				break;
			case SweptHiemenzFlow : {
				Scalar kappa       = 0.;
				Scalar sweep_angle = 0.;
				Scalar pi = 4.*std::atan(1.);
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
						space()->getCoordinatesGlobal()->getX( ECoord::X, EField::S ),
						space()->getCoordinatesGlobal()->getX( ECoord::X, EField::U ),
						space()->getCoordinatesLocal()->getX( ECoord::Z, EField::W ),
						space()->getInterpolateV2S()->getC( ECoord::X ),
						space()->getDomainSize()->getRe(),
						0,                   // nonDim  
						kappa,               // kappa
						sweep_angle,         // sweep_angle_degrees,  
						sweep_angle*pi/180., // sweep_angle,          
						0.,                  // angle_attack,         
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
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
						space()->getCoordinatesLocal()->getX( ECoord::X, EField::U ),
						space()->getCoordinatesLocal()->getX( ECoord::X, EField::S ),
						space()->getCoordinatesLocal()->getX( ECoord::Y, EField::S ),
						space()->getCoordinatesLocal()->getX( ECoord::Z, EField::W ),
						space()->getCoordinatesLocal()->getX( ECoord::Z, EField::S ),
						3, // dist_type,          
						0.15, // vortex_ampli_prim,  
						3., // vortex_x1pos,       
						3., // vortex_x3pos,       
						3., // vortex_radius,      
						10, // vortex_band,        
						sFields_[U]->getRawPtr(),
						sFields_[V]->getRawPtr(),
						sFields_[W]->getRawPtr() );
				break;
			}
			case ScalarFields : 
			 break;
		}
		changed();
	}


	/// dirty hack(necessary for MG)
  void level() const {}

  /// @}


  /// \brief Print the vector.  To be used for debugging only.
  void print( std::ostream& out=std::cout ) const {
    for(int i=0; i<space()->dim(); ++i) {
      sFields_[i]->print( out );
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
		for( int i=0; i<space()->dim(); ++i )
			getConstFieldPtr(i)->write( count, restart );
	}


public:

	/// \name storage methods.
	/// \brief highly dependent on underlying storage should only be used by
	/// Operator or on top field implementer.
	/// @{ 

	constexpr Ordinal getStorageSize() const { return( sFields_[0]->getStorageSize()*3 ); }

  void setStoragePtr( Scalar*  array ) {
    s_ = array;
    Ordinal n = sFields_[0]->getStorageSize();
    for( int i=0; i<3; ++i )
      sFields_[i]->setStoragePtr( s_+i*n );
  }

	constexpr Scalar* getRawPtr() { return( s_ ); }

	constexpr const Scalar* getConstRawPtr() const { return( s_ ); }

  constexpr Scalar* getRawPtr ( int i )       { return( sFields_[i]->getRawPtr() ); }

	constexpr const Scalar* getConstRawPtr ( int i )  const  { return( sFields_[i]->getConstRawPtr() ); }

	///  @} 

  constexpr Teuchos::RCP<SF> getFieldPtr( int i ) { return(  sFields_[i] ); }
  constexpr SF& getField   ( int i ) { return( *sFields_[i] ); }

  constexpr Teuchos::RCP<const SF> getConstFieldPtr( int i ) const { return(  sFields_[i] ); }
  constexpr const SF&  getConstField   ( int i ) const { return( *sFields_[i] ); }

  constexpr const Teuchos::RCP<const SpaceT>& space() const { return( AbstractField<SpaceT>::space_ ); }

  constexpr const MPI_Comm& comm() const { return(space()->comm()); }


  /// \name comunication methods.
  /// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
  ///
  /// @{
	
  void changed( const int& vel_dir, const int& dir ) const {
    getConstFieldPtr( vel_dir )->changed( dir );
  }

  void changed() const {
    for( int vel_dir=0; vel_dir<space()->dim(); ++vel_dir )
      for( int dir=0; dir<space()->dim(); ++dir ) {
        changed( vel_dir, dir );
      }
  }

  void exchange( const int& vel_dir, const int& dir ) const {
    getConstFieldPtr(vel_dir)->exchange(dir);
  }

  void exchange() const {
    for( int vel_dir=0; vel_dir<space()->dim(); ++vel_dir )
      for( int dir=0; dir<space()->dim(); ++dir )
        exchange( vel_dir, dir );
  }

  void setExchanged( const int& vel_dir, const int& dir ) const {
    getConstFieldPtr( vel_dir )->setExchanged( dir );
  }

	void setExchanged() const {
		for( int vel_dir=0; vel_dir<space()->dim(); ++vel_dir )
			for( int dir=0; dir<space()->dim(); ++dir ) {
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
