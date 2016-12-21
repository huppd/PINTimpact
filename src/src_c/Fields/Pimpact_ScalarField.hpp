#pragma once
#ifndef PIMPACT_SCALARFIELD_HPP
#define PIMPACT_SCALARFIELD_HPP


#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Tuple.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_AbstractField.hpp"
#include "Pimpact_extern_ScalarField.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact {



/// \brief important basic Vector class
/// vector for a scalar field, e.g.: pressure,
/// \ingroup Field
/// \todo make owning_ fType_ template parameter
template<class SpaceType>
class ScalarField : private AbstractField< SpaceType > {

public:

  using SpaceT = SpaceType;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

	//static const int sdim = SpaceT::sdim;

  using ScalarArray = Scalar*;
  using FieldT = ScalarField< SpaceT >;
  using State = Teuchos::Tuple<bool,3>;


  ScalarArray s_;

  bool owning_;

  State exchangedState_;

  const EField fType_;

	void allocate() {
		Ordinal n = getStorageSize();
		s_ = new Scalar[n];
	}

public:

  ScalarField( const Teuchos::RCP<const SpaceT>& space, bool owning=true, EField fType=EField::S ):
    AbstractField<SpaceT>( space ),
    owning_(owning),
    exchangedState_( Teuchos::tuple( true, true, true ) ),
    fType_(fType) {

    if( owning_ ) {
			allocate();
			initField();
    }
  };

  /// \brief copy constructor.
  ///
	/// \note copyType is 
  /// \param sF ScalarField which is copied
  /// \param copyType by default a ECopy::Deep is done but also allows to ECopy::Shallow
	ScalarField( const ScalarField& sF, ECopy copyType=ECopy::Deep ):
		AbstractField<SpaceT>( sF.space() ),
		owning_( sF.owning_ ),
		exchangedState_( sF.exchangedState_ ),
		fType_( sF.fType_ ) {

    if( owning_ ) {

			allocate();

			switch( copyType ) {
				case ECopy::Shallow:
					initField();
					break;
				case ECopy::Deep:
					*this = sF;
					break;
			}
		}
		};


	~ScalarField() { if( owning_ ) delete[] s_; }


  Teuchos::RCP<FieldT> clone( ECopy copyType=ECopy::Deep ) const {

		Teuchos::RCP<FieldT> mv = Teuchos::rcp( new FieldT( space(), true, this->fType_ ) );

		switch( copyType ) {
			case ECopy::Shallow:
				break;
			case ECopy::Deep:
				*mv = *this;
				break;
		}
		return( mv );
	}

  /// \name Attribute methods
  /// \{

  /// \brief returns the length of Field.
	constexpr Ordinal getLength() const {

		Teuchos::RCP<const BoundaryConditionsGlobal<SpaceT::dimension> > bc =
			space()->getBCGlobal();

		Ordinal vl = 1;

		for( int dir = 0; dir<SpaceT::sdim; ++dir ) {
			vl *= space()->nGlo(dir) +
				( (PeriodicBC==bc->getBCL(dir))?
					-1:
					(fType_==dir)?1:0);
		}

		return( vl );
	}



  /// @}
  /// \name Update methods
  /// @{

  /// \brief Replace \c this with \f$\alpha A + \beta B\f$.
	/// \todo make checks for spaces and k
	void add( const Scalar& alpha, const FieldT& A, const Scalar& beta, const
			FieldT& B, const With& wB=With::B ) {

#ifndef NDEBUG
		for( int dir=0; dir<3; ++dir ) {
			bool same_space = space()->nLoc(dir)>A.space()->nLoc(dir) || 
				space()->nLoc(dir)>B.space()->nLoc(dir);
			assert( !same_space );
			bool consistent_space = (
					(A.space()->nLoc(dir)-1)%(space()->nLoc(dir)-1) )!=0 || (
					(B.space()->nLoc(dir)-1)%(space()->nLoc(dir)-1) )!=0 ;
			assert( !consistent_space );
		}
#endif
		Teuchos::Tuple<Ordinal,3> da;
		Teuchos::Tuple<Ordinal,3> db;

		for( int dir=0; dir<3; ++dir ) {
			da[dir] = ( A.space()->nLoc(dir)-1 )/( space()->nLoc(dir)-1 );
			db[dir] = ( B.space()->nLoc(dir)-1 )/( space()->nLoc(dir)-1 );
		}

		for( Ordinal k=space()->begin(fType_,Z,wB); k<=space()->end(fType_,Z,wB); ++k )
			for( Ordinal j=space()->begin(fType_,Y,wB); j<=space()->end(fType_,Y,wB); ++j )
				for( Ordinal i=space()->begin(fType_,X,wB); i<=space()->end(fType_,X,wB); ++i )
					at(i,j,k) = alpha*A.at( (i-1)*da[0]+1, (j-1)*da[1]+1,(k-1)*da[2]+1 )
						         + beta*B.at( (i-1)*db[0]+1, (j-1)*db[1]+1,(k-1)*db[2]+1 );

		changed();
	}


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
	void abs( const FieldT& y, const With& bcYes=With::B ) {

		for( int dir=0; dir<3; ++dir )
			assert( space()->nLoc(dir)==y.space()->nLoc(dir) );

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					at(i,j,k) = std::fabs( y.at(i,j,k) );

		changed();
	}


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal( const FieldT& y, const With& bcYes=With::B ) {

#ifndef NDEBUG
		for( int dir=0; dir<3; ++dir ) {
			bool same_space = space()->nLoc(dir)!=y.space()->nLoc(dir);
			assert( !same_space );
		}
#endif

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					at(i,j,k) = Teuchos::ScalarTraits<Scalar>::one()/ y.at(i,j,k);

    changed();
  }


  /// \brief Scale each element of the vector with \c alpha.
	void scale( const Scalar& alpha, const With& bcYes=With::B ) {

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					at(i,j,k) *= alpha;

		changed();
	}


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot y_i \quad \mbox{for } i=1,\dots,n \f]
	void scale( const FieldT& y, const With& bcYes=With::B ) {

#ifndef NDEBUG
		for( int dir=0; dir<3; ++dir ) {
			bool same_space = space()->nLoc(dir)!=y.space()->nLoc(dir);
			assert( !same_space );
		}
#endif

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					at(i,j,k) *= y.at(i,j,k);
		changed();
	}

  /// @}
  /// \name Norm method(reductions)
  /// @{

	/// \brief Compute a local scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
	constexpr Scalar dotLoc( const FieldT& y, const With& bcYes=With::B ) const {

#ifndef NDEBUG
		for( int dir=0; dir<3; ++dir ) {
			bool same_space = space()->nLoc(dir)!=y.space()->nLoc(dir);
			assert( !same_space );
		}
#endif

		Scalar b = Teuchos::ScalarTraits<Scalar>::zero();

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					b += at(i,j,k)*y.at(i,j,k);

		return( b );
	}

	/// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y
	/// and \c this, i.e.\f$b = y^H this\f$.
	constexpr Scalar dot( const FieldT& y, const With& bcYes=With::B ) const {
		return( this->reduce( comm(), dotLoc( y, bcYes ) ) );
	}

  constexpr Scalar normLoc1( const With& bcYes=With::B ) const {

    Scalar normvec = Teuchos::ScalarTraits<Scalar>::zero();

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					normvec += std::fabs( at(i,j,k) );

    return( normvec );
  }


  constexpr Scalar normLoc2( const With& bcYes=With::B ) const {

    Scalar normvec = Teuchos::ScalarTraits<Scalar>::zero();

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					normvec += std::pow( at(i,j,k), 2 );

    return( normvec );
  }


  constexpr Scalar normLocInf( const With& bcYes=With::B ) const {

    Scalar normvec = Teuchos::ScalarTraits<Scalar>::zero();

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					normvec = std::fmax( std::fabs(at(i,j,k)), normvec );

    return( normvec );
  }

	constexpr Scalar normLoc( Belos::NormType type = Belos::TwoNorm, const With& bcYes=With::B ) const {

		return(
				( Belos::OneNorm==type)?
					normLoc1(bcYes):
					(Belos::TwoNorm==type)?
						normLoc2(bcYes):
						normLocInf(bcYes) );
	}


  /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  constexpr Scalar norm( Belos::NormType type = Belos::TwoNorm, const With& bcYes=With::B ) const {

		Scalar normvec = this->reduce(
				comm(),
				normLoc( type,bcYes ),
				(Belos::InfNorm==type)?MPI_MAX:MPI_SUM );

		normvec =
			(Belos::TwoNorm==type) ?
				std::sqrt(normvec) :
				normvec;

    return( normvec );
  }


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr Scalar normLoc( const FieldT& weights, const With& bcYes=With::B ) const {

		for( int dir=0; dir<3; ++dir )
			assert( space()->nLoc(dir)==weights.space()->nLoc(dir) );

    Scalar normvec = Teuchos::ScalarTraits<Scalar>::zero();

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					normvec += at(i,j,k)*at(i,j,k)*weights.at(i,j,k)*weights.at(i,j,k);

    return( normvec );
  }


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr Scalar norm( const FieldT& weights, const With& bcYes=With::B ) const {
		return( std::sqrt( this->reduce( comm(), normLoc( weights, bcYes ) ) ) );
	}


  //\}
  /// \name Initialization methods
  //\{

  /// \brief *this := a 
  ///
  /// Assign (deep copy) \c a into \c this.
  /// total deep, boundaries and everything.
  /// \note the \c StencilWidths is not take care of assuming every field is generated with one
	ScalarField& operator=( const ScalarField& a ) {

		assert( getType()==a.getType() );
		assert( getStorageSize()==a.getStorageSize() );

		std::copy_n( a.s_, getStorageSize(), s_ );

		for( int dir=0; dir<SpaceT::sdim; ++dir )
			exchangedState_[dir] = a.exchangedState_[dir];
		return *this;
	}

  /// \brief Replace the vectors with a random vectors.
  /// Depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is required
  void random( bool useSeed = false, const With& bcYes=With::B , int seed = 1 ) {

		std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis( -0.5, 0.5 );

		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					at(i,j,k) = dis(gen);

		if( !space()->getProcGrid()->participating() )
			for( Ordinal k=space()->begin(fType_,Z,With::B); k<=space()->end(fType_,Z,With::B); ++k )
				for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j )
					for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i )
						at(i,j,k) = Teuchos::ScalarTraits<Scalar>::zero();
		changed();
  }


  /// \brief Replace each element of the vector  with \c alpha.
	/// \param alpha init value
	/// \param bcYes also initializing the boundary values
	void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero(), const With& bcYes=With::B ) {

		if( With::B==bcYes ){
			std::fill_n( s_, getStorageSize(), alpha );
			exchangedState_[X] = true;
			exchangedState_[Y] = true;
			exchangedState_[Z] = true;
		}
		else{
			for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
				for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
					for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
						at(i,j,k) = alpha;

			if( !space()->getProcGrid()->participating() )
				for( Ordinal k=space()->begin(fType_,Z,With::B); k<=space()->end(fType_,Z,With::B); ++k )
					for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j )
						for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i )
							at(i,j,k) = Teuchos::ScalarTraits<Scalar>::zero();
			changed();
		}
	}

protected:

	/// \brief helper function getting \c EScalarField for switch statement from name 
	///
	/// \param name input name
	/// \return according int number
	EScalarField string2enum( const std::string& name ) {

		std::string lcName = name;
		std::transform(lcName.begin(), lcName.end(), lcName.begin(), ::tolower);

		if( "constant" == lcName) return( ConstField );
		else if( "grad in x" == lcName ) return( Grad2D_inX );
		else if( "grad in y" == lcName ) return( Grad2D_inY );
		else if( "grad in z" == lcName ) return( Grad2D_inZ );
		else if( "poiseuille" == lcName ) return( Poiseuille2D_inX );
		else if( "poiseuille in x" == lcName ) return( Poiseuille2D_inX );
		else if( "poiseuille in y" == lcName ) return( Poiseuille2D_inY );
		else if( "poiseuille in z" == lcName ) return( Poiseuille2D_inZ );
		else if( "point" == lcName ) return( FPoint );
		else {
			const bool& Flow_Type_not_known = false; 
			assert( Flow_Type_not_known );
		}
		return( ConstField ); // just to please the compiler
	}

public:

	/// \brief initializes field from a lambda function
	///
	/// \tparam Functor need a to have an Operator (x,y,z)->u
	/// \param func  note that x,y,z, should be defined between 0 and one the scaling happens here
	template<typename Functor>
	void initFromFunction( Functor&& func ) {

		Teuchos::RCP<const CoordinatesLocal<Scalar,Ordinal,SpaceT::dimension,SpaceT::dimNC> > coord =
			space()->getCoordinatesLocal();
		Teuchos::RCP<const DomainSize<Scalar,SpaceT::sdim> > domain = space()->getDomainSize();

		const With& bcYes = With::B;
		for( Ordinal k=space()->begin(fType_,Z,bcYes); k<=space()->end(fType_,Z,bcYes); ++k )
			for( Ordinal j=space()->begin(fType_,Y,bcYes); j<=space()->end(fType_,Y,bcYes); ++j )
				for( Ordinal i=space()->begin(fType_,X,bcYes); i<=space()->end(fType_,X,bcYes); ++i )
					at(i,j,k) = func(
							( coord->getX(fType_,X,i)-domain->getOrigin(X) )/domain->getSize(X),
							( coord->getX(fType_,Y,j)-domain->getOrigin(Y) )/domain->getSize(Y),
							( coord->getX(fType_,Z,k)-domain->getOrigin(Z) )/domain->getSize(Z) );
	}


	///  \brief initializes including boundaries to zero 
	void initField( Teuchos::ParameterList& para ) {

		EScalarField type =
			string2enum( para.get<std::string>( "Type", "constant" ) );

		switch( type ) {
			case ConstField :
				initFromFunction( [&para](Scalar,Scalar,Scalar)->const Scalar&{
						return( para.get<Scalar>( "C", Teuchos::ScalarTraits<Scalar>::zero() ) ); } );
				break;
			case Grad2D_inX :
				SF_init_2DGradX(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( X ),
						space()->getCoordinatesLocal()->getX( X, fType_ ),
						s_,
						para.get<Scalar>( "dx", Teuchos::ScalarTraits<Scalar>::one() )	);
				break;
			case Grad2D_inY :
				SF_init_2DGradY(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( Y ),
						space()->getCoordinatesLocal()->getX( Y, fType_ ),
						s_ ,
						para.get<Scalar>( "dy", Teuchos::ScalarTraits<Scalar>::one() )	);
				break;
			case Grad2D_inZ :
				SF_init_2DGradZ(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( Z ),
						space()->getCoordinatesLocal()->getX( Z, fType_ ),
						s_ ,
						para.get<Scalar>( "dz", Teuchos::ScalarTraits<Scalar>::one() )	);
				break;
			case Poiseuille2D_inX :
				SF_init_2DPoiseuilleX(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( X ),
						space()->getCoordinatesLocal()->getX( X, fType_ ),
						s_ );
				break;
			case Poiseuille2D_inY :
				SF_init_2DPoiseuilleY(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( Y ),
						space()->getCoordinatesLocal()->getX( Y, fType_ ),
						s_ );
				break;
			case Poiseuille2D_inZ :
				SF_init_2DPoiseuilleZ(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( Z ),
						space()->getCoordinatesLocal()->getX( Z, fType_ ),
						s_ );
				break;
			case FPoint :
				Scalar xc[3] = { 
					para.get<Scalar>( "c_x", Teuchos::ScalarTraits<Scalar>::one() ),
					para.get<Scalar>( "c_y", space()->getDomainSize()->getSize( Y )/2. ),
					para.get<Scalar>( "c_z", space()->getDomainSize()->getSize( Z )/2. ) };
				Scalar amp = para.get<Scalar>( "amp", Teuchos::ScalarTraits<Scalar>::one() );
				Scalar sig[3] = {
					para.get<Scalar>( "sig_x", 0.2 ),
					para.get<Scalar>( "sig_y", 0.2 ),
					para.get<Scalar>( "sig_z", 0.2 ) };
				SF_init_Vpoint(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getCoordinatesLocal()->getX( X, fType_ ),
						space()->getCoordinatesLocal()->getX( Y, fType_ ),
						space()->getCoordinatesLocal()->getX( Z, fType_ ),
						xc,
						amp,
						sig,
						s_ );
				break;
		}

		if( !space()->getProcGrid()->participating() ) // not sure why?
			initFromFunction( [](Scalar,Scalar,Scalar)->Scalar{
					return( Teuchos::ScalarTraits<Scalar>::zero() ); } );

		changed();
	}


	/// \brief initializes VectorField with the initial field defined in Fortran
	/// \deprecated
	void initField( EScalarField fieldType = ConstField, Scalar alpha=Teuchos::ScalarTraits<Scalar>::zero() ) {

		switch( fieldType ) {
			case ConstField :
				std::fill_n( s_, getStorageSize(), alpha );
				//initFromFunction( [&alpha](Scalar,Scalar,Scalar)->Scalar&{
						//return( alpha ); } );
				break;
			case Grad2D_inX :
				SF_init_2DGradX(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( X ),
						space()->getCoordinatesLocal()->getX( X, fType_ ),
						s_,
						(std::fabs(alpha)<Teuchos::ScalarTraits<Scalar>::eps())?Teuchos::ScalarTraits<Scalar>::one():alpha	);
				break;
			case Grad2D_inY :
				SF_init_2DGradY(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( Y ),
						space()->getCoordinatesLocal()->getX( Y, fType_ ),
						s_ ,
						(std::fabs(alpha)<Teuchos::ScalarTraits<Scalar>::eps())?Teuchos::ScalarTraits<Scalar>::one():alpha	);
				break;
			case Grad2D_inZ :
				SF_init_2DGradZ(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( Z ),
						space()->getCoordinatesLocal()->getX( Z, fType_ ),
						s_ ,
						(std::fabs(alpha)<Teuchos::ScalarTraits<Scalar>::eps())?Teuchos::ScalarTraits<Scalar>::one():alpha	);
				break;
			case Poiseuille2D_inX :
				SF_init_2DPoiseuilleX(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( X ),
						space()->getCoordinatesLocal()->getX( X, fType_ ),
						s_ );
				break;
			case Poiseuille2D_inY :
				SF_init_2DPoiseuilleY(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( Y ),
						space()->getCoordinatesLocal()->getX( Y, fType_ ),
						s_ );
				break;
			case Poiseuille2D_inZ :
				SF_init_2DPoiseuilleZ(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getDomainSize()->getSize( Z ),
						space()->getCoordinatesLocal()->getX( Z, fType_ ),
						s_ );
				break;
			case FPoint :
				Scalar xc[3] =
				{ 
					Teuchos::ScalarTraits<Scalar>::one(),
					//				space()->getDomainSize()->getSize( X )/4.,
					space()->getDomainSize()->getSize( Y )/2.,
					space()->getDomainSize()->getSize( Z )/2. };
				Scalar amp = alpha; //2./space()->getDomainSize()->getRe();
				Scalar sig[3] = { 0.2, 0.2, 0.2 };
				SF_init_Vpoint(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						space()->getCoordinatesLocal()->getX( X, fType_ ),
						space()->getCoordinatesLocal()->getX( Y, fType_ ),
						space()->getCoordinatesLocal()->getX( Z, fType_ ),
						xc,
						amp,
						sig,
						s_ );
				break;
		}

		if( !space()->getProcGrid()->participating() )
			initFromFunction( [&alpha](Scalar,Scalar,Scalar)->Scalar{
					return( Teuchos::ScalarTraits<Scalar>::zero() ); } );
		changed();
	}


	/// \brief extrapolate the velocity points outside the domain such that the
	///  interpolated value on the boundary is zero
	///
	/// \param trans transposed
	void extrapolateBC( const Belos::ETrans& trans=Belos::NOTRANS ) {

		switch( trans ) {
			case( Belos::NOTRANS ): {

				switch( fType_ ) {
					case( U ):  {
						if( space()->getBCLocal()->getBCL(X) > 0 ) {
							Ordinal i = space()->begin(fType_,X,With::B);
							for( Ordinal k=space()->begin(fType_,Z, With::B); k<=space()->end(fType_,Z,With::B); ++k )
								for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j ) {
									at(i,j,k) = 0.;
									for( Ordinal ii=0; ii<=space()->du(X); ++ii )
										at(i,j,k) -= at(1+ii,j,k)*space()->getInterpolateV2S()->getC(X,1,ii)/space()->getInterpolateV2S()->getC(X,1,-1);
								}
						}
						if( space()->getBCLocal()->getBCU(X) > 0 ) {
							Ordinal i = space()->end(fType_,X,With::B);
							for( Ordinal k=space()->begin(fType_,Z, With::B); k<=space()->end(fType_,Z,With::B); ++k )
								for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j ) {
									at(i,j,k) = 0.;
									for( Ordinal ii=space()->dl(X); ii<=-1; ++ii )
										at(i,j,k) -= space()->getInterpolateV2S()->getC(X,i,ii)*at(i+ii,j,k)/space()->getInterpolateV2S()->getC(X,i,0);
								}
						}
						break;
					}
					case( V ) : {
						if( space()->getBCLocal()->getBCL(Y) > 0 ) {
							Ordinal j = space()->begin(fType_,Y,With::B);
							for( Ordinal k=space()->begin(fType_,Z, With::B); k<=space()->end(fType_,Z,With::B); ++k )
								for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i ) {
									at(i,j,k) = 0.;
									for( Ordinal jj=0; jj<=space()->du(Y); ++jj )
										at(i,j,k) -= at(i,1+jj,k)*space()->getInterpolateV2S()->getC(Y,1,jj)/space()->getInterpolateV2S()->getC(Y,1,-1);  
								}
						}
						if( space()->getBCLocal()->getBCU(Y) > 0 ) {
							Ordinal j = space()->end(fType_,Y,With::B);
							for( Ordinal k=space()->begin(fType_,Z, With::B); k<=space()->end(fType_,Z,With::B); ++k )
								for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i ) {
									at(i,j,k) = 0.;
									for( Ordinal jj=space()->dl(Y); jj<=-1; ++jj )
										at(i,j,k) -= space()->getInterpolateV2S()->getC(Y,j,jj)*at(i,j+jj,k)/space()->getInterpolateV2S()->getC(Y,j,0);
								}
						}
						break;
					}
					case( W ) : {
						if( space()->getBCLocal()->getBCL(Z) > 0 ) {
							Ordinal k = space()->begin(fType_,Z,With::B);
							for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j )
								for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i ) {
									at(i,j,k) = 0.;
									for( Ordinal kk=0; kk<=space()->du(Z); ++kk )
										at(i,j,k) -= space()->getInterpolateV2S()->getC(Z,1,kk)*at(i,j,1+kk)/space()->getInterpolateV2S()->getC(Z,1,-1);  
								}
						}
						if( space()->getBCLocal()->getBCU(Z) > 0 ) {
							Ordinal k = space()->end(fType_,Z,With::B);
							for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j )
								for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i ) {
									at(i,j,k) = 0.;
									for( Ordinal kk=space()->dl(Z); kk<=-1; ++kk )
										at(i,j,k) -= space()->getInterpolateV2S()->getC(Z,k,kk)*at(i,j,k+kk)/space()->getInterpolateV2S()->getC(Z,k,0);
								}
						}
						break;
					}
					case( S ) : break;
				}
				break;
			}
			case Belos::TRANS : {

				switch( fType_ ) {
					case( U ) : {
						if( space()->getBCLocal()->getBCL(X) > 0 ) {
							Ordinal i = space()->begin(fType_,X,With::B);
							for( Ordinal k=space()->begin(fType_,Z, With::B); k<=space()->end(fType_,Z,With::B); ++k )
								for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j ) {
									for( Ordinal ii=0; ii<=space()->du(X); ++ii )
										at(i+ii+1,j,k) -= at(i,j,k)*space()->getInterpolateV2S()->getC(X,1,ii)/space()->getInterpolateV2S()->getC(X,1,-1);  
									at(i,j,k) = 0.;
								}
						}
						if( space()->getBCLocal()->getBCU(X) > 0 ) {
							Ordinal i = space()->end(fType_,X,With::B);
							for( Ordinal k=space()->begin(fType_,Z, With::B); k<=space()->end(fType_,Z,With::B); ++k )
								for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j ) {
									for( Ordinal ii=space()->dl(X); ii<=-1; ++ii )
										at(i+ii,j,k) -= space()->getInterpolateV2S()->getC(X,i,ii)*at(i,j,k)/space()->getInterpolateV2S()->getC(X,i,0); 
									at(i,j,k) = 0.;
								}
						}
						break;
					}
					case( V ) : {
						if( space()->getBCLocal()->getBCL(Y) > 0 ) {
							Ordinal j = space()->begin(fType_,Y,With::B);
							for( Ordinal k=space()->begin(fType_,Z, With::B); k<=space()->end(fType_,Z,With::B); ++k )
								for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i ) {
									for( Ordinal jj=0; jj<=space()->du(Y); ++jj )
										at(i,j+jj+1,k) -= at(i,j,k)*space()->getInterpolateV2S()->getC(Y,1,jj)/space()->getInterpolateV2S()->getC(Y,1,-1);  
									at(i,j,k) = 0.;
								}
						}
						if( space()->getBCLocal()->getBCU(Y) > 0 ) {
							Ordinal j = space()->end(fType_,Y,With::B);
							for( Ordinal k=space()->begin(fType_,Z, With::B); k<=space()->end(fType_,Z,With::B); ++k )
								for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i ) {
									for( Ordinal jj=space()->dl(Y); jj<=-1; ++jj )
										at(i,j+jj,k) -= space()->getInterpolateV2S()->getC(Y,j,jj)*at(i,j,k)/space()->getInterpolateV2S()->getC(Y,j,0); 
									at(i,j,k) = 0.;
								}
						}
						break;
					}
					case( W ) : {
						if( space()->getBCLocal()->getBCL(Z) > 0 ) {
							Ordinal k = space()->begin(fType_,Z,With::B);
							for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j )
								for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i ) {
									at(i,j,k) /= space()->getInterpolateV2S()->getC(Z,1,-1);
									for( Ordinal kk=0; kk<=space()->du(Z); ++kk )
										at(i,j,k+kk+1) -= space()->getInterpolateV2S()->getC(Z,1,kk)*at(i,j,k);  
									at(i,j,k) = 0.;
								}
						}
						if( space()->getBCLocal()->getBCU(Z) > 0 ) {
							Ordinal k = space()->end(fType_,Z,With::B);
							for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j )
								for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i ) {
									at(i,j,k) /= space()->getInterpolateV2S()->getC(Z,k,0);
									for( Ordinal kk=space()->dl(Z); kk<=-1; ++kk )
										at(i,j,k+kk) -= space()->getInterpolateV2S()->getC(Z,k,kk)*at(i,j,k); 
									at(i,j,k) = 0.;
								}
						}
						break;
					}
					case( S ) : break;
				}
			}
			case Belos::CONJTRANS : break;
		}
	}


	/// \brief levels field if scalar field
	void level() {

		if( EField::S == fType_ ) {

			Scalar pre0 = Teuchos::ScalarTraits<Scalar>::zero();

			for( Ordinal k=space()->begin(fType_,Z); k<=space()->end(fType_,Z); ++k )
				for( Ordinal j=space()->begin(fType_,Y); j<=space()->end(fType_,Y); ++j )
					for( Ordinal i=space()->begin(fType_,X); i<=space()->end(fType_,X); ++i )
						pre0 += at(i,j,k);

			pre0 = this->reduce( space()->comm(), pre0 );
			pre0 /= static_cast<Scalar>( getLength() );

			for( Ordinal k=space()->begin(fType_,Z); k<=space()->end(fType_,Z); ++k )
				for( Ordinal j=space()->begin(fType_,Y); j<=space()->end(fType_,Y); ++j )
					for( Ordinal i=space()->begin(fType_,X); i<=space()->end(fType_,X); ++i )
						at(i,j,k) -= pre0;
		}
	}

  /// \}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& out=std::cout )  const {

    out << "--- FieldType: " << fType_ << "--- \n";
    out << "--- StorageSize: " << getStorageSize() << "---\n";
    out << "--- owning: " << owning_ << "---\n";
    out << "--- exchangedState: " << exchangedState_ << "--\n\n";
		out << "i,\tj,\tk,\tphi(i,j,k)\n";

		Teuchos::Tuple<Ordinal,3> cw;
		for(int i=0; i<3; ++i)
			cw[i] = space()->nLoc(i) + space()->bu(i) - space()->bl(i) + 1;

		for( Ordinal k=space()->begin(fType_,Z,With::B); k<=space()->end(fType_,Z,With::B); ++k )
			for( Ordinal j=space()->begin(fType_,Y,With::B); j<=space()->end(fType_,Y,With::B); ++j )
				for( Ordinal i=space()->begin(fType_,X,With::B); i<=space()->end(fType_,X,With::B); ++i )
					out << i << "\t" << j << "\t" << k << "\t" << at(i,j,k) << "\n";
  }


  /// Write the ScalarField to an hdf5 file, the velocities are interpolated to the pressure points
  /// \todo add restart
	void write( int count=0 , bool restart=false ) const {

		if( 0==space()->rankS() )
			switch(fType_) {
				case U:
					std::cout << "writing velocity field x(" << count << ") ...\n";
					break;
				case V:
					std::cout << "writing velocity field y(" << count << ") ...\n";
					break;
				case W:
					std::cout << "writing velocity field z(" << count << ") ...\n";
					break;
				case EField::S:
					std::cout << "writing pressure field  (" << count << ") ...\n";
					Teuchos::Tuple<Ordinal,3> N;
					for( int i=0; i<3; ++i ) {
						N[i] = space()->nGlo(i);
						if( space()->getBCGlobal()->getBCL(i)==Pimpact::PeriodicBC )
							N[i] = N[i]-1;
					}
					std::ofstream xfile;
					std::ostringstream ss;
					ss << std::setw( 5 ) << std::setfill( '0' ) << count;
					std::string fname = "pre_"+ss.str();
					xfile.open( fname+".xmf", std::ofstream::out );
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
					xfile << "\t\t\t\t\t" << fname << ".h5:/VectorX\n";
					xfile << "\t\t\t\t</DataItem>\n";
					xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
					xfile << "\t\t\t\t\tDimensions=\""<< N[1] << "\"\n";
					xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
					xfile << "\t\t\t\t\tPrecision=\"8\"\n";
					xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
					xfile << "\t\t\t\t\t" << fname << ".h5:/VectorY\n";
					xfile << "\t\t\t\t</DataItem>\n";
					xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
					xfile << "\t\t\t\t\tDimensions=\""<< N[2] << "\"\n";
					xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
					xfile << "\t\t\t\t\tPrecision=\"8\"\n";
					xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
					xfile << "\t\t\t\t\t" << fname << ".h5:/VectorZ\n";
					xfile << "\t\t\t\t</DataItem>\n";
					xfile << "\t\t\t</Geometry>\n";
					xfile << "\t\t\t<Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n";
					xfile << "\t\t\t\t<DataItem Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
					xfile << "\t\t\t\t\t" << fname << ".h5:/pre\n";
					xfile << "\t\t\t\t</DataItem>\n";
					xfile << "\t\t\t</Attribute>\n";
					xfile << "\t\t</Grid>\n";
					xfile << "\t</Domain>\n";
					xfile << "</Xdmf>\n";
					xfile.close();
					break;
			}

		if( !restart ) {
			Teuchos::RCP< ScalarField<SpaceT> > temp;

			if( EField::S != fType_ ) {
				temp = Teuchos::rcp(
						new ScalarField<SpaceT>( space(), true, EField::S ) );
				space()->getInterpolateV2S()->apply( *this, *temp );
			}

			if( 2==SpaceT::sdim ) {

				write_hdf5_2D(
            space()->rankS(),
            MPI_Comm_c2f( space()->comm() ),
            space()->nGlo(),
            space()->getBCGlobal()->getBCL(),
            space()->getBCGlobal()->getBCU(),
            space()->nLoc(),
            space()->bl(),
            space()->bu(),
            space()->sInd(EField::S),
            space()->eInd(EField::S),
            space()->getStencilWidths()->getLS(),
            space()->np(),
            space()->ib(),
            space()->getShift(),
            (int)fType_,
            count,
            (EField::S==fType_)?9:10,
						(EField::S==fType_)?s_:temp->s_,
						space()->getCoordinatesGlobal()->getX(0,EField::S),
						space()->getCoordinatesGlobal()->getX(1,EField::S),
						space()->getDomainSize()->getRe(),
						space()->getDomainSize()->getAlpha2() );
      }
      else if( 3==SpaceT::sdim ) {

        int stride[3] = {1,1,1};

        write_hdf_3D(
            space()->rankS(),
            MPI_Comm_c2f( space()->comm() ),
            space()->nGlo(),
            space()->getBCGlobal()->getBCL(),
            space()->getBCGlobal()->getBCU(),
            space()->nLoc(),
            space()->bl(),
            space()->bu(),
            space()->sInd(EField::S),
            space()->eInd(EField::S),
            space()->getStencilWidths()->getLS(),
            space()->np(),
            space()->ib(),
            space()->getShift(),
						(int)fType_+1,
						(int)EField::S+1,
            count,
            (EField::S==fType_)?9:10,
            stride,
            (EField::S==fType_)?s_:temp->s_,
            space()->getCoordinatesGlobal()->getX(0,EField::S),
            space()->getCoordinatesGlobal()->getX(1,EField::S),
            space()->getCoordinatesGlobal()->getX(2,EField::S),
            space()->getCoordinatesGlobal()->getX(0,EField::U),
            space()->getCoordinatesGlobal()->getX(1,EField::V),
            space()->getCoordinatesGlobal()->getX(2,EField::W),
            space()->getDomainSize()->getRe(),
            space()->getDomainSize()->getAlpha2() );

      }
    }
    else {

      int stride[3] = {1,1,1};

      write_hdf_3D(
          space()->rankS(),
          MPI_Comm_c2f( space()->comm() ),
          space()->nGlo(),
          space()->getBCGlobal()->getBCL(),
          space()->getBCGlobal()->getBCU(),
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType_),
          space()->eInd(fType_),
          space()->getStencilWidths()->getLS(),
          space()->np(),
          space()->ib(),
          space()->getShift(),
          (int)fType_+1,
          (int)fType_+1,
          count,
          (EField::S==fType_)?9:10,
          stride,
          s_,
          space()->getCoordinatesGlobal()->getX(0,EField::S),
          space()->getCoordinatesGlobal()->getX(1,EField::S),
          space()->getCoordinatesGlobal()->getX(2,EField::S),
          space()->getCoordinatesGlobal()->getX(0,EField::U),
          space()->getCoordinatesGlobal()->getX(1,EField::V),
          space()->getCoordinatesGlobal()->getX(2,EField::W),
          space()->getDomainSize()->getRe(),
          space()->getDomainSize()->getAlpha2() );
    }
  }



public:

  constexpr const EField& getType() const { return( fType_ ); }

   /// \name storage methods.
	 /// \brief highly dependent on underlying storage should only be used by
	 /// Operator or on top field implementer.  
   /// @{

  Ordinal getStorageSize() const {

    Ordinal n = 1;
    for(int i=0; i<3; ++i)
      n *= space()->nLoc(i)+space()->bu(i)-space()->bl(i)+1; // seems wrong: there a one was added for AMG, but it is not neede error seem to be in Impact there it should be (B1L+1:N1+B1U) probably has to be changed aganin for 3D

    return( n );
  }

	void setStoragePtr( Scalar*  array ) { s_ = array; }

	ScalarArray getRawPtr() { return( s_ ); }

	const Scalar* getConstRawPtr() const { return( s_ ); }


  /// @}

  constexpr const Teuchos::RCP<const SpaceT>& space() const { return( AbstractField<SpaceT>::space_ ); }

  constexpr const MPI_Comm& comm() const { return( space()->comm() ); }

  /// \name comunication methods.
	/// \brief highly dependent on underlying storage should only be used by
	/// Operator or on top field implementer.
  /// \{

  void changed( const int& dir ) const {
    exchangedState_[dir] = false;
  }
  void changed() const {
    for( int dir=0; dir<SpaceT::sdim; ++dir )
      changed( dir );
  }


  bool is_exchanged( const int& dir ) const {
    return( exchangedState_[dir] );
  }
  bool is_exchanged() const {
    bool all_exchanged = true;
    for( int dir=0; dir<SpaceT::sdim; ++dir )
      all_exchanged = all_exchanged && is_exchanged(dir);
    return( all_exchanged );
  }

  /// \brief updates ghost layers
  void exchange( const int& dir ) const {
		int ones[3] = {0,0,0};
    if( !exchangedState_[dir] ) {
      F_exchange(
          static_cast<int>( SpaceT::sdim ),
          MPI_Comm_c2f( space()->getProcGrid()->getCommWorld() ),
          space()->getProcGrid()->getRankL(),
          space()->getProcGrid()->getRankU(),
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->getBCLocal()->getBCL(),
          space()->getBCLocal()->getBCU(),
          space()->sInd(EField::S),
          space()->eInd(EField::S),
//				space()->sIndB(fType_), // should it work
//        space()->eIndB(fType_),
					ones,
          space()->nLoc(),
          1+dir,
          1+(int)fType_,
          s_);
      exchangedState_[dir] = true;
    }
  }
	void exchange() const {
		for( int dir=0; dir<SpaceT::sdim; ++dir )
			exchange( dir );
	}

  void setExchanged( const int& dir ) const {
    exchangedState_[dir] = true;
  }
  void setExchanged(  ) const {
    for( int dir=0; dir<SpaceT::sdim; ++dir )
      changed( dir );
  }

  /// \}

	/// \name indexing
	/// @{ 

protected:

	/// \brief stride in X direction
	constexpr Ordinal stride0() const {
		return( 1 );
	}

	/// \brief stride in Y direction
	constexpr Ordinal stride1() const {
		return( space()->nLoc(0)+space()->bu(0)-space()->bl(0)+1 );
	}

	/// \brief stride in Z direction
	constexpr Ordinal stride2() const {
		return(
				( space()->nLoc(0)+space()->bu(0)-space()->bl(0)+1 )*(
					space()->nLoc(1)+space()->bu(1)-space()->bl(1)+1 ) );
	}


	/// \brief stride
	///
	/// \param dir direction of stride
	constexpr Ordinal stride( const int& dir ) const {
		return(
				(0==dir)?
					stride0():
					( (1==dir)?
						stride1():
						stride2()
					)
				);
	}


	/// \brief computed index
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	inline constexpr Ordinal index( const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {
		return( (i-space()->bl(0)) +
				    (j-space()->bl(1))*stride1() +
				    (k-space()->bl(2))*stride2() );
	}

	/// \brief field access
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	///
	/// \return const reference
	inline constexpr const Scalar& at( const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {
		return( s_[ index(i,j,k) ] );
	}

	/// \brief field access
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	///
	/// \return reference
	inline Scalar& at( const Ordinal& i, const Ordinal& j, const Ordinal& k )  {
		return( s_[ index(i,j,k) ] );
	}

	/// \brief field access
	///
	/// \param i index coordinate 
	inline constexpr const Scalar& at( const Ordinal* const i ) const {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}
	/// \brief field access
	///
	/// \param i index coordinate 
	inline Scalar& at( const Ordinal* const i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}

	/// \brief field access
	///
	/// \param i index coordinate 
	inline Scalar& at( const Teuchos::Tuple<const Ordinal,3>& i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}
	/// \brief field access
	///
	/// \param i index coordinate 
	inline constexpr const Scalar& at( const Teuchos::Tuple<const Ordinal,3>& i ) const {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}

public:

	/// \brief field access
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	///
	/// \return const reference
	inline constexpr const Scalar& operator()( const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {
		return( s_[ index(i,j,k) ] );
	}

	/// \brief field access
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	///
	/// \return reference
	inline Scalar& operator()( const Ordinal& i, const Ordinal& j, const Ordinal& k )  {
		return( s_[ index(i,j,k) ] );
	}

	/// \brief field access
	///
	/// \param i index coordinate 
	inline constexpr const Scalar& operator()( const Ordinal* const i ) const {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}
	/// \brief field access
	///
	/// \param i index coordinate 
	inline Scalar& operator()( const Ordinal* const i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}

	/// \brief field access
	///
	/// \param i index coordinate 
	inline Scalar& operator()( const Teuchos::Tuple<const Ordinal,3>& i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}
	/// \brief field access
	///
	/// \param i index coordinate 
	inline constexpr const Scalar& operator()( const Teuchos::Tuple<const Ordinal,3>& i ) const {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}

}; // end of class ScalarField




/// \brief creates a scalar field(vector) belonging to a space
///
/// \param space scalar Vector Space to which returned vector belongs
/// \param fType
/// \return scalar vector
/// \relates ScalarField
template<class SpaceT>
Teuchos::RCP< ScalarField<SpaceT> >
createScalarField(
		const Teuchos::RCP<const SpaceT >& space,
		EField fType=EField::S ) {

	return( Teuchos::rcp(
			new ScalarField<SpaceT>( space, true, fType ) ) );
}

///  @} 

} // end of namespace Pimpact


#ifdef COMPILE_ETI
#include "Pimpact_Space.hpp"
extern template class Pimpact::ScalarField< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::ScalarField< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::ScalarField< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::ScalarField< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_SCALARFIELD_HPP
