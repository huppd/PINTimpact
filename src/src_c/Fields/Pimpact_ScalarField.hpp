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
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief important basic Vector class
/// vector for a scalar field, e.g.: pressure,
/// \ingroup Field
template<class SpaceType>
class ScalarField : private AbstractField< SpaceType > {

public:

  using SpaceT = SpaceType;

protected:

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  using ScalarArray = ST*;
  using State = Teuchos::Tuple<bool,3>;

	using SW = typename SpaceT::SW;

  ScalarArray s_;

  const bool owning_;

  State exchangedState_;

  const F fType_; /// < make template parameter (default:=S)

	void allocate() {
		OT n = getStorageSize();
		s_ = new ST[n];
	}

public:

  ScalarField( const Teuchos::RCP<const SpaceT>& space, bool owning=true, F fType=F::S ):
    AbstractField<SpaceT>( space ),
    owning_(owning),
    exchangedState_( Teuchos::tuple( true, true, true ) ),
    fType_(fType) {

    if( owning_ ) {
			allocate();
			init();
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
					init();
					break;
				case ECopy::Deep:
					*this = sF;
					break;
			}
		}
		};


	~ScalarField() { if( owning_ ) delete[] s_; }


  Teuchos::RCP<ScalarField> clone( ECopy copyType=ECopy::Deep ) const {

		Teuchos::RCP<ScalarField> mv = Teuchos::rcp( new ScalarField( space(), true, this->fType_ ) );

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
	constexpr OT getLength() {

		Teuchos::RCP<const BoundaryConditionsGlobal<SpaceT::dimension> > bc =
			space()->getBCGlobal();

		OT vl = 1;

		for( int dir = 0; dir<SpaceT::sdim; ++dir ) {
			vl *= space()->nGlo(dir) +
				( (BC::Periodic==bc->getBCL(dir))?
					-1:
					(fType_==dir)?1:0);
		}

		return( vl );
	}



  /// @}
  /// \name Update methods
  /// @{

  /// \brief Replace \c this with \f$\alpha a + \beta B\f$.
	/// \todo make checks for spaces and k
	void add( const ST& alpha, const ScalarField& a, const ST& beta, const
			ScalarField& b, const B& wb=B::Y ) {

		assert( a.getType()==b.getType() );
		assert( getType()==b.getType() );
#ifndef NDEBUG
		for( int dir=0; dir<3; ++dir ) {
			bool same_space = space()->nLoc(dir)>a.space()->nLoc(dir) || 
				space()->nLoc(dir)>b.space()->nLoc(dir);
			assert( !same_space );
			bool consistent_space = (
					(a.space()->nLoc(dir)-1)%(space()->nLoc(dir)-1) )!=0 || (
					(b.space()->nLoc(dir)-1)%(space()->nLoc(dir)-1) )!=0 ;
			assert( !consistent_space );
		}
#endif
		Teuchos::Tuple<OT,3> da;
		Teuchos::Tuple<OT,3> db;

		for( int dir=0; dir<3; ++dir ) {
			da[dir] = ( a.space()->nLoc(dir)-1 )/( space()->nLoc(dir)-1 );
			db[dir] = ( b.space()->nLoc(dir)-1 )/( space()->nLoc(dir)-1 );
		}

		for( OT k=space()->si(fType_,Z,wb); k<=space()->ei(fType_,Z,wb); ++k )
			for( OT j=space()->si(fType_,Y,wb); j<=space()->ei(fType_,Y,wb); ++j )
				for( OT i=space()->si(fType_,X,wb); i<=space()->ei(fType_,X,wb); ++i )
					at(i,j,k) = alpha*a.at( (i-1)*da[0]+1, (j-1)*da[1]+1,(k-1)*da[2]+1 )
						         + beta*b.at( (i-1)*db[0]+1, (j-1)*db[1]+1,(k-1)*db[2]+1 );

		changed();
	}


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
	void abs( const ScalarField& y, const B& bcYes=B::Y ) {

		for( int dir=0; dir<3; ++dir )
			assert( space()->nLoc(dir)==y.space()->nLoc(dir) );

		for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
			for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
				for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
					at(i,j,k) = std::fabs( y.at(i,j,k) );

		changed();
	}


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal( const ScalarField& y, const B& bcYes=B::Y ) {

#ifndef NDEBUG
		for( int dir=0; dir<3; ++dir ) {
			bool same_space = space()->nLoc(dir)!=y.space()->nLoc(dir);
			assert( !same_space );
		}
#endif

		for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
			for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
				for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
					at(i,j,k) = Teuchos::ScalarTraits<ST>::one()/ y.at(i,j,k);

    changed();
  }


  /// \brief Scale each element of the vector with \c alpha.
	void scale( const ST& alpha, const B& wB=B::Y ) {

		for( OT k=space()->si(fType_,Z,wB); k<=space()->ei(fType_,Z,wB); ++k )
			for( OT j=space()->si(fType_,Y,wB); j<=space()->ei(fType_,Y,wB); ++j )
				for( OT i=space()->si(fType_,X,wB); i<=space()->ei(fType_,X,wB); ++i )
					at(i,j,k) *= alpha;

		changed();
	}


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot y_i \quad \mbox{for } i=1,\dots,n \f]
	void scale( const ScalarField& y, const B& bcYes=B::Y ) {

#ifndef NDEBUG
		for( int dir=0; dir<3; ++dir ) {
			bool same_space = space()->nLoc(dir)!=y.space()->nLoc(dir);
			assert( !same_space );
		}
#endif

		for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
			for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
				for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
					at(i,j,k) *= y.at(i,j,k);
		changed();
	}

  /// @}
  /// \name Norm method(reductions)
  /// @{

	/// \brief Compute a local scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
	constexpr ST dotLoc( const ScalarField& y, const B& bcYes=B::Y ) {

#ifndef NDEBUG
		for( int dir=0; dir<3; ++dir ) {
			bool same_space = space()->nLoc(dir)!=y.space()->nLoc(dir);
			assert( !same_space );
		}
#endif

		ST b = Teuchos::ScalarTraits<ST>::zero();

		for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
			for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
				for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
					b += at(i,j,k)*y.at(i,j,k);

		return( b );
	}

	/// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y
	/// and \c this, i.e.\f$b = y^H this\f$.
	constexpr ST dot( const ScalarField& y, const B& bcYes=B::Y ) {
		return( this->reduce( comm(), dotLoc( y, bcYes ) ) );
	}

  constexpr ST normLoc1( const B& bcYes=B::Y ) {

    ST normvec = Teuchos::ScalarTraits<ST>::zero();

		for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
			for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
				for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
					normvec += std::fabs( at(i,j,k) );

    return( normvec );
  }


  constexpr ST normLoc2( const B& bcYes=B::Y ) {

    ST normvec = Teuchos::ScalarTraits<ST>::zero();

		for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
			for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
				for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
					normvec += std::pow( at(i,j,k), 2 );

    return( normvec );
  }


  constexpr ST normLocInf( const B& bcYes=B::Y ) {

    ST normvec = Teuchos::ScalarTraits<ST>::zero();

		for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
			for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
				for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
					normvec = std::fmax( std::fabs(at(i,j,k)), normvec );

    return( normvec );
  }

	constexpr ST normLoc( Belos::NormType type = Belos::TwoNorm, const B& bcYes=B::Y ) {

		return(
				( Belos::OneNorm==type)?
					normLoc1(bcYes):
					(Belos::TwoNorm==type)?
						normLoc2(bcYes):
						normLocInf(bcYes) );
	}


  /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  constexpr ST norm( Belos::NormType type = Belos::TwoNorm, const B& bcYes=B::Y ) {

		ST normvec = this->reduce(
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
  constexpr ST normLoc( const ScalarField& weights, const B& bcYes=B::Y ) {

		for( int dir=0; dir<3; ++dir )
			assert( space()->nLoc(dir)==weights.space()->nLoc(dir) );

    ST normvec = Teuchos::ScalarTraits<ST>::zero();

		for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
			for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
				for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
					normvec += at(i,j,k)*at(i,j,k)*weights.at(i,j,k)*weights.at(i,j,k);

    return( normvec );
  }


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST norm( const ScalarField& weights, const B& bcYes=B::Y ) {
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
  void random( bool useSeed = false, const B& bcYes=B::Y , int seed = 1 ) {

		std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis( -0.5, 0.5 );

		for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
			for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
				for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
					at(i,j,k) = dis(gen);

		if( !space()->getProcGrid()->participating() )
			for( OT k=space()->si(fType_,Z,B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
				for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j )
					for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i )
						at(i,j,k) = Teuchos::ScalarTraits<ST>::zero();
		changed();
  }


  /// \brief Replace each element of the vector  with \c alpha.
	/// \param alpha init value
	/// \param bcYes also initializing the boundary values
	void init( const ST& alpha = Teuchos::ScalarTraits<ST>::zero(), const B& bcYes=B::Y ) {

		if( B::Y==bcYes ){
			std::fill_n( s_, getStorageSize(), alpha );
			exchangedState_[X] = true;
			exchangedState_[Y] = true;
			exchangedState_[Z] = true;
		}
		else{
			for( OT k=space()->si(fType_,Z,bcYes); k<=space()->ei(fType_,Z,bcYes); ++k )
				for( OT j=space()->si(fType_,Y,bcYes); j<=space()->ei(fType_,Y,bcYes); ++j )
					for( OT i=space()->si(fType_,X,bcYes); i<=space()->ei(fType_,X,bcYes); ++i )
						at(i,j,k) = alpha;

			if( !space()->getProcGrid()->participating() )
				for( OT k=space()->si(fType_,Z,B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
					for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j )
						for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i )
							at(i,j,k) = Teuchos::ScalarTraits<ST>::zero();
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
	void initFromFunction( Functor&& func, const Add& add=Add::N ) {

		Teuchos::RCP<const CoordinatesLocal<ST,OT,SpaceT::dimension,SpaceT::dimNC> > coord =
			space()->getCoordinatesLocal();
		Teuchos::RCP<const DomainSize<ST,SpaceT::sdim> > domain = space()->getDomainSize();

		const B& bY = B::Y;

		for( OT k=space()->si(fType_,Z,bY); k<=space()->ei(fType_,Z,bY); ++k )
			for( OT j=space()->si(fType_,Y,bY); j<=space()->ei(fType_,Y,bY); ++j )
				for( OT i=space()->si(fType_,X,bY); i<=space()->ei(fType_,X,bY); ++i ) {
					if( Add::Y==add )
						at(i,j,k) += func(
								( coord->getX(fType_,X,i)-domain->getOrigin(X) )/domain->getSize(X),
								( coord->getX(fType_,Y,j)-domain->getOrigin(Y) )/domain->getSize(Y),
								( coord->getX(fType_,Z,k)-domain->getOrigin(Z) )/domain->getSize(Z) );
					else
						at(i,j,k) = func(
								( coord->getX(fType_,X,i)-domain->getOrigin(X) )/domain->getSize(X),
								( coord->getX(fType_,Y,j)-domain->getOrigin(Y) )/domain->getSize(Y),
								( coord->getX(fType_,Z,k)-domain->getOrigin(Z) )/domain->getSize(Z) );
				}
		changed();
	}


	///  \brief initializes including boundaries to zero 
	void initField( Teuchos::ParameterList& para, const Add& add=Add::N ) {

		EScalarField type =
			string2enum( para.get<std::string>( "Type", "constant" ) );

		switch( type ) {
			case ConstField :
				{
					if( Add::N==add ) init();
					break;
				}
			case Grad2D_inX :
				{
					ST a = para.get<ST>( "dx", Teuchos::ScalarTraits<ST>::one() );
					initFromFunction(
							[&a] (ST x, ST y, ST z)->ST { return( a*(x-0.5) ); },
							add );
					break;
				}
			case Grad2D_inY :
				{
					ST a = para.get<ST>( "dy", Teuchos::ScalarTraits<ST>::one() );
					initFromFunction(
							[&a] (ST x, ST y, ST z)->ST { return( a*(y-0.5) ); },
							add );
					break;
				}
			case Grad2D_inZ :
				{
					ST a = para.get<ST>( "dz", Teuchos::ScalarTraits<ST>::one() );
					initFromFunction(
							[&a] (ST x, ST y, ST z)->ST { return( a*(z-0.5) ); },
							add );
					break;
				}
			case Poiseuille2D_inX :
				{
					initFromFunction(
							[] (ST x, ST y, ST z)->ST { return( 4.*x*(1.-x) ); },
							add );
					break;
				}
			case Poiseuille2D_inY :
				{
					initFromFunction(
							[] (ST x, ST y, ST z)->ST { return( 4.*y*(1.-y) ); },
							add );
					break;
				}
			case Poiseuille2D_inZ :
				{
					initFromFunction(
							[] (ST x, ST y, ST z)->ST { return( 4.*z*(1.-z) ); },
							add );
					break;
				}
			case FPoint :
				{
					ST xc[3] = { 
						para.get<ST>( "c_x", Teuchos::ScalarTraits<ST>::one() ),
						para.get<ST>( "c_y", space()->getDomainSize()->getSize( Y )/2. ),
						para.get<ST>( "c_z", space()->getDomainSize()->getSize( Z )/2. ) };
					ST amp = para.get<ST>( "amp", Teuchos::ScalarTraits<ST>::one() );
					ST sig[3] = {
						para.get<ST>( "sig_x", 0.2 ),
						para.get<ST>( "sig_y", 0.2 ),
						para.get<ST>( "sig_z", 0.2 ) };
					initFromFunction(
							[&xc,&amp,&sig] (ST x, ST y, ST z)->ST {
								return( amp*std::exp(
										-std::pow( (x-xc[0])/sig[0], 2 )
										-std::pow( (x-xc[1])/sig[1], 2 )
										-std::pow( (x-xc[2])/sig[2], 2 ) ) ); },
							add );
					break;
				}
		}

		if( !space()->getProcGrid()->participating() ) // not sure why?
			init( Teuchos::ScalarTraits<ST>::zero() );

		changed();
	}


	/// \brief initializes VectorField with the initial field defined in Fortran
	/// \deprecated
	void initField( EScalarField fieldType, ST alpha=Teuchos::ScalarTraits<ST>::zero() ) {

		switch( fieldType ) {
			case ConstField :
				{
					init( alpha );
					break;
				}
			case Grad2D_inX :
				{
					ST a = (std::fabs(alpha)<Teuchos::ScalarTraits<ST>::eps())?Teuchos::ScalarTraits<ST>::one():alpha;
					initFromFunction( [&a] (ST x, ST y, ST z)->ST { return( a*(x-0.5) ); } );
				break;
				}
			case Grad2D_inY :
				{
					ST a = (std::fabs(alpha)<Teuchos::ScalarTraits<ST>::eps())?Teuchos::ScalarTraits<ST>::one():alpha;
					initFromFunction( [&a] (ST x, ST y, ST z)->ST { return( a*(y-0.5) ); } );
					break;
				}
			case Grad2D_inZ :
				{
					ST a = (std::fabs(alpha)<Teuchos::ScalarTraits<ST>::eps())?Teuchos::ScalarTraits<ST>::one():alpha;
					initFromFunction( [&a] (ST x, ST y, ST z)->ST { return( a*(z-0.5) ); } );
					break;
				}
			case Poiseuille2D_inX :
				{
					initFromFunction( [] (ST x, ST y, ST z)->ST { return( 4.*x*(1.-x) ); } );
					break;
				}
			case Poiseuille2D_inY :
				{
					initFromFunction( [] (ST x, ST y, ST z)->ST { return( 4.*y*(1.-y) ); } );
					break;
				}
			case Poiseuille2D_inZ :
				{
					initFromFunction( [] (ST x, ST y, ST z)->ST { return( 4.*z*(1.-z) ); } );
					break;
				}
			case FPoint :
				{
					ST xc[3] = { 0.5, 0.5, 0.5 };
					ST amp = alpha; 
					ST sig[3] = { 0.2, 0.2, 0.2 };
					initFromFunction( [&xc,&amp,&sig] (ST x, ST y, ST z)->ST {
							return( amp*std::exp(
									-std::pow( (x-xc[0])/sig[0], 2 )
									-std::pow( (x-xc[1])/sig[1], 2 )
									-std::pow( (x-xc[2])/sig[2], 2 ) ) ); }
							);
					break;
				}
		}

		if( !space()->getProcGrid()->participating() )
			init( Teuchos::ScalarTraits<ST>::zero() );
		changed();
	}


	/// \brief for Dirichlet BC extrapolate the velocity points outside the domain such that the
	///  interpolated value on the boundary is zero
	///
	/// \test Neumann BC
	/// \param trans transposed
	void extrapolateBC( const Belos::ETrans& trans=Belos::NOTRANS ) {

		switch( trans ) {
			case( Belos::NOTRANS ): {
				switch( fType_ ) {
					case( F::U ):  {
						using StencD = Stencil< ST, OT, 0, SW::DL(0), SW::DU(0) >;

						StencD c_( space()->nLoc(X) );

						if( space()->bcl(X)==BC::Neumann || space()->bcu(X)==BC::Neumann )
							FD_getDiffCoeff(
									1,
									space()->nLoc(X),
									space()->bl(X),
									space()->bu(X),
									space()->dl(X),
									space()->du(X),
									space()->getBCLocal()->getBCL(X),
									space()->getBCLocal()->getBCU(X),
									space()->getShift(X),
									3,
									1,
									1,
									0,
									false, // mapping
									space()->getStencilWidths()->getDimNcbD(X),
									space()->getStencilWidths()->getNcbD(X),
									space()->getCoordinatesLocal()->getX( F::U, X ),
									space()->getCoordinatesLocal()->getX( F::S, X ),
									c_.get() );

						if( 0 < space()->bcl(X) ) {
							OT i = space()->si(fType_,X,B::Y);
							for( OT k=space()->si(fType_,Z, B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
								for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j ) {
									at(i,j,k) = 0.;
									if( BC::Dirichlet==space()->bcl(X) ) {
										for( OT ii=0; ii<=SW::DU(X); ++ii )
											at(i,j,k) -= at(1+ii,j,k)*space()->getInterpolateV2S()->getC(X,1,ii)/space()->getInterpolateV2S()->getC(X,1,-1);
									}
									else if( BC::Neumann==space()->bcl(X) ) {
										for( OT ii=0; ii<=SW::DU(X); ++ii )
											at(i,j,k) -= at(1+ii,j,k)*c_(1,ii)/c_(1,-1);
									}
								}
						}
						if( 0 < space()->bcu(X) ) {

							OT i = space()->ei(fType_,X,B::Y);
							for( OT k=space()->si(fType_,Z, B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
								for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j ) {
									at(i,j,k) = 0.;
									if( BC::Dirichlet==space()->bcu(X) ) {
										for( OT ii=SW::DL(X); ii<=-1; ++ii )
											at(i,j,k) -= space()->getInterpolateV2S()->getC(X,i,ii)*at(i+ii,j,k)/space()->getInterpolateV2S()->getC(X,i,0);
									}
									else if( BC::Neumann==space()->bcu(X) ) {
										for( OT ii=SW::DL(X); ii<=-1; ++ii )
											at(i,j,k) -= c_(i,ii)*at(i+ii,j,k)/c_(i,0);
									}
								}
						}
						break;
					}
					case( F::V ) : {
						using StencD = Stencil< ST, OT, 0, SW::DL(0), SW::DU(0) >;

						StencD c_( space()->nLoc(Y) );

						if( space()->bcl(Y)==BC::Neumann || space()->bcu(Y)==BC::Neumann )
							FD_getDiffCoeff(
									1,
									space()->nLoc(Y),
									space()->bl(Y),
									space()->bu(Y),
									space()->dl(Y),
									space()->du(Y),
									space()->getBCLocal()->getBCL(Y),
									space()->getBCLocal()->getBCU(Y),
									space()->getShift(Y),
									3,
									1,
									1,
									0,
									false, // mapping
									space()->getStencilWidths()->getDimNcbD(Y),
									space()->getStencilWidths()->getNcbD(Y),
									space()->getCoordinatesLocal()->getX( F::V, Y ),
									space()->getCoordinatesLocal()->getX( F::S, Y ),
									c_.get() );

						if( 0 < space()->bcl(Y) ) {
							OT j = space()->si(fType_,Y,B::Y);
							for( OT k=space()->si(fType_,Z, B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
								for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i ) {
									at(i,j,k) = 0.;
									if( BC::Dirichlet==space()->bcl(Y) ) {
										for( OT jj=0; jj<=SW::DU(Y); ++jj )
											at(i,j,k) -= at(i,1+jj,k)*space()->getInterpolateV2S()->getC(Y,1,jj)/space()->getInterpolateV2S()->getC(Y,1,-1);  
									}
									else if( BC::Neumann==space()->bcl(Y) ) {
										for( OT jj=0; jj<=SW::DU(Y); ++jj )
											at(i,j,k) -= at(i,1+jj,k)*c_(1,jj)/c_(1,-1);  
									}
								}
						}
						if( 0 < space()->bcu(Y) ) {
							OT j = space()->ei(fType_,Y,B::Y);
							for( OT k=space()->si(fType_,Z, B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
								for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i ) {
									at(i,j,k) = 0.;
									if( BC::Dirichlet==space()->bcu(Y) ) {
										for( OT jj=SW::DL(Y); jj<=-1; ++jj )
											at(i,j,k) -= space()->getInterpolateV2S()->getC(Y,j,jj)*at(i,j+jj,k)/space()->getInterpolateV2S()->getC(Y,j,0);
									}
									else if( BC::Neumann==space()->bcu(Y) ) {
										for( OT jj=SW::DL(Y); jj<=-1; ++jj )
											at(i,j,k) -= c_(j,jj)*at(i,j+jj,k)/c_(j,0);
									}
								}
						}
						break;
					}
					case( F::W ) : {
						using StencD = Stencil< ST, OT, 0, SW::DL(0), SW::DU(0) >;

						StencD c_( space()->nLoc(Z) );

						if( space()->bcl(Z)==BC::Neumann || space()->bcu(Z)==BC::Neumann )
							FD_getDiffCoeff(
									1,
									space()->nLoc(Z),
									space()->bl(Z),
									space()->bu(Z),
									space()->dl(Z),
									space()->du(Z),
									space()->getBCLocal()->getBCL(Z),
									space()->getBCLocal()->getBCU(Z),
									space()->getShift(Z),
									3,
									1,
									1,
									0,
									false, // mapping
									space()->getStencilWidths()->getDimNcbD(Z),
									space()->getStencilWidths()->getNcbD(Z),
									space()->getCoordinatesLocal()->getX( F::W, Z ),
									space()->getCoordinatesLocal()->getX( F::S, Z ),
									c_.get() );

						if( space()->bcl(Z) > 0 ) {
							OT k = space()->si(fType_,Z,B::Y);
							for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j )
								for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i ) {
									at(i,j,k) = 0.;
									if( BC::Dirichlet==space()->bcl(Z) ) {
										for( OT kk=0; kk<=SW::DU(Z); ++kk )
											at(i,j,k) -= space()->getInterpolateV2S()->getC(Z,1,kk)*at(i,j,1+kk)/space()->getInterpolateV2S()->getC(Z,1,-1);  
									}
									else if( BC::Neumann==space()->bcl(Z) ) {
										for( OT kk=0; kk<=SW::DU(Z); ++kk )
											at(i,j,k) -= c_(1,kk)*at(i,j,1+kk)/c_(1,-1);  
									}
								}
						}
						if( space()->bcu(Z) > 0 ) {
							OT k = space()->ei(fType_,Z,B::Y);
							for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j )
								for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i ) {
									at(i,j,k) = 0.;
									if( BC::Dirichlet==space()->bcu(Z) ) {
										for( OT kk=SW::DL(Z); kk<=-1; ++kk )
											at(i,j,k) -= space()->getInterpolateV2S()->getC(Z,k,kk)*at(i,j,k+kk)/space()->getInterpolateV2S()->getC(Z,k,0);
									}
									else if( BC::Neumann==space()->bcu(Z) ) {
										for( OT kk=SW::DL(Z); kk<=-1; ++kk )
											at(i,j,k) -= c_(k,kk)*at(i,j,k+kk)/c_(k,0);
									}
								}
						}
						break;
					}
					case( F::S ) : break;
					//case( F::end ) : break;
				}
				break;
			}
			case Belos::TRANS : {

				switch( fType_ ) {
					case( F::U ) : {
						if( space()->bcl(X) > 0 ) {
							OT i = space()->si(fType_,X,B::Y);
							for( OT k=space()->si(fType_,Z, B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
								for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j ) {
									for( OT ii=0; ii<=SW::DU(X); ++ii )
										at(i+ii+1,j,k) -= at(i,j,k)*space()->getInterpolateV2S()->getC(X,1,ii)/space()->getInterpolateV2S()->getC(X,1,-1);  
									at(i,j,k) = 0.;
								}
						}
						if( space()->bcu(X) > 0 ) {
							OT i = space()->ei(fType_,X,B::Y);
							for( OT k=space()->si(fType_,Z, B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
								for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j ) {
									for( OT ii=SW::DL(X); ii<=-1; ++ii )
										at(i+ii,j,k) -= space()->getInterpolateV2S()->getC(X,i,ii)*at(i,j,k)/space()->getInterpolateV2S()->getC(X,i,0); 
									at(i,j,k) = 0.;
								}
						}
						break;
					}
					case( F::V ) : {
						if( space()->bcl(Y) > 0 ) {
							OT j = space()->si(fType_,Y,B::Y);
							for( OT k=space()->si(fType_,Z, B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
								for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i ) {
									for( OT jj=0; jj<=SW::DU(Y); ++jj )
										at(i,j+jj+1,k) -= at(i,j,k)*space()->getInterpolateV2S()->getC(Y,1,jj)/space()->getInterpolateV2S()->getC(Y,1,-1);  
									at(i,j,k) = 0.;
								}
						}
						if( space()->bcu(Y) > 0 ) {
							OT j = space()->ei(fType_,Y,B::Y);
							for( OT k=space()->si(fType_,Z, B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
								for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i ) {
									for( OT jj=SW::DL(Y); jj<=-1; ++jj )
										at(i,j+jj,k) -= space()->getInterpolateV2S()->getC(Y,j,jj)*at(i,j,k)/space()->getInterpolateV2S()->getC(Y,j,0); 
									at(i,j,k) = 0.;
								}
						}
						break;
					}
					case( F::W ) : {
						if( space()->bcl(Z) > 0 ) {
							OT k = space()->si(fType_,Z,B::Y);
							for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j )
								for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i ) {
									at(i,j,k) /= space()->getInterpolateV2S()->getC(Z,1,-1);
									for( OT kk=0; kk<=SW::DU(Z); ++kk )
										at(i,j,k+kk+1) -= space()->getInterpolateV2S()->getC(Z,1,kk)*at(i,j,k);  
									at(i,j,k) = 0.;
								}
						}
						if( space()->bcu(Z) > 0 ) {
							OT k = space()->ei(fType_,Z,B::Y);
							for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j )
								for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i ) {
									at(i,j,k) /= space()->getInterpolateV2S()->getC(Z,k,0);
									for( OT kk=SW::DL(Z); kk<=-1; ++kk )
										at(i,j,k+kk) -= space()->getInterpolateV2S()->getC(Z,k,kk)*at(i,j,k); 
									at(i,j,k) = 0.;
								}
						}
						break;
					}
					case( F::S ) : break;
				}
			}
			case Belos::CONJTRANS : break;
		}
	}


	/// \brief levels field if scalar field
	void level() const {

		if( F::S == fType_ ) {
			
			ST pre0 = Teuchos::ScalarTraits<ST>::zero();

			for( OT k=space()->si(fType_,Z); k<=space()->ei(fType_,Z); ++k )
				for( OT j=space()->si(fType_,Y); j<=space()->ei(fType_,Y); ++j )
					for( OT i=space()->si(fType_,X); i<=space()->ei(fType_,X); ++i )
						pre0 += at(i,j,k);

			pre0 = this->reduce( space()->comm(), pre0 );
			pre0 /= static_cast<ST>( getLength() );

			for( OT k=space()->si(fType_,Z); k<=space()->ei(fType_,Z); ++k )
				for( OT j=space()->si(fType_,Y); j<=space()->ei(fType_,Y); ++j )
					for( OT i=space()->si(fType_,X); i<=space()->ei(fType_,X); ++i )
						const_cast<ScalarField*>(this)->at(i,j,k) -= pre0;
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

		Teuchos::Tuple<OT,3> cw;
		for(int i=0; i<3; ++i)
			cw[i] = space()->nLoc(i) + SW::BU(i) - SW::BL(i) + 1;

		for( OT k=space()->si(fType_,Z,B::Y); k<=space()->ei(fType_,Z,B::Y); ++k )
			for( OT j=space()->si(fType_,Y,B::Y); j<=space()->ei(fType_,Y,B::Y); ++j )
				for( OT i=space()->si(fType_,X,B::Y); i<=space()->ei(fType_,X,B::Y); ++i )
					out << i << "\t" << j << "\t" << k << "\t" << at(i,j,k) << "\n";
  }


  /// Write the ScalarField to an hdf5 file, the velocities are interpolated to the pressure points
  /// \todo add restart
	void write( int count=0 , bool restart=false ) const {

		if( 0==space()->rankS() )
			switch(fType_) {
				case F::U:
					std::cout << "writing velocity field x(" << count << ") ...\n";
					break;
				case F::V:
					std::cout << "writing velocity field y(" << count << ") ...\n";
					break;
				case F::W:
					std::cout << "writing velocity field z(" << count << ") ...\n";
					break;
				case F::S:
					std::cout << "writing pressure field  (" << count << ") ...\n";
					Teuchos::Tuple<OT,3> N;
					for( int i=0; i<3; ++i ) {
						N[i] = space()->nGlo(i);
						if( space()->getBCGlobal()->getBCL(i)==Pimpact::BC::Periodic )
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

			if( F::S != fType_ ) {
				temp = Teuchos::rcp(
						new ScalarField<SpaceT>( space(), true, F::S ) );
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
            space()->sInd(F::S),
            space()->eInd(F::S),
            space()->getStencilWidths()->getLS(),
            space()->np(),
            space()->ib(),
            space()->getShift(),
            (int)fType_,
            count,
            (F::S==fType_)?9:10,
						(F::S==fType_)?s_:temp->s_,
						space()->getCoordinatesGlobal()->getX(F::S,0),
						space()->getCoordinatesGlobal()->getX(F::S,1),
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
            space()->sInd(F::S),
            space()->eInd(F::S),
            space()->getStencilWidths()->getLS(),
            space()->np(),
            space()->ib(),
            space()->getShift(),
						(int)fType_+1,
						(int)F::S+1,
            count,
            (F::S==fType_)?9:10,
            stride,
            (F::S==fType_)?s_:temp->s_,
            space()->getCoordinatesGlobal()->getX(F::S,0),
            space()->getCoordinatesGlobal()->getX(F::S,1),
            space()->getCoordinatesGlobal()->getX(F::S,2),
            space()->getCoordinatesGlobal()->getX(F::U,0),
            space()->getCoordinatesGlobal()->getX(F::V,1),
            space()->getCoordinatesGlobal()->getX(F::W,2),
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
          (F::S==fType_)?9:10,
          stride,
          s_,
          space()->getCoordinatesGlobal()->getX(F::S,0),
          space()->getCoordinatesGlobal()->getX(F::S,1),
          space()->getCoordinatesGlobal()->getX(F::S,2),
          space()->getCoordinatesGlobal()->getX(F::U,0),
          space()->getCoordinatesGlobal()->getX(F::V,1),
          space()->getCoordinatesGlobal()->getX(F::W,2),
          space()->getDomainSize()->getRe(),
          space()->getDomainSize()->getAlpha2() );
    }
  }

	void read( int count=0 ) const {
	}


public:

  constexpr const F& getType() { return( fType_ ); }

   /// \name storage methods.
	 /// \brief highly dependent on underlying storage should only be used by
	 /// Operator or on top field implementer.  
   /// @{

  OT getStorageSize() const {

    OT n = 1;
    for(int i=0; i<3; ++i)
      n *= space()->nLoc(i)+SW::BU(i)-SW::BL(i)+1; // seems wrong: there a one was added for AMG, but it is not neede error seem to be in Impact there it should be (B1L+1:N1+B1U) probably has to be changed aganin for 3D

    return( n );
  }

	void setStoragePtr( ST*  array ) { s_ = array; }

	constexpr ScalarArray getRawPtr() { return( s_ ); }

	constexpr const ST* getConstRawPtr() { return( s_ ); }


  /// @}

  constexpr const Teuchos::RCP<const SpaceT>& space() { return( AbstractField<SpaceT>::space_ ); }

  constexpr const MPI_Comm& comm() { return( space()->comm() ); }

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
          space()->sInd(F::S),
          space()->eInd(F::S),
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
	constexpr OT stride0() {
		return( 1 );
	}

	/// \brief stride in Y direction
	constexpr OT stride1() {
		return( space()->nLoc(0)+SW::BU(0)-SW::BL(0)+1 );
	}

	/// \brief stride in Z direction
	constexpr OT stride2() {
		return(
				( space()->nLoc(0)+SW::BU(0)-SW::BL(0)+1 )*(
					space()->nLoc(1)+SW::BU(1)-SW::BL(1)+1 ) );
	}


	/// \brief stride
	///
	/// \param dir direction of stride
	constexpr OT stride( const int& dir ) {
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
	constexpr OT index( const OT& i, const OT& j, const OT& k ) {
		return( (i-SW::BL(0)) +
				    (j-SW::BL(1))*stride1() +
				    (k-SW::BL(2))*stride2() );
	}

	/// \brief field access
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	///
	/// \return const reference
	constexpr const ST& at( const OT& i, const OT& j, const OT& k ) {
		return( s_[ index(i,j,k) ] );
	}

	/// \brief field access
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	///
	/// \return reference
	ST& at( const OT& i, const OT& j, const OT& k )  {
		return( s_[ index(i,j,k) ] );
	}

	/// \brief field access
	///
	/// \param i index coordinate 
	constexpr const ST& at( const OT* const i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}
	/// \brief field access
	///
	/// \param i index coordinate 
	ST& at( const OT* const i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}

	/// \brief field access
	///
	/// \param i index coordinate 
	ST& at( const Teuchos::Tuple<const OT,3>& i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}
	/// \brief field access
	///
	/// \param i index coordinate 
	constexpr const ST& at( const Teuchos::Tuple<const OT,3>& i ) {
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
	constexpr const ST& operator()( const OT& i, const OT& j, const OT& k ) {
		return( s_[ index(i,j,k) ] );
	}

	/// \brief field access
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	///
	/// \return reference
	ST& operator()( const OT& i, const OT& j, const OT& k )  {
		return( s_[ index(i,j,k) ] );
	}

	/// \brief field access
	///
	/// \param i index coordinate 
	constexpr const ST& operator()( const OT* const i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}
	/// \brief field access
	///
	/// \param i index coordinate 
	ST& operator()( const OT* const i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}

	/// \brief field access
	///
	/// \param i index coordinate 
	ST& operator()( const Teuchos::Tuple<const OT,3>& i ) {
		return( s_[ index(i[0],i[1],i[2]) ] );
	}
	/// \brief field access
	///
	/// \param i index coordinate 
	constexpr const ST& operator()( const Teuchos::Tuple<const OT,3>& i ) {
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
		F fType=F::S ) {

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
