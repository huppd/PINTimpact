#pragma once
#ifndef PIMPACT_SCALARFIELD_HPP
#define PIMPACT_SCALARFIELD_HPP


#include <vector>
#include <iostream>

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
/// \note all indexing is done in Fortran
/// \ingroup Field
/// \todo think about using Teuchos::ArrayRCP instead of Scalar* should make delete
template<class SpaceType>
class ScalarField : private AbstractField< SpaceType > {

public:

  using SpaceT = SpaceType;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  static const int dimension = SpaceT::dimension;

protected:

  using ScalarArray = Scalar*;
  using MV = ScalarField< SpaceT >;
  using State = Teuchos::Tuple<bool,3>;

  ScalarArray s_;

  bool owning_;

  State exchangedState_;

  EField fType_;

private:

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
  /// \param copyType by default a DeepCopy is done but also allows to ShallowCopy
  ScalarField( const ScalarField& sF, ECopyType copyType=DeepCopy ):
    AbstractField<SpaceT>( sF.space() ),
    owning_( sF.owning_ ),
    exchangedState_( sF.exchangedState_ ),
    fType_( sF.fType_ ) {

    if( owning_ ) {

			allocate();

      switch( copyType ) {
      case ShallowCopy:
				initField();
        break;
      case DeepCopy:
        for( int i=0; i<getStorageSize(); ++i)
          s_[i] = sF.s_[i];
        break;
      }
    }

  };


	~ScalarField() { if( owning_ ) delete[] s_; }


  Teuchos::RCP<MV> clone( ECopyType copyType=DeepCopy ) const {

		Teuchos::RCP<MV> mv = Teuchos::rcp( new MV( space(), true, this->fType_ ) );

		switch( copyType ) {
			case ShallowCopy:
				break;
			case DeepCopy:
				for( int i=0; i<getStorageSize(); ++i)
					mv->getRawPtr()[i] = s_[i];
				break;
		}
    return( mv );

  }

  /// \name Attribute methods
  /// \{

  /// \brief returns the length of Field.
	Ordinal getLength( bool dummy=false ) const {

		Teuchos::RCP<const BoundaryConditionsGlobal<dimension> > bc =
			space()->getBCGlobal();

		Ordinal vl = 1;

		switch( fType_ ) {
			case EField::S : {

				Teuchos::Tuple<Ordinal,3> n;

				for(int i = 0; i<space()->dim(); ++i) {
					if( PeriodicBC==bc->getBCL(i) )
						n[i] = space()->nGlo(i)-1;
					else
						n[i] = space()->nGlo(i);
					vl *= n[i];
				}

				const int* bcl =  space()->getBCGlobal()->getBCL();
				const int* bcu =  space()->getBCGlobal()->getBCU();

				if( 2==space()->dim() ) {
					if( bcl[0]>0 && bcl[1]>0 ) vl -= 1;
					if( bcl[0]>0 && bcu[1]>0 ) vl -= 1;
					if( bcu[0]>0 && bcl[1]>0 ) vl -= 1;
					if( bcu[0]>0 && bcu[1]>0 ) vl -= 1;
				}
				else{
					if( bcl[0]>0 && bcl[1]>0 ) vl -= n[2] - ((bcl[2]==-1)?0:2);
					if( bcl[0]>0 && bcu[1]>0 ) vl -= n[2] - ((bcl[2]==-1)?0:2);
					if( bcu[0]>0 && bcl[1]>0 ) vl -= n[2] - ((bcl[2]==-1)?0:2);
					if( bcu[0]>0 && bcu[1]>0 ) vl -= n[2] - ((bcl[2]==-1)?0:2);

					if( bcl[0]>0 && bcl[2]>0 ) vl -= n[1] - ((bcl[1]==-1)?0:2);
					if( bcl[0]>0 && bcu[2]>0 ) vl -= n[1] - ((bcl[1]==-1)?0:2);
					if( bcu[0]>0 && bcl[2]>0 ) vl -= n[1] - ((bcl[1]==-1)?0:2);
					if( bcu[0]>0 && bcu[2]>0 ) vl -= n[1] - ((bcl[1]==-1)?0:2);

					if( bcl[1]>0 && bcl[2]>0 ) vl -= n[0] - ((bcl[0]==-1)?0:2);
					if( bcl[1]>0 && bcu[2]>0 ) vl -= n[0] - ((bcl[0]==-1)?0:2);
					if( bcu[1]>0 && bcl[2]>0 ) vl -= n[0] - ((bcl[0]==-1)?0:2);
					if( bcu[1]>0 && bcu[2]>0 ) vl -= n[0] - ((bcl[0]==-1)?0:2);

					if( bcl[0]>0 && bcl[1]>0 && bcl[2]>0 ) vl -= 1;
					if( bcu[0]>0 && bcl[1]>0 && bcl[2]>0 ) vl -= 1;
					if( bcl[0]>0 && bcu[1]>0 && bcl[2]>0 ) vl -= 1;
					if( bcl[0]>0 && bcl[1]>0 && bcu[2]>0 ) vl -= 1;

					if( bcu[0]>0 && bcu[1]>0 && bcu[2]>0 ) vl -= 1;
					if( bcl[0]>0 && bcu[1]>0 && bcu[2]>0 ) vl -= 1;
					if( bcu[0]>0 && bcl[1]>0 && bcu[2]>0 ) vl -= 1;
					if( bcu[0]>0 && bcu[1]>0 && bcl[2]>0 ) vl -= 1;
				}
				break;
			}
			default: {
				for( int j=0; j<space()->dim(); ++j) {
					if( fType_==j ) {
						vl *= space()->nGlo(j)-1;
					}
					else {
						if( PeriodicBC==bc->getBCL(j) )
							vl *= space()->nGlo(j)-2+1;
						else
							vl *= space()->nGlo(j)-2;
					}
				}
				break;
			}
		}
		return( vl );
	}


  /// \brief get number of stored Field's
  int getNumberVecs() const { return( 1 ); }


  /// @}
  /// \name Update methods
  /// @{

  /// \brief Replace \c this with \f$\alpha A + \beta B\f$.
  void add( const Scalar& alpha, const MV& A, const Scalar& beta, const MV& B ) {
    // add test for consistent VectorSpaces in debug mode
    if( s_==A.s_ && s_==B.s_ )
      scale( alpha+beta );
    else if( s_==A.s_ && s_!=B.s_ )
      SF_add2(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType_),
          space()->eInd(fType_),
          s_, B.s_,
          alpha, beta );
    else if( s_!=A.s_ && s_==B.s_ )
      SF_add2(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType_),
          space()->eInd(fType_),
          s_, A.s_,
          beta, alpha );
    else if( s_!=A.s_ && s_!=B.s_ )
      SF_add(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType_),
          space()->eInd(fType_),
          s_, A.s_, B.s_,
          alpha, beta );
    changed();
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void abs(const MV& y) {
    // add test for consistent VectorSpaces in debug mode
    SF_abs(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_),
        space()->eInd(fType_),
        s_,
        y.s_ );
		changed();
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal(const MV& y){
    // add test for consistent VectorSpaces in debug mode
    SF_reciprocal(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_),
        space()->eInd(fType_),
        s_,
        y.s_ );
    changed();
  }


  /// \brief Scale each element of the vector with \c alpha.
	void scale( const Scalar& alpha ) {
		SF_scale(
				space()->nLoc(),
				space()->bl(),
				space()->bu(),
				space()->sInd(fType_),
				space()->eInd(fType_),
				s_,
				alpha);
		changed();
	}


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
	void scale(const MV& a) {
		// add test for consistent VectorSpaces in debug mode
		SF_scale2(
				space()->nLoc(),
				space()->bl(),
				space()->bu(),
				space()->sInd(fType_),
				space()->eInd(fType_),
				s_, a.s_ );
		changed();
	}

  /// @}
  /// \name Norm method(reductions)
  /// @{

	/// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  Scalar dot ( const MV& a, bool global=true ) const {

    Scalar b = 0.;

		setCornersZero();
		a.setCornersZero();

    SF_dot(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_),
        space()->eInd(fType_),
        s_,
        a.s_,
        b );

    if( global ) this->reduceNorm( comm(), b );

    return( b );
  }




  /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
	/// \todo include scaled norm
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = 0.;

		setCornersZero();

    switch(type) {
    case Belos::OneNorm:
      SF_comp1Norm(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType_),
          space()->eInd(fType_),
          s_,
          normvec );
      break;
    case Belos::TwoNorm:
      SF_comp2Norm(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType_),
          space()->eInd(fType_),
          s_,
          normvec );
      break;
    case Belos::InfNorm:
      SF_compInfNorm(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType_),
          space()->eInd(fType_),
          s_,
          normvec );
      break;
    }

    if( global ) this->reduceNorm( comm(), normvec, type );

    return( normvec );
  }


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  double norm(const MV& weights, bool global=true ) const {

    Scalar normvec = 0.;

		setCornersZero();

    SF_weightedNorm(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_ ),
        space()->eInd(fType_),
        s_,
				weights.s_,
        normvec );

    if( global ) this->reduceNorm( comm(), normvec, Belos::TwoNorm );

    return( normvec );

  }


  //\}
  /// \name Initialization methods
  //\{

  /// \brief mv := A
  ///
  /// Assign (deep copy) \c a into \c this.
  /// total deep, boundaries and everythin.
  /// \note the \c StencilWidths is not take care of assuming every field is generated with one
  /// \note "indexing" is done c++
  void assign( const MV& a ) {

//		SF_assign(
//				space()->nLoc(),
//				space()->bl(),
//				space()->bu(),
//				space()->sInd(fType_ ),
//				space()->eInd(fType_),
//				s_, a.s_ );
//
//		changed();

	 for(int i=0; i<getStorageSize(); ++i)
		 s_[i] = a.s_[i];

	 for( int dir=0; dir<space()->dim(); ++dir )
		 exchangedState_[dir] = a.exchangedState_[dir];
  }


  /// \brief Replace the vectors with a random vectors.
  /// depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is required
  void random( bool useSeed = false, int seed = 1 ) {
    SF_random(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_),
        space()->eInd(fType_),
        s_);
		if( !space()->getProcGrid()->participating() )
			SF_init(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sIndB(fType_),
          space()->eIndB(fType_),
          s_,
          0. );
    changed();
  }


  /// \brief Replace each element of the vector  with \c alpha.
	/// \deprecated
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    SF_init(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_),
        space()->eInd(fType_),
        s_, alpha);
		if( !space()->getProcGrid()->participating() )
			SF_init(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sIndB(fType_),
          space()->eIndB(fType_),
          s_,
          0. );
    changed();
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
			const bool& Flow_Type_not_known = true; 
			TEUCHOS_TEST_FOR_EXCEPT( Flow_Type_not_known );
		}
		return( ConstField ); // just to please the compiler
	}

public:

	///  \brief initializes including boundaries to zero 
	void initField( Teuchos::ParameterList& para ) {

		EScalarField type =
			string2enum( para.get<std::string>( "Type", "constant" ) );

		switch( type ) {
			case ConstField :
				SF_init(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						s_,
						para.get<Scalar>( "C", 0. ) );
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
						para.get<Scalar>( "dx", 1. )	);
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
						para.get<Scalar>( "dy", 1. )	);
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
						para.get<Scalar>( "dz", 1. )	);
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
					para.get<Scalar>( "c_x", 1. ),
					para.get<Scalar>( "c_y", space()->getDomainSize()->getSize( Y )/2. ),
					para.get<Scalar>( "c_z", space()->getDomainSize()->getSize( Z )/2. ) };
				Scalar amp = para.get<Scalar>( "amp", 1. );
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
			SF_init(
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->sIndB(fType_),
					space()->eIndB(fType_),
					s_,
					0. );

		changed();
	}


	/// \brief initializes VectorField with the initial field defined in Fortran
	/// \deprecated
	void initField( EScalarField fieldType = ConstField, Scalar alpha=0. ) {
		switch( fieldType ) {
			case ConstField :
				SF_init(
						space()->nLoc(),
						space()->bl(),
						space()->bu(),
						space()->sIndB(fType_),
						space()->eIndB(fType_),
						s_,
						alpha );
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
						(std::abs(alpha)<1.e-16)?1.:alpha	);
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
						(std::abs(alpha)<1.e-16)?1.:alpha	);
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
						(std::abs(alpha)<1.e-16)?1.:alpha	);
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
					1.,
					//				1.,
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
			SF_init(
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->sIndB(fType_),
					space()->eIndB(fType_),
					s_,
					0. );
		changed();
	}


	/// \brief set corners of Dirichlet boundary conditions to zero
	void setCornersZero() const {

		if( EField::S == fType_ ) {
			SF_handle_corner(
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->getBCLocal()->getBCL(),
					space()->getBCLocal()->getBCU(),
					s_ );

			changed();
		}
	}


	void level() const {

		if( EField::S == fType_ ) {

			// set corners to zero, such that level depends only on inner field
			setCornersZero();

			SF_level(
					MPI_Comm_c2f( space()->comm() ),
					getLength(),
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->sIndB(fType_),
					space()->eIndB(fType_),
					s_ );

			//changed(); already done in setCornersZero
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

		for( Ordinal k=space()->sIndB(fType_,2); k<=space()->eIndB(fType_,2); ++k ) {
			for( Ordinal j=space()->sIndB(fType_,1); j<=space()->eIndB(fType_,1); ++j ) {
				for( Ordinal i=space()->sIndB(fType_,0); i<=space()->eIndB(fType_,0); ++i ) {
					out << i << "\t" << j << "\t" << k << "\t"
						<< at(i,j,k) << "\n";
						//<< s_[  (i-space()->bl(0)) +
										//(j-space()->bl(1))*cw[0] +
									 // (k-space()->bl(2))*cw[0]*cw[1] ] << "\n";
				}
			}
		}
  }


  /// Write the ScalarField to an hdf5 file, the velocities are interpolated to the pressure points
  /// \todo add 3d case here
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

			if( 2==space()->dim() ) {

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
						space()->getCoordinatesGlobal()->get(0,EField::S),
						space()->getCoordinatesGlobal()->get(1,EField::S),
						space()->getDomainSize()->getRe(),
						space()->getDomainSize()->getAlpha2() );
      }
      else if( 3==space()->dim() ) {

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
            space()->getCoordinatesGlobal()->get(0,EField::S),
            space()->getCoordinatesGlobal()->get(1,EField::S),
            space()->getCoordinatesGlobal()->get(2,EField::S),
            space()->getCoordinatesGlobal()->get(0,EField::U),
            space()->getCoordinatesGlobal()->get(1,EField::V),
            space()->getCoordinatesGlobal()->get(2,EField::W),
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
          space()->getCoordinatesGlobal()->get(0,EField::S),
          space()->getCoordinatesGlobal()->get(1,EField::S),
          space()->getCoordinatesGlobal()->get(2,EField::S),
          space()->getCoordinatesGlobal()->get(0,EField::U),
          space()->getCoordinatesGlobal()->get(1,EField::V),
          space()->getCoordinatesGlobal()->get(2,EField::W),
          space()->getDomainSize()->getRe(),
          space()->getDomainSize()->getAlpha2() );

    }

  }



public:

  const EField& getType() const { return( fType_ ); }

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

  const Teuchos::RCP<const SpaceT>& space() const { return( AbstractField<SpaceT>::space_ ); }

  const MPI_Comm& comm() const { return( space()->comm() ); }

  /// \name comunication methods.
	/// \brief highly dependent on underlying storage should only be used by
	/// Operator or on top field implementer.
  /// \{

  void changed( const int& dir ) const {
    exchangedState_[dir] = false;
  }
  void changed() const {
    for( int dir=0; dir<space()->dim(); ++dir )
      changed( dir );
  }


  bool is_exchanged( const int& dir ) const {
    return( exchangedState_[dir] );
  }
  bool is_exchanged() const {
    bool all_exchanged = true;
    for( int dir=0; dir<space()->dim(); ++dir )
      all_exchanged = all_exchanged && is_exchanged(dir);
    return( all_exchanged );
  }

  /// \brief updates ghost layers
  void exchange( const int& dir ) const {
		int ones[3] = {0,0,0};
    if( !exchangedState_[dir] ) {
      F_exchange(
          space()->dim(),
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
		for( int dir=0; dir<space()->dim(); ++dir )
			exchange( dir );
	}

  void setExchanged( const int& dir ) const {
    exchangedState_[dir] = true;
  }
  void setExchanged(  ) const {
    for( int dir=0; dir<space()->dim(); ++dir )
      changed( dir );
  }

  /// \}


	/// \brief indexing
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	///
	/// \return const reference
	inline const Scalar& at( const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {
		return( s_[ (i-space()->bl(0)) +
				        (j-space()->bl(1))*(space()->nLoc(0)+space()->bu(0)-space()->bl(0)+1) +
				        (k-space()->bl(2))*(space()->nLoc(0)+space()->bu(0)-space()->bl(0)+1)*(space()->nLoc(1)+space()->bu(1)-space()->bl(1)+1) ] );
	}

	/// \brief indexing
	///
	/// \param i index in x-direction
	/// \param j index in y-direction
	/// \param k index in z-direction
	///
	/// \return reference
	inline Scalar& at( const Ordinal& i, const Ordinal& j, const Ordinal& k )  {
		return( s_[ (i-space()->bl(0)) +
				        (j-space()->bl(1))*(space()->nLoc(0)+space()->bu(0)-space()->bl(0)+1) +
				        (k-space()->bl(2))*(space()->nLoc(0)+space()->bu(0)-space()->bl(0)+1)*(space()->nLoc(1)+space()->bu(1)-space()->bl(1)+1) ] );
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


} // end of namespace Pimpact


#ifdef COMPILE_ETI
#include "Pimpact_Space.hpp"
extern template class Pimpact::ScalarField< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::ScalarField< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::ScalarField< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::ScalarField< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_SCALARFIELD_HPP
