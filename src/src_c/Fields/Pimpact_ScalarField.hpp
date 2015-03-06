#pragma once
#ifndef PIMPACT_SCALARFIELD_HPP
#define PIMPACT_SCALARFIELD_HPP

#include <vector>
#include <iostream>
#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "BelosTypes.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_extern_ScalarField.hpp"

#include "Pimpact_AbstractField.hpp"




namespace Pimpact {


/// \brief important basic Vector class
/// vector for a scalar field, e.g.: pressure,
/// \note all indexing is done in Fortran
/// \ingroup Field
/// \todo think about using Teuchos::ArrayRCP instead of Scalar* should make delete
template<class SpaceType>
class ScalarField : private AbstractField< SpaceType > {

  template<class SpaceTT>
  friend class DivGradO2JSmoother;
  template<class OperatorTT>
  friend class ConvectionDiffusionJSmoother;
  template<class Field>
  friend class TimeField;

public:

  typedef SpaceType SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

protected:

  typedef Scalar* ScalarArray;
  typedef ScalarField< SpaceT > MV;
  typedef Teuchos::Tuple<bool,3> State;

  ScalarArray s_;

  bool owning_;

  State exchangedState_;

  EField fType_;

public:


  ScalarField( const Teuchos::RCP<const SpaceT>& space, bool owning=true, EField fType=EField::S ):
    AbstractField<SpaceT>( space ),
    owning_(owning),
    exchangedState_( Teuchos::tuple(true,true,true) ),
    fType_(fType) {

    if( owning_ ) {

      Ordinal N = getStorageSize();

      s_ = new Scalar[N];

      for(int i=0; i<N; ++i)
        s_[i] = 0.;
    }
  };


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param sF
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
  ScalarField( const ScalarField& sF, ECopyType copyType=DeepCopy ):
    AbstractField<SpaceT>( sF.space() ),
    owning_( sF.owning_ ),
    exchangedState_( sF.exchangedState_ ),
    fType_( sF.fType_ ) {

    if( owning_ ) {

      Ordinal N = getStorageSize();

      s_ = new Scalar[N];

      switch( copyType ) {
      case ShallowCopy:
        for(int i=0; i<N; ++i)
          s_[i] = 0;
        break;
      case DeepCopy:
        for( int i=0; i<N; ++i)
          s_[i] = sF.s_[i];
        break;
      }
    }

  };

  ~ScalarField() {
    if( owning_ ) delete[] s_;
  }


  Teuchos::RCP<MV> clone( ECopyType ctype=DeepCopy ) const {
    return( Teuchos::rcp( new MV(*this, ctype) ) );
  }

  /// \name Attribute methods
  /// \{

  /// \brief returns the length of Field.
  Ordinal getLength( bool dummy=false ) const {

    //    auto bc = space_->getDomain()->getBCGlobal();
    //    auto bc = AbstractField<S,O,d>::space_->getDomain()->getBCGlobal();
    //        auto bc = this->space_->getDomain()->getBCGlobal();
    auto bc = space()->getDomain()->getBCGlobal();

    Ordinal vl = 1;

    switch( fType_ ) {
    case EField::S: {
      for(int i = 0; i<space()->dim(); ++i)
        if( PeriodicBC==bc->getBCL(i) )
          vl *= space()->nGlo(i)-1;
        else
          vl *= space()->nGlo(i);
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


  /// \}
  /// \name Update methods
  /// \{

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


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  Scalar dot ( const MV& a, bool global=true ) const {
    Scalar b = 0.;

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


  ///\}
  /// \name Norm method
  ///\{


  /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = 0.;

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

    SF_weightedNorm(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_ ),
        space()->eInd(fType_),
        s_, weights.s_,
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

//    SF_assign(
//        space()->nLoc(),
//        space()->bl(),
//        space()->bu(),
//        space()->sInd(fType_ ),
//        space()->eInd(fType_),
//        s_, a.s_ );
//
//    changed();


    for(int i=0; i<getStorageSize(); ++i)
      s_[i] = a.s_[i];

    for( int dir=0; dir<space()->dim(); ++dir )
      exchangedState_[dir] = a.exchangedState_[dir];
  }


  /// \brief Replace the vectors with a random vectors.
  /// depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is requiered
  void random( bool useSeed = false, int seed = 1 ) {
    SF_random(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_),
        space()->eInd(fType_),
        s_);
    changed();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    SF_init(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_),
        space()->eInd(fType_),
        s_, alpha);
    changed();
  }


  ///  \brief initializes VectorField with the initial field defined in Fortran
  void initField( EScalarField fieldType = ConstField, Scalar re=0. ) {
    switch( fieldType ) {
    case ConstField :
      SF_init(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sIndB(fType_),
          space()->eIndB(fType_),
          s_,
          re );
      break;
    case Poiseuille2D_inX :
      SF_init_2DPoiseuilleX(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sIndB(fType_),
          space()->eIndB(fType_),
          space()->getDomain()->getDomainSize()->getSize( Y ),
          space()->getCoordinatesLocal()->getX(Y,fType_),
          s_ );
      break;
    case Poiseuille2D_inY :
      SF_init_2DPoiseuilleY(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sIndB(fType_),
          space()->eIndB(fType_),
          space()->getDomain()->getDomainSize()->getSize( X ),
          space()->getCoordinatesLocal()->getX(X,fType_),
          s_ );
      break;
    case Grad2D_inX :
      SF_init_2DGradX(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sIndB(fType_),
          space()->eIndB(fType_),
          space()->getDomain()->getDomainSize()->getSize(),
          space()->getCoordinatesLocal()->getX( X, fType_ ),
          s_ );
      break;
    case Grad2D_inY :
      SF_init_2DGradY(
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sIndB(fType_),
          space()->eIndB(fType_),
          space()->getDomain()->getDomainSize()->getSize(),
          space()->getCoordinatesLocal()->getX( Y, fType_ ),
          s_ );
      break;
    }
    changed();
  }


	void level() {
		// set corners to zero, such that level depends only on inner field
//		SF_handle_corner(
//				space()->nLoc(),
//				space()->bl(),
//				space()->bu(),
//				space()->getDomain()->getBCLocal()->getBCL(),
//				space()->getDomain()->getBCLocal()->getBCU(),
//				s_ );

		SF_level(
				space()->commf(),
				getLength(),
				space()->nLoc(),
				space()->bl(),
				space()->bu(),
				space()->sIndB(fType_),
				space()->eIndB(fType_),
				s_ );

//		SF_handle_corner(
//				space()->nLoc(),
//				space()->bl(),
//				space()->bu(),
//				space()->getDomain()->getBCLocal()->getBCL(),
//				space()->getDomain()->getBCLocal()->getBCU(),
//				s_ );
	}


  /// \}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& out=std::cout )  const {

    out << "--- FieldType: " << fType_ << "--- \n";
    out << "--- StorageSize: " << getStorageSize() << "---\n";
    out << "--- owning: " << owning_ << "---\n";
    out << "--- exchangedState: " << exchangedState_ << "--\n";

    SF_print(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sIndB(fType_),
        space()->eIndB(fType_),
        s_ );

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
            space()->commf(),
            space()->nGlo(),
            space()->getDomain()->getBCGlobal()->getBCL(),
            space()->getDomain()->getBCGlobal()->getBCU(),
            space()->nLoc(),
            space()->bl(),
            space()->bu(),
            space()->sInd(EField::S),
            space()->eInd(EField::S),
            space()->getStencilWidths()->getLS(),
            space()->getProcGridSize()->get(),
            space()->getProcGrid()->getIB(),
            space()->getProcGrid()->getShift(),
            (int)fType_,
            count,
            (EField::S==fType_)?9:10,
                (EField::S==fType_)?s_:temp->s_,
                    space()->getCoordinatesGlobal()->get(0,EField::S),
                    space()->getCoordinatesGlobal()->get(1,EField::S),
                    space()->getDomain()->getDomainSize()->getRe(),
                    space()->getDomain()->getDomainSize()->getAlpha2() );
      }
      else if( 3==space()->dim() ) {

        int stride[3] = {1,1,1};

        write_hdf_3D(
            space()->rankS(),
            space()->commf(),
            space()->nGlo(),
            space()->getDomain()->getBCGlobal()->getBCL(),
            space()->getDomain()->getBCGlobal()->getBCU(),
            space()->nLoc(),
            space()->bl(),
            space()->bu(),
            space()->sInd(EField::S),
            space()->eInd(EField::S),
            space()->getStencilWidths()->getLS(),
            space()->getProcGridSize()->get(),
            space()->getProcGrid()->getIB(),
            space()->getProcGrid()->getShift(),
            (int)fType_+1,
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
            space()->getDomain()->getDomainSize()->getRe(),
            space()->getDomain()->getDomainSize()->getAlpha2() );

      }
    }
    else {

      int stride[3] = {1,1,1};

      write_hdf_3D(
          space()->rankS(),
          space()->commf(),
          space()->nGlo(),
          space()->getDomain()->getBCGlobal()->getBCL(),
          space()->getDomain()->getBCGlobal()->getBCU(),
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->sInd(fType_),
          space()->eInd(fType_),
          space()->getStencilWidths()->getLS(),
          space()->getProcGridSize()->get(),
          space()->getProcGrid()->getIB(),
          space()->getProcGrid()->getShift(),
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
          space()->getDomain()->getDomainSize()->getRe(),
          space()->getDomain()->getDomainSize()->getAlpha2() );

    }

  }



public:

  const EField& getType() const { return( fType_ ); }

   /// \name storage methods.
   /// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
   ///
   ///\{

  Ordinal getStorageSize() const {

    Ordinal n = 1;
    for(int i=0; i<3; ++i)
      n *= space()->nLoc(i)+space()->bu(i)-space()->bl(i)+1; // seems wrong: there a one was added for AMG, but it is not neede error seem to be in Impact there it should be (B1L+1:N1+B1U) probably has to be changed aganin for 3D

    return( n );
  }

  void setStoragePtr( Scalar*  array ) {
    s_ = array;
  }

  ScalarArray getRawPtr() {
    return( s_ );
  }

  const Scalar* getConstRawPtr() const {
    return( s_ );
  }


  ///\}

  Teuchos::RCP<const SpaceT> space() const { return( AbstractField<SpaceT>::space_ ); }

  const MPI_Comm& comm() const { return(space()->comm()); }

  /// \name comunication methods.
  /// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
  ///
  ///\{
//protected:

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
//    int ones[3] = {1,1,1};
    if( !exchangedState_[dir] ) {
      F_exchange(
          space()->dim(),
          space()->commf(),
          space()->getProcGrid()->getRankL(),
          space()->getProcGrid()->getRankU(),
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->getDomain()->getBCLocal()->getBCL(),
          space()->getDomain()->getBCLocal()->getBCU(),
          space()->sInd(EField::S),
          space()->eInd(EField::S),
          space()->sIndB(fType_), // should it work
//          space()->eIndB(fType_),
//          ones,
          space()->nLoc(),
          1+dir,
          1+(int)fType_,
          s_);
      exchangedState_[dir] = true;
    }
  }

  void exchange( bool forward=true ) const {
    if(forward)
      for( int dir=0; dir<space()->dim(); ++dir )
        exchange( dir );
    else
      for( int dir=space()->dim()-1; dir>=0; --dir )
        exchange( dir );
  }

  ///\}

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


#endif // end of #ifndef PIMPACT_SCALARFIELD_HPP
