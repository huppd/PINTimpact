#pragma once
#ifndef PIMPACT_SCALARFIELD_HPP
#define PIMPACT_SCALARFIELD_HPP

#include <vector>
#include <iostream>
#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "BelosTypes.hpp"
//#include "Teuchos_ScalarTraitsDecl.hpp"
//#include "Teuchos_SerialDenseMatrix.hpp"

#include "Pimpact_Types.hpp"
#include "Pimpact_Space.hpp"

#include "Pimpact_extern_ScalarField.hpp"

#include "Pimpact_AbstractField.hpp"





namespace Pimpact {


/// \brief important basic Vector class
/// vector for a scalar field, e.g.: pressure,
/// \note all indexing is done in Fortran
/// \ingroup Field
template< class S=double, class O=int, int d=3 >
class ScalarField : private AbstractField<S,O,d> {

  template<class S1,class O1,int dimension1>
  friend class VectorField;
  template<class S1,class O1,int dimension1>
  friend class GradOp;
  template<class S1,class O1,int dimension1>
  friend class DivOp;
  template<class S1,class O1,int dimension1>
  friend class InterpolateV2S;
  template<class S1,class O1,int dimension1>
  friend class InterpolateS2V;
  template<class S1,class O1,int dimension1>
  friend class ConvectionSOp;
  template<class S1,class O1,int dimension1>
  friend class DivGradOp;
  template<class S1,class O1,int dimension1>
  friend class MGVDivGradOp;
  template<class S1,class O1,int dimension1>
  friend class RestrictionOp;
  template<class S1,class O1,int dimension1>
  friend class InterpolationOp;
  template<class Field>
  friend class TimeField;

public:

  typedef S Scalar;
  typedef O Ordinal;

  static const int dimension = d;

  typedef typename AbstractField<S,O,d>::SpaceT SpaceT;


protected:

  typedef Scalar* ScalarArray;
  typedef ScalarField<Scalar,Ordinal,dimension> MV;
  typedef Teuchos::Tuple<bool,3> State;

  ScalarArray s_;

  bool owning_;

  State exchangedState_;

  EField fType_;

public:


  ScalarField( const Teuchos::RCP<SpaceT>& space, bool owning, EField fType=EField::S ):
    AbstractField<S,O,d>( space ),
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
    AbstractField<S,O,d>( sF.space() ),
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
  ///@{

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


  //@}
  /// @name Update methods
  //@{

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
  /// \todo implement me
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
  /// \todo implement me
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
  /// \todo implement me
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
  /// \todo add test in debuging mode for testing equality of VectorSpaces
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

    if( global ) this->reduceNorm( space()->comm(), b );

    return( b );
  }


  ///@}
  /// @name Norm method
  ///@{


  /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  /// \todo implement OneNorm
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

    if( global ) this->reduceNorm( space()->comm(), normvec, type );

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

    if( global ) this->reduceNorm( space()->comm(), normvec, Belos::TwoNorm );

    return( normvec );

  }


  //@}
  /// @name Initialization methods
  //@{

  /// \brief mv := A
  ///
  /// Assign (deep copy) \c a into \c this.
  /// total deep, boundaries and everythin.
  /// \note the \c StencilWidths is not take care of assuming every field is generated with one
  /// \note "indexing" is done c++
  void assign( const MV& a ) {
    SF_assign(
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->sInd(fType_ ),
        space()->eInd(fType_),
        s_, a.s_ );
    changed();
    //    Ordinal N = 1;
    //    for(int i=0; i<3; ++i)
    //      N *= space()->nLoc(i)+space()->bu(i)-space()->bl(i);
    //
    //    for(int i=0; i<N; ++i)
    //      s_[i] = a.s_[i];
    //
    //    for( int dir=0; dir<dim(); ++dir )
    //      exchangedState_[dir] = a.exchangedState_[dir];
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
  void initField( EScalarField fieldType = Poiseuille2D_inX, double re=0. ) {
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


  //@}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os=std::cout )  const {
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
  void write( int count=0 ) {

    if( 0==space()->rankST() )
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

    Teuchos::RCP< ScalarField<S,O,d> > temp;

    if( EField::S != fType_ )
      temp = Teuchos::rcp(
          new ScalarField<S,O,d>( space(), true, EField::S ) );

    if( 2==space()->dim() ) {
      if( EField::S==fType_ ) {
        write_hdf5_2D(
            space()->rankST(),
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
            fType_,
            count,
            9,
            s_,
            space()->getCoordinatesGlobal()->get(0,EField::S),
            space()->getCoordinatesGlobal()->get(1,EField::S),
            space()->getDomain()->getDomainSize()->getRe(),
            space()->getDomain()->getDomainSize()->getAlpha2() );
      }
      else {

        space()->getInterpolateV2S()->apply( *this, *temp );

        write_hdf5_2D(
            space()->rankST(),
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
            10,
            temp->s_,
            space()->getCoordinatesGlobal()->get(0,EField::S),
            space()->getCoordinatesGlobal()->get(1,EField::S),
            space()->getDomain()->getDomainSize()->getRe(),
            space()->getDomain()->getDomainSize()->getAlpha2() );
      }
    }
    else if( 3==space()->dim() ) {
      if( EField::S==fType_ ) {

      }
      else {

      }
//      SF_write3D(
//          space()->nLoc(),
//          space()->bl(),
//          space()->bu(),
//          space()->sInd(fType_),
//          space()->eInd(fType_),
//          s_, count );
    }

  }



public:

  Ordinal getStorageSize() const {

    Ordinal n = 1;
    for(int i=0; i<3; ++i)
      n *= space()->nLoc(i)+space()->bu(i)-space()->bl(i)+1; // there a one was added for AMG, but it is not neede error seem to be in Impact there it should be (B1L+1:N1+B1U) probably has to be changed aganin for 3D

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

  Teuchos::RCP<SpaceT> space() const { return( AbstractField<S,O,d>::space_ ); }

protected:

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
    int ones[3] = {1,1,1};
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

}; // end of class ScalarField




/// \brief creates a scalar field(vector) belonging to a StencilWidths
///
/// \param fS scalar Vector Space to which returned vector belongs
/// \return scalar vector
/// \relates ScalarField
template<class S=double, class O=int, int d=3>
Teuchos::RCP< ScalarField<S,O,d> >
createScalarField(
    const Teuchos::RCP<const Space<S,O,d> >& space,
    EField fType=EField::S ) {

  return( Teuchos::rcp(
      new ScalarField<S,O,d>( space, true, fType ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SCALARFIELD_HPP
