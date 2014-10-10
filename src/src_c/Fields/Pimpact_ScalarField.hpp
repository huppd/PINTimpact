#pragma once
#ifndef PIMPACT_SCALARFIELD_HPP
#define PIMPACT_SCALARFIELD_HPP

#include <vector>
#include <iostream>
#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

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
class ScalarField : private AbstractField<S,O> {

  template<class S1,class O1,int dimension1>
  friend class VectorField;
  template<class S1,class O1,int dimension1>
  friend class GradOp;
  template<class S1,class O1,int dimension1>
  friend class DivOp;
  template<class S1,class O1,int dimension1>
  friend class InterpolateV2S;
  template<class S1,class O1,int dimension1>
  friend class DivGradOp;
  template<class S1,class O1,int dimension1>
  friend class MGVDivGradOp;
  template<class S1,class O1,int dimension1>
  friend class RestrictionOp;
  template<class S1,class O1,int dimension1>
  friend class InterpolationOp;

public:

  typedef S Scalar;
  typedef O Ordinal;

  static const int dimension = d;

  typedef Space<Scalar,Ordinal,dimension> SpaceT;


protected:

  typedef Scalar* ScalarArray;
  typedef ScalarField<Scalar,Ordinal,dimension> MV;
  typedef Teuchos::Tuple<bool,3> State;

  Teuchos::RCP<const SpaceT > space_;

  ScalarArray s_;

  bool owning_;

  State exchangedState_;

  EField fType_;

public:

  ScalarField( EField fType=EField::S ):
    s_(0),
    space_(Teuchos::null),
    owning_(true),
    exchangedState_( Teuchos::tuple(true,true,true) ),
    fType_(fType) {};

  ScalarField( const Teuchos::RCP<const SpaceT >& space, bool owning=true, EField fType=EField::S ):
    space_(space),
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
    space_(sF.space_),
    owning_(sF.owning_),
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

    auto bc = space_->getDomain()->getBCGlobal();

    Ordinal vl = 1;

    switch( fType_ ) {
    case EField::S: {
      for(int i = 0; i<dim(); ++i)
        if( PeriodicBC==bc->getBCL(i) )
          vl *= nGlo(i)-1;
        else
          vl *= nGlo(i);
      break;
    }
    default: {
      for( int j=0; j<dim(); ++j) {
        if( fType_==j ) {
          vl *= nGlo(j)-1;
        }
        else {
          if( PeriodicBC==bc->getBCL(j) )
            vl *= nGlo(j)-2+1;
          else
            vl *= nGlo(j)-2;
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
          nLoc(), bl(), bu(),
          sInd(), eInd(),
          s_, B.s_,
          alpha, beta );
    else if( s_!=A.s_ && s_==B.s_ )
      SF_add2(
          nLoc(), bl(), bu(),
          sInd(), eInd(),
          s_, A.s_,
          beta, alpha );
    else if( s_!=A.s_ && s_!=B.s_ )
      SF_add(
          nLoc(), bl(), bu(),
          sInd(), eInd(),
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
        nLoc(),
        bl(), bu(),
        sInd(), eInd(),
        s_, y.s_ );
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
        nLoc(),
        bl(), bu(),
        sInd(), eInd(),
        s_, y.s_ );
    changed();
  }


  /// \brief Scale each element of the vector with \c alpha.
  void scale( const Scalar& alpha ) {
    SF_scale(
        nLoc(),
        bl(), bu(),
        sInd(), eInd(),
        s_, alpha);
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
        nLoc(),
        bl(), bu(),
        sInd(), eInd(),
        s_, a.s_ );
    changed();
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  /// \todo add test in debuging mode for testing equality of VectorSpaces
  Scalar dot ( const MV& a, bool global=true ) const {
    Scalar b = 0.;

    SF_dot(
        nLoc(),
        bl(), bu(),
        sInd(), eInd(),
        s_, a.s_, b );

    if( global ) this->reduceNorm( comm(), b );

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
          nLoc(),
          bl(), bu(),
          sInd(), eInd(),
          s_,
          normvec );
      break;
    case Belos::TwoNorm:
      SF_comp2Norm(
          nLoc(),
          bl(), bu(),
          sInd(), eInd(),
          s_,
          normvec );
      break;
    case Belos::InfNorm:
      SF_compInfNorm(
          nLoc(),
          bl(), bu(),
          sInd(), eInd(),
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
        nLoc(),
        bl(), bu(),
        sInd(), eInd(),
        s_, weights.s_,
        normvec );

    if( global ) this->reduceNorm( comm(), normvec, Belos::TwoNorm );

    return( normvec );

  }


  //@}
  /// @name Initialization methods
  //@{

  /// \brief mv := A
  ///
  /// Assign (deep copy) \c a into \c this.
  /// total deep, boundaries and everythin.
  /// \note the \c FieldSpace is not take care of assuming every field is generated with one
  /// \note "indexing" is done c++
  void assign( const MV& a ) {
    SF_assign(
        nLoc(),
        bl(), bu(),
        sInd(), eInd(),
        s_, a.s_ );
    changed();
    //    Ordinal N = 1;
    //    for(int i=0; i<3; ++i)
    //      N *= nLoc(i)+bu(i)-bl(i);
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
        nLoc(),
        bl(), bu(),
        sInd(), eInd(),
        s_);
    changed();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    SF_init(
        nLoc(),
        bl(), bu(),
        sInd(), eInd(),
        s_, alpha);
    changed();
  }


  ///  \brief initializes VectorField with the initial field defined in Fortran
  void initField( EScalarField fieldType = Poiseuille2D_inX, double re=1., double om=1., double px = 1. ) {
    switch( fieldType ) {
    case ZeroField :
      SF_init(
          nLoc(),
          bl(),
          bu(),
          sIndB(),
          eIndB(),
          s_,
          0. );
      break;
    case Poiseuille2D_inX :
      SF_init_2DPoiseuilleX(
          nLoc(),
          bl(),
          bu(),
          sIndB(),
          eIndB(),
          space_->getDomain()->getDomainSize()->getSize( Y ),
          space_->getCoordinatesLocal()->getX(Y,fType_),
          s_ );
      break;
    case Poiseuille2D_inY :
      SF_init_2DPoiseuilleY(
          nLoc(),
          bl(),
          bu(),
          sIndB(),
          eIndB(),
          space_->getDomain()->getDomainSize()->getSize( X ),
          space_->getCoordinatesLocal()->getX(X,fType_),
          s_ );
      break;
    case Grad2D_inX :
      SF_init_2DGradX(
          nLoc(),
          bl(),
          bu(),
          sIndB(),
          eIndB(),
          space_->getDomain()->getDomainSize()->getSize(),
          space_->getCoordinatesLocal()->getX( X, fType_ ),
          s_ );
      break;
    case Grad2D_inY :
      SF_init_2DGradY(
          nLoc(),
          bl(),
          bu(),
          sIndB(),
          eIndB(),
          space_->getDomain()->getDomainSize()->getSize(),
          space_->getCoordinatesLocal()->getX( Y, fType_ ),
          s_ );
      break;
    }
    changed();
  }


  //@}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os=std::cout )  const {
    SF_print(
        nLoc(),
        bl(), bu(),
        sIndB(),
        eIndB(),
        s_ );

  }


  void write( int count=0 ) {
    if( 0==space_->rankST() )
      std::cout << "writing pressure field (" << count << ") ...\n";
    // exchange?
    if( 2==space_->dim() )
      write_hdf5_2D(
          space_->commf(),
          space_->nGlo(),
          space_->getDomain()->getBCGlobal()->getBCL(),
          space_->getDomain()->getBCGlobal()->getBCU(),
          space_->nLoc(),
          space_->bl(),
          space_->bu(),
          sInd(),
          eInd(),
          space_->getFieldSpace()->getLS(),
          space_->getProcGridSize()->get(),
          space_->getProcGrid()->getIB(),
          space_->getProcGrid()->getShift(),
          fType_,
          count,
          s_,
          space_->getCoordinatesGlobal()->get(0,EField::S),
          space_->getCoordinatesGlobal()->get(1,EField::S),
          space_->getCoordinatesGlobal()->get(2,EField::S),
          space_->getDomain()->getDomainSize()->getRe(),
          space_->getDomain()->getDomainSize()->getAlpha2() );
//      SF_write2D(
//          nLoc(),
//          bl(), bu(),
//          sInd(),
//          eInd(),
//          s_, count );
    if( 3==space_->dim() )
      SF_write3D(
          nLoc(),
          bl(), bu(),
          sInd(),
          eInd(),
          s_, count );
  }



public:

  const MPI_Fint& commf() const { return( space_->commf() ); }
  MPI_Comm        comm()  const { return( space_->comm() ); }
  const int&      dim()   const { return( space_->dim()   ); }

  Teuchos::RCP<const SpaceT> getSpace() const { return( space_ ); };

  Ordinal getStorageSize() const {

    Ordinal n = 1;
    for(int i=0; i<3; ++i)
      n *= nLoc()[i]+bu()[i]-bl()[i]+1; // there a one was added for AMG, but it is not neede error seem to be in Impact there it should be (B1L+1:N1+B1U) probably has to be changed aganin for 3D

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

protected:

//  const Ordinal* nGlo()  const { return( space_->nGlo() ); }
  const Ordinal& nGlo(int i)  const { return( space_->nGlo()[i] ); }

  const Ordinal* nLoc()      const { return( space_->nLoc() ) ; }
  const Ordinal& nLoc(int i) const { return( space_->nLoc(i) ) ; }

  const Ordinal* bl() const { return( space_->bl() ); }
  const Ordinal& bl(int i) const { return( space_->bl(i) ); }

  const Ordinal* bu() const { return( space_->bu() ); }
  const Ordinal& bu(int i) const { return( space_->bu(i) ); }

  const Ordinal* sInd() const { return( space_->sInd( (int)fType_ ) ); }
  const Ordinal& sInd(int i) const { return( space_->sInd( (int)fType_ )[i] ); }

  const Ordinal* eInd() const { return( space_->eInd( (int)fType_ ) ); }
  const Ordinal& eInd(int i) const { return( space_->eInd( (int)fType_ )[i] ); }

  const Ordinal* sIndB() const { return( space_->sIndB( (int)fType_ ) ); }
  const Ordinal& sIndB(int i) const { return( space_->sIndB( (int)fType_ )[i] ); }

  const Ordinal* eIndB() const { return( space_->eIndB( (int)fType_ ) ); }
  const Ordinal& eIndB(int i) const { return( space_->eIndB( (int)fType_ )[i] ); }

  const int*     bcL() const { return( space_->getDomain()->getBCLocal()->getBCL() ); }
  const int*     bcU() const { return( space_->getDomain()->getBCLocal()->getBCU() ); }

  const int* rankL() const { return( space_->getProcGrid()->getRankL() ); }
  const int* rankU() const { return( space_->getProcGrid()->getRankU() ); }

  void changed( const int& dir ) const {
    exchangedState_[dir] = false;
  }

public:

  void changed() const {
    for( int dir=0; dir<dim(); ++dir )
      changed( dir );
  }

protected:

  bool is_exchanged( const int& dir ) const {
    return( exchangedState_[dir] );
  }
  bool is_exchanged() const {
    bool all_exchanged = true;
    for( int dir=0; dir<dim(); ++dir )
      all_exchanged = all_exchanged && is_exchanged(dir);
    return( all_exchanged );
  }

  /// \brief updates ghost layers
  void exchange( const int& dir ) const {
    int ones[3] = {1,1,1};
    if( !exchangedState_[dir] ) {
      F_exchange(
          dim(),
          commf(),
          rankL(), rankU(),
          nLoc(),
          bl(), bu(),
          bcL(), bcU(),
          sInd(), eInd(),
          ones,
          nLoc(),
          dir+1, 0,
          s_);
      exchangedState_[dir] = true;
    }
  }
  void exchange() const {
    for( int dir=0; dir<dim(); ++dir )
      exchange( dir );
  }

}; // end of class ScalarField




/// \brief creates a scalar field(vector) belonging to a FieldSpace
///
/// \param fS scalar Vector Space to which returned vector belongs
/// \return scalar vector
/// \relates ScalarField
template<class S=double, class O=int, int d=3>
Teuchos::RCP< ScalarField<S,O,d> >
createScalarField(
    const Teuchos::RCP<const Space<S,O,d> >& space=Teuchos::null,
    EField fType=EField::S ) {

  return( Teuchos::rcp(
      new ScalarField<S,O,d>( space, fType ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SCALARFIELD_HPP
