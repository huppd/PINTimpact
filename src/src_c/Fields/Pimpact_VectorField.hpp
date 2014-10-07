#pragma once
#ifndef PIMPACT_VECTORFIELD_HPP
#define PIMPACT_VECTORFIELD_HPP

#include <vector>
#include <iostream>
#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Pimpact_Space.hpp"

#include "Pimpact_extern_ScalarField.hpp"
#include "Pimpact_extern_VectorField.hpp"

#include "Pimpact_AbstractField.hpp"

#include "Pimpact_ScalarField.hpp"


namespace Pimpact {


/// \brief important basic Vector class
/// vector for a vector field, e.g.: velocity,
/// here also happens the fortran wrapping
/// \ingroup Field
/// \todo There is one issue: update methods changes state but could be
///implemented, such that they keep the state, but then the boundary conditions
///have to be taken care of
/// \todo move exhange to \c ScalarField
template<class S=double, class O=int, int d=3 >
class VectorField : AbstractField<S,O> {

  template<class S1, class O1, int dimension1>
  friend class Grad;
  template<class S1,class O1, int dimension1>
  friend class Div;
  template<class S1,class O1, int dimension1>
  friend class HelmholtzOp;
  template<class S1,class O1, int d1>
  friend class Nonlinear;
  template<class S1,class O1, int dimension1>
  friend class NonlinearJacobian;
  template<class S1,class O1>
  friend class DtLapOp;
  template<class S1,class O1>
  friend class MLHelmholtzOp;
  template<class S1,class O1>
  friend class InverseHelmholtzOp;
  template<class S1,class O1,int dimension1>
  friend class MGVHelmholtzOp;
  template< class S1, class O1, bool CNY >
  friend class TimeNonlinearJacobian;


public:

  typedef S Scalar;
  typedef O Ordinal;

  static const int dimension = d;

  typedef Space<Scalar,Ordinal,dimension> SpaceT;

  typedef Teuchos::Tuple<Teuchos::Tuple<bool,3>,3> State; // obsolte in ScalarField

protected:

  typedef Scalar* ScalarArray;
  typedef VectorField<Scalar,Ordinal,dimension> VF;
  typedef ScalarField<Scalar,Ordinal,dimension> SF;

  Teuchos::RCP< const SpaceT > space_;

  //  Teuchos::Tuple<ScalarArray,3> vec_;
  ScalarArray vec_;

  const bool owning_;

  State exchangedState_; // obsolete > SF

  Teuchos::Tuple< Teuchos::RCP<SF>, 3 > sFields_;

public:

  VectorField():
    space_(Teuchos::null),
    vec_(0),
    owning_(true),
    exchangedState_(Teuchos::tuple(Teuchos::tuple(true,true,true),Teuchos::tuple(true,true,true),Teuchos::tuple(true,true,true))),
    sFields_( Teuchos::tuple(Teuchos::null,Teuchos::null,Teuchos::null) )
{};

  VectorField( const Teuchos::RCP< const SpaceT >& space, bool owning=true ):
    space_(space),
    owning_(owning),
    exchangedState_(Teuchos::tuple(Teuchos::tuple(true,true,true),Teuchos::tuple(true,true,true),Teuchos::tuple(true,true,true)))//,
  //    sFields_( Teuchos::tuple(Teuchos::null,Teuchos::null,Teuchos::null) )
  {

    for( int i=0; i<3; ++i )
      sFields_[i] = Teuchos::rcp( new SF( space_, false, EField(i) ) );

    if( owning_ ) {

      Ordinal N = getStorageSize()/3;

      vec_ = new Scalar[3*N];
      //      vec_[0] = new Scalar[3*N];
      //      vec_[1] = vec_[0]+N; // shoudl become obsolete
      //      vec_[2] = vec_[1]+N;

      //      for(int i=0; i<3; ++i)
      //        for(int j=0; j<N; ++j)
      //          vec_[i][j] = 0.;
      for( int i=0; i<3*N; ++i )
        vec_[i] = 0.;

      for( int i=0; i<3; ++i )
        sFields_[i]->setStoragePtr( vec_+i*N );
    }

  };


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
  VectorField(const VectorField& vF, ECopyType copyType=DeepCopy):
    space_(vF.space_),
    owning_(vF.owning_),
    exchangedState_(vF.exchangedState_) {

    for( int i=0; i<3; ++i )
      sFields_[i] = Teuchos::rcp( new SF( space_, false, EField(i) ) );

    if( owning_ ) {

      Ordinal n = getStorageSize()/3;

      vec_ = new Scalar[3*n];

      for( int i=0; i<3; ++i )
        sFields_[i]->setStoragePtr( vec_+i*n );

      switch( copyType ) {
      case ShallowCopy:
        for( int i=0; i<3*n; ++i )
          vec_[i] = 0.;
        break;
      case DeepCopy:
        for( int i=0; i<dim(); ++i )
          sFields_[i]->assign( *vF.sFields_[i] );
        changed();
        break;
      }
    }
  };


  ~VectorField() {
    if( owning_ ) delete[] vec_;
  }

  Teuchos::RCP<VF> clone( ECopyType ctype=DeepCopy ) const {
    return( Teuchos::rcp( new VF( *this, ctype ) ) );
  }

  /// \name Attribute methods
  //@{


  /// \brief returns the length of Field.
  ///
  /// the vector length is withregard to the inner points such that
  /// \f[ N_u = (N_x-1)(N_y-2)(N_z-2) \f]
  /// \f[ N_v = (N_x-2)(N_y-1)(N_z-2) \f]
  /// \f[ N_w = (N_x-2)(N_y-2)(N_z-1) \f]
  /// \return vect length \f[= N_u+N_v+N_w\f]
  Ordinal getLength( bool dummy=false ) const {

    auto bc = space_->getDomain()->getBCGlobal();

    Ordinal n = 0;
    for( int i=0; i<dim(); ++i )
      n += sFields_[i]->getLength( dummy );

    return( n );
  }


  /// \brief get number of stored Field's
  int getNumberVecs() const { return( 1 ); }


  //@}
  /// \name Update methods
  //@{

  /// \brief Replace \c this with \f$\alpha A + \beta B\f$.
  ///
  /// only inner points
  void add( const Scalar& alpha, const VF& A, const Scalar& beta, const VF& B ) {
    // add test for consistent VectorSpaces in debug mode
    for( int i=0; i<dim(); ++i ) {
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
  void abs(const VF& y) {
    for( int i=0; i<dim(); ++i )
      sFields_[i]->abs( *y.sFields_[i] );
    changed();
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal(const VF& y){
    // add test for consistent VectorSpaces in debug mode
    for( int i=0; i<dim(); ++i)
      sFields_[i]->reciprocal( *y.sFields_[i] );
    changed();
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale( const Scalar& alpha ) {
    for(int i=0; i<dim(); ++i)
      sFields_[i]->scale( alpha );
    changed();
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void scale(const VF& a) {
    // add test for consistent VectorSpaces in debug mode
    for(int i=0; i<dim(); ++i)
      sFields_[i]->scale( *a.sFields_[i] );
    changed();
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  /// \todo add test in debuging mode for testing equality of VectorSpaces
  Scalar dot ( const VF& a, bool global=true ) const {
    Scalar b = 0.;

    for( int i=0; i<dim(); ++i )
      b += sFields_[i]->dot( *a.sFields_[i], false );

    if( global ) this->reduceNorm( comm(), b );

    return( b );

  }


  //@}
  /// @name Norm method
  //@{

  /// \brief Compute the norm of each individual vector.
  ///
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2^2\f$, the \c i-th column of \c this.
  /// \attention the two norm is not the real two norm but its square
  /// \todo implement OneNorm
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = 0.;

    for( int i=0; i<dim(); ++i )
      switch(type) {
      default:
        normvec += sFields_[i]->norm(type,false);
        break;
      case Belos::InfNorm:
        normvec = std::max( sFields_[i]->norm(type,false), normvec ) ;
        break;
      }

    if( global ) this->reduceNorm( comm(), normvec, type );

    return( normvec );

  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  double norm( const VF& weights, bool global=true ) const {
    Scalar normvec = 0.;

    for( int i=0; i<dim(); ++i )
      normvec += sFields_[i]->norm( *weights.sFields_[i], false);

    if( global ) this->reduceNorm( comm(), normvec, Belos::TwoNorm );

    return( normvec );

  }


  //@}
  /// \name Initialization methods
  //@{


  /// \brief mv := A
  ///
  /// Assign (deep copy) A into mv.
  void assign( const VF& a ) {

    for( int i=0; i<dim(); ++i)
      sFields_[i]->assign( *a.sFields_[i] );
    changed();
  }


  /// \brief Replace the vectors with a random vectors.
  ///
  /// depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is requiered
  void random(bool useSeed = false, int seed = 1) {

    for( int i=0; i<dim(); ++i )
      sFields_[i]->random( useSeed, seed );

    changed();
  }

  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    for( int i=0; i<dim(); ++i )
      sFields_[i]->init( alpha );
    changed();
  }


  /// \brief Replace each element of the vector \c vec[i] with \c alpha[i].
  void init( const Teuchos::Tuple<Scalar,3>& alpha ) {
    for( int i=0; i<dim(); ++i )
      sFields_[i]->init( alpha[i] );
    changed();
  }


  ///  \brief initializes VectorField with the initial field defined in Fortran
  void initField( EFlowProfile flowType = Poiseuille2D_inX, double re=1., double om=1., double px = 1. ) {
    switch( flowType) {
    case ZeroProf :
      for( int i=0; i<dim(); ++i )
        SF_init(
            nLoc(),
            bl(),
            bu(),
            sIndB(i),
            eIndB(i),
            vec(i),
            0. );
      break;
    case Poiseuille2D_inX :
      for( int i=0; i<dim(); ++i )
        if( U==i )
          VF_init_2DPoiseuilleX(
              nLoc(),
              bl(),
              bu(),
              sIndB(i),
              eIndB(i),
              space_->getDomain()->getDomainSize()->getSize( Y ),
              space_->getCoordinatesLocal()->getX(Y,i),
              vec(i) );
        else
          SF_init(
              nLoc(),
              bl(),
              bu(),
              sIndB(i),
              eIndB(i),
              vec(i),
              0. );
      break;
    case Poiseuille2D_inY :
      for( int i=0; i<dim(); ++i )
        if( V==i )
          VF_init_2DPoiseuilleY(
              nLoc(),
              bl(),
              bu(),
              sIndB(i),
              eIndB(i),
              space_->getDomain()->getDomainSize()->getSize( X ),
              space_->getCoordinatesLocal()->getX(X,i),
              vec(i) );
        else
          SF_init(
              nLoc(),
              bl(),
              bu(),
              sIndB(i),
              eIndB(i),
              vec(i),
              0. );
      break;
    case Pulsatile2D_inXC :
      VF_init_2DPulsatileXC(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re, om, px);
      break;
    case Pulsatile2D_inYC :
      VF_init_2DPulsatileYC(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re, om, px);
      break;
    case Pulsatile2D_inXS :
      VF_init_2DPulsatileXS(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re, om, px);
      break;
    case Pulsatile2D_inYS:
      VF_init_2DPulsatileYS(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re, om, px);
      break;
    case Streaming2D:
      VF_init_StreamingS(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re );
      break;
    case Streaming2DC:
      VF_init_StreamingC(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re );
      break;
    case Streaming2DS:
      VF_init_StreamingS(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re );
      break;
    case Circle2D:
      VF_init_Circle(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr() );
      break;
    case RankineVortex2D:
      VF_init_RankineVortex(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr() );
      break;
    case GaussianForcing1D:
      VF_init_GaussianForcing1D(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr() );
      break;
    case BoundaryFilter1D:
      VF_init_BoundaryFilter1D(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr() );
      break;
    case GaussianForcing2D:
      VF_init_GaussianForcing2D(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr() );
      break;
    case BoundaryFilter2D:
      VF_init_BoundaryFilter2D(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr() );
      break;
    case VPoint2D:
      VF_init_Vpoint(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re );
      break;
    case Disc2D:
      VF_init_Disc(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re, om, px );
      break;
    case RotationDisc2D:
      VF_init_RotatingDisc(
          nLoc(0), nLoc(1), nLoc(2),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(),
          re, om, px );
      break;
    }
    //    }
    changed();
  }


  //@}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os=std::cout )  {
    int rank;
    MPI_Comm_rank(comm(),&rank);
    for(int i=0; i<dim(); ++i) {
      os << "rank: " << rank << " :dir: " << i << "\n";
      os << "rank: " << rank << " :dirs: " << sFields_[i]->fType_ << "\n";
      os << "rank: " << rank << " :nGlo: " << nGlo(i) << "\n";
      os << "rank: " << rank << " :nLoc: " << nLoc(i) << "\n";
      for( int j=0; j<3; ++j ) {
        os << "rank: " << rank << "field: " << j << " :sInd: " << sInd(i,j) << "\n";
        os << "rank: " << rank << "field: " << j << " :eInd: " << eInd(i,j) << "\n";
      }
      os << "rank: " << rank << " :bl: " << bl(i) << "\n";
      os << "rank: " << rank << " :bu: " << bu(i) << "\n\n";
    }
    std::cout << "rank: " << rank << "\n";
    //    for( int i=0; i<dim(); ++i ) {
    //      std::cout << "field: " << i << "\n";
    //      SF_print(
    //          nLoc(),
    //          bl(), bu(),
    //          sInd(i), eInd(i),
    //          vec_[i] );
    //    }
  }


  void write( int count=0 ) {
    exchange();
    VF_write( sFields_[U]->getRawPtr(), sFields_[V]->getRawPtr(), sFields_[W]->getRawPtr(), count );
  }


public:

  /// \todo add good documetnation here
  /// @return
  const MPI_Fint& commf() const { return( space_->commf() ); }
  MPI_Comm        comm()  const { return( space_->comm() ); }
  const int&      dim()   const { return( space_->dim()   ); }

  Ordinal getStorageSize() const {

    return( sFields_[0]->getStorageSize()*3 );
    //    Ordinal n = 1;
    //    for(int i=0; i<3; ++i)
    //      n *= nLoc(i)+bu(i)-bl(i)+1;
    //
    //    return( 3*n );
  }

  void setStoragePtr( Scalar*  array ) {

    Ordinal n = getStorageSize()/3;

    vec_ = array;
    //    vec_[0] = array;
    //    vec_[1] = array+n;
    //    vec_[2] = array+2*n;

    for( int i=0; i<3; ++i )
      sFields_[i]->setStoragePtr( vec_+i*n );
  }

  Scalar* getRawPtr() {
    //    return( vec_[0] );
    return( vec_ );
  }

protected:

  const Ordinal& nGlo(int i)                 const { return( space_->nGlo()[i] ); }
  const Ordinal* nGlo()                 const { return( space_->nGlo() ); }

  const Ordinal& nLoc(int i)                 const { return( space_->nLoc()[i]) ; }
  const Ordinal* nLoc()                 const { return( space_->nLoc() ) ; }

  const Ordinal& sInd(int i, int fieldType)  const { return( space_->sInd(fieldType)[i] ); }
  const Ordinal& eInd(int i, int fieldType)  const { return( space_->eInd(fieldType)[i] ); }


  const Ordinal& sIndB(int i, int fieldType) const { return( space_->sIndB(fieldType)[i] ); }
  const Ordinal& eIndB(int i, int fieldType) const { return( space_->eIndB(fieldType)[i] ); }

  const Ordinal& bl(int i)                   const { return( space_->bl()[i] ); }
  const Ordinal& bu(int i)                   const { return( space_->bu()[i] ); }

  const Ordinal* bl()                   const { return( space_->bl() ); }
  const Ordinal* bu()                   const { return( space_->bu() ); }

  //  const Ordinal* sInd() const { return( space_->sInd() ); }
  //  const Ordinal* eInd() const { return( space_->eInd() ); }

  const Ordinal* sInd(  int fieldType ) const { return( space_->sInd(fieldType)  ); }
  const Ordinal* eInd(  int fieldType ) const { return( space_->eInd(fieldType) ); }

  const Ordinal* sIndB( int fieldType ) const { return( space_->sIndB(fieldType) ); }
  const Ordinal* eIndB( int fieldType ) const { return( space_->eIndB(fieldType) ); }

  const int*     bcL() const { return( space_->getDomain()->getBCLocal()->getBCL() ); }
  const int*     bcU() const { return( space_->getDomain()->getBCLocal()->getBCU() ); }

  const int* rankL() const { return( space_->getProcGrid()->getRankL() ); }
  const int* rankU() const { return( space_->getProcGrid()->getRankU() ); }

  Scalar* vec ( int i )       { return( sFields_[i]->s_ ); }
  Scalar* vecC( int i ) const { return( sFields_[i]->s_ ); }

  void changed( const int& vel_dir, const int& dir ) const {
    exchangedState_[vel_dir][dir] = false;
  }

public:

  void changed() const {
    for( int vel_dir=0; vel_dir<dim(); ++vel_dir )
      for( int dir=0; dir<dim(); ++dir )
        changed( vel_dir, dir );
  }

protected:

  bool is_exchanged( const int& vel_dir, const int& dir ) const {
    return( exchangedState_[vel_dir][dir] );
  }
  bool is_exchanged() const {
    bool all_exchanged = true;
    for( int vel_dir=0; vel_dir<dim(); ++vel_dir )
      for( int dir=0; dir<dim(); ++dir )
        all_exchanged = all_exchanged && is_exchanged(vel_dir,dir);
    return( all_exchanged );
  }

  /// \brief updates ghost layers
  void exchange( const int& vel_dir, const int& dir ) const {
    int ones[3] = {1,1,1};

    if( !exchangedState_[vel_dir][dir] ) {
      F_exchange(
          dim(),
          commf(),
          rankL(), rankU(),
          nLoc(),
          bl(), bu(),
          bcL(), bcU(),
          sInd(EField::S), eInd(EField::S),
          ones,
          nLoc(),
          dir+1, vel_dir+1,
          //          vec_[vel_dir]);
          sFields_[vel_dir]->getRawPtr() );
      exchangedState_[vel_dir][dir] = true;
    }
  }
  void exchange() const {
    for( int vel_dir=0; vel_dir<dim(); ++vel_dir )
      for( int dir=0; dir<dim(); ++dir )
        exchange( vel_dir, dir );
  }



}; // end of class VectorField



/// \brief creates a vector field belonging to a \c Space
/// \relates VectorField
template<class S=double, class O=int, int d=3>
Teuchos::RCP< VectorField<S,O,d> > createVectorField( const Teuchos::RCP< const Space<S,O,d> >& space ) {

  //  return( create< VectorField<S,O,d> >() );
  return( Teuchos::rcp(
      new VectorField<S,O,d>( space ) ) );

}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_VECTORFIELD_HPP
