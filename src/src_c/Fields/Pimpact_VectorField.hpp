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

//#include "Pimpact_FieldSpace.hpp"
//#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_Space.hpp"

#include "Pimpact_extern_ScalarField.hpp"
#include "Pimpact_extern_VectorField.hpp"

#include "Pimpact_AbstractField.hpp"


namespace Pimpact {


/// \brief important basic Vector class
/// vector for a vector field, e.g.: velocity,
/// here also happens the fortran wrapping
/// \ingroup Field
/// \todo There is one issue: update methods changes state but could be
///implemented, such that they keep the state, but then the boundary conditions
///have to be taken care of
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


  typedef Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > >  IndexSpaces;
  typedef Teuchos::Tuple<Teuchos::Tuple<bool,3>,3> State;

protected:

  typedef Scalar* ScalarArray;
  typedef VectorField<Scalar,Ordinal,dimension> VF;

  Teuchos::RCP< const Space<Scalar,Ordinal,dimension> > space_;

  Teuchos::Tuple<ScalarArray,3> vec_;

  const bool owning_;

  State exchangedState_;

public:

  VectorField():
    space_(Teuchos::null),
    vec_(0),
    owning_(true),
    exchangedState_(Teuchos::tuple(Teuchos::tuple(true,true,true),Teuchos::tuple(true,true,true),Teuchos::tuple(true,true,true)))
    {};

  VectorField( const Teuchos::RCP< const Space<Scalar,Ordinal,dimension> >& space, bool owning=true ):
    space_(space),
    owning_(owning),
    exchangedState_(Teuchos::tuple(Teuchos::tuple(true,true,true),Teuchos::tuple(true,true,true),Teuchos::tuple(true,true,true))) {

    if( owning_ ) {
      Ordinal N = 1;
      for(int i=0; i<3; ++i)
        N *= nLoc(i)+bu(i)-bl(i)+1;


      vec_[0] = new Scalar[3*N];
      vec_[1] = vec_[0]+N;
      vec_[2] = vec_[1]+N;

      for(int i=0; i<3; ++i)
        for(int j=0; j<N; ++j)
          vec_[i][j] = 0.;
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

    if( owning_ ) {
      Ordinal n = 1;
      for(int i=0; i<3; ++i)
        n *= nLoc(i)+bu(i)-bl(i);

      vec_[0] = new Scalar[3*n];
      vec_[1] = vec_[0]+n;
      vec_[2] = vec_[1]+n;

      switch( copyType ) {
      case ShallowCopy:
        for( int i=0; i<dim(); ++i )
          for( int j=0; j<n; ++j)
            vec_[i][j] = 0.;
//        changed();
        break;
      case DeepCopy:
        for( int i=0; i<dim(); ++i )
          for( int j=0; j<n; ++j)
            vec_[i][j] = vF.vec_[i][j];
        break;
      }
    }
  };


  ~VectorField() {
    if( owning_ ) delete[] vec_[0];
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
    for( int i=0; i<dim(); ++i ) {
      Ordinal vl = 1;
      for( int j=0; j<dim(); ++j) {
        if( i==j ) {
//          if( PeriodicBC==bc->getBCL(i) )
//            vl *= nGlo(j)-1-1;
//          else
            vl *= nGlo(j)-1;
        }
        else {
          if( PeriodicBC==bc->getBCL(j) )
            vl *= nGlo(j)-2+1;
          else
            vl *= nGlo(j)-2;
        }
      }
      n += vl;
    }
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
      if( vec_[i]==A.vec_[i] && vec_[i]==B.vec_[i] )
        SF_scale(
            nLoc(),
            bl(), bu(),
            sInd(i), eInd(i),
            vec_[i], alpha+beta );
      else if( vec_[i]==A.vec_[i] && vec_[i]!=B.vec_[i] )
        SF_add2(
            nLoc(), bl(), bu(),
            sInd(i), eInd(i),
            vec_[i], B.vec_[i],
            alpha, beta );
      else if( vec_[i]!=A.vec_[i] && vec_[i]==B.vec_[i] )
        SF_add2(
            nLoc(),
            bl(), bu(),
            sInd(i), eInd(i),
            vec_[i], A.vec_[i],
            beta, alpha );
      else if( vec_[i]!=A.vec_[i] && vec_[i]!=B.vec_[i] )
        SF_add(
            nLoc(),
            bl(), bu(),
            sInd(i), eInd(i),
            vec_[i], A.vec_[i], B.vec_[i],
            alpha, beta );
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
      SF_abs(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i], y.vec_[i] );
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
      SF_reciprocal(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i], y.vec_[i] );
    changed();
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale( const Scalar& alpha ) {
    for(int i=0; i<dim(); ++i)
      SF_scale(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i], alpha);
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
      SF_scale2(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i], a.vec_[i] );
    changed();
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  /// \todo add test in debuging mode for testing equality of VectorSpaces
  Scalar dot ( const VF& a, bool global=true ) const {
    Scalar b = 0.;

    for( int i=0; i<dim(); ++i )
      SF_dot(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i], a.vec_[i],
          b );

    if( global ) reduceNorm( comm(), b );

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
      case Belos::OneNorm:
        SF_comp1Norm(
            nLoc(),
            bl(), bu(),
            sInd(i), eInd(i),
            vec_[i],
            normvec );
        break;
      case Belos::TwoNorm:
        SF_comp2Norm(
            nLoc(),
            bl(), bu(),
            sInd(i), eInd(i),
            vec_[i],
            normvec );
        break;
      case Belos::InfNorm:
        SF_compInfNorm(
            nLoc(),
            bl(), bu(),
            sInd(i), eInd(i),
            vec_[i],
            normvec );
        break;
      }

    if( global ) reduceNorm( comm(), normvec, type );

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
      SF_weightedNorm(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i], weights.vec_[i],
          normvec );

     if( global ) reduceNorm( comm(), normvec, Belos::TwoNorm );

     return( normvec );

  }


  //@}
  /// \name Initialization methods
  //@{


  /// \brief mv := A
  ///
  /// Assign (deep copy) A into mv.
  void assign( const VF& a ) {
//#ifdef DEBUG
//    for(int i=0; i<3; ++i) {
//      TEST_EQUALITY( nLoc(i), a.Nloc(i) )
//			                            TEST_EQUALITY( bu(i), a.bu(i) )
//			                            TEST_EQUALITY( bl(i), a.bl(i) )
//    }
//#endif

    for( int dir=0; dir<dim(); ++dir)
      SF_assign(
          nLoc(),
          bl(), bu(),
          sIndB(dir), eIndB(dir),
          vec_[dir], a.vec_[dir] );
     changed();
//    Ordinal N = 1;
//    for(int i=0; i<3; ++i)
//      N *= nLoc(i)+bu(i)-bl(i);
//
//    for( int dir=0; dir<dim(); ++dir)
//      for(int i=0; i<N; ++i) {
//        vec_[dir][i] = a.vec_[dir][i];
//      }
//
//    for( int vel_dir=0; vel_dir<dim(); ++vel_dir )
//      for( int dir=0; dir<dim(); ++dir )
//        exchangedState_[vel_dir][dir] = a.exchangedState_[vel_dir][dir];
  }


  /// \brief Replace the vectors with a random vectors.
  ///
  /// depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is requiered
  void random(bool useSeed = false, int seed = 1) {
    for( int i=0; i<dim(); ++i )
      SF_random(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i]);
    changed();
  }

  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    for( int i=0; i<dim(); ++i )
      SF_init(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i], alpha);
    changed();
  }


  /// \brief Replace each element of the vector \c vec[i] with \c alpha[i].
  void init( const Teuchos::Tuple<Scalar,3>& alpha ) {
    for( int i=0; i<dim(); ++i )
      SF_init(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i], alpha[i]);
    changed();
  }


  ///  \brief initializes VectorField with the initial field defined in Fortran
  void initField( EFlowProfile flowType = Poiseuille2D_inX, double re=1., double om=1., double px = 1. ) {
    switch( flowType) {
    case ZeroProf :
      VF_init_Zero(
          nLoc(),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2] );
      break;
    case Poiseuille2D_inX :
      VF_init_2DPoiseuilleX(
          nLoc(),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2] );
      break;
    case Poiseuille2D_inY :
      VF_init_2DPoiseuilleY(
          //          nLoc(0), nLoc(1), nLoc(2),
          nLoc(),
          sIndB(0,0), sIndB(1,0), sIndB(2,0),
          eIndB(0,0), eIndB(1,0), eIndB(2,0),
          sIndB(0,1), sIndB(1,1), sIndB(2,1),
          eIndB(0,1), eIndB(1,1), eIndB(2,1),
          sIndB(0,2), sIndB(1,2), sIndB(2,2),
          eIndB(0,2), eIndB(1,2), eIndB(2,2),
          bl(0),   bl(1),   bl(2),
          bu(0),   bu(1),   bu(2),
          vec_[0], vec_[1], vec_[2] );
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
          vec_[0], vec_[1], vec_[2], re, om, px);
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
          vec_[0], vec_[1], vec_[2], re, om, px);
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
          vec_[0], vec_[1], vec_[2], re, om, px);
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
          vec_[0], vec_[1], vec_[2], re, om, px);
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
          vec_[0], vec_[1], vec_[2],
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
          vec_[0], vec_[1], vec_[2],
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
          vec_[0], vec_[1], vec_[2],
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
          vec_[0], vec_[1], vec_[2] );
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
          vec_[0], vec_[1], vec_[2] );
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
          vec_[0], vec_[1], vec_[2] );
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
          vec_[0], vec_[1], vec_[2] );
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
          vec_[0], vec_[1], vec_[2] );
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
          vec_[0], vec_[1], vec_[2] );
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
          vec_[0], vec_[1], vec_[2],
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
          vec_[0], vec_[1], vec_[2],
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
          vec_[0], vec_[1], vec_[2],
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
    for( int i=0; i<dim(); ++i ) {
      std::cout << "field: " << i << "\n";
      SF_print(
          nLoc(),
          bl(), bu(),
          sInd(i), eInd(i),
          vec_[i] );
    }
  }


  void write( int count=0 ) {
    exchange();
    VF_write( vec_[0], vec_[1], vec_[2], count );
  }


public:

  /// \todo add good documetnation here
  /// @return
  const MPI_Fint& commf() const { return( space_->commf() ); }
  MPI_Comm        comm()  const { return( space_->comm() ); }
  const int&      dim()   const { return( space_->dim()   ); }

  Ordinal getStorageSize() const {

    Ordinal N = 1;
    for(int i=0; i<3; ++i)
      N *= nLoc(i)+bu(i)-bl(i);

    return( 3*N );
  }

  void setStoragePtr( Scalar*  array ) {
    Ordinal N = 1;
    for(int i=0; i<3; ++i)
      N *= nLoc(i)+bu(i)-bl(i);

    vec_[0] = array;
    vec_[1] = array+N;
    vec_[2] = array+2*N;
  }
  Scalar* getStoragePtr() {
    return( vec_[0] );
  }

protected:

  const Ordinal& nGlo(int i)                 const { return( space_->nGlo()[i] ); }
  const Ordinal& nLoc(int i)                 const { return( space_->nLoc()[i]) ; }

  const Ordinal& sInd(int i, int fieldType)  const { return( space_->sInd(fieldType)[i] ); }
  const Ordinal& eInd(int i, int fieldType)  const { return( space_->eInd(fieldType)[i] ); }

  const Ordinal& sIndB(int i, int fieldType) const { return( space_->sIndB(fieldType)[i] ); }
  const Ordinal& eIndB(int i, int fieldType) const { return( space_->eIndB(fieldType)[i] ); }

  const Ordinal& bl(int i)                   const { return( space_->bl()[i] ); }
  const Ordinal& bu(int i)                   const { return( space_->bu()[i] ); }

  const Ordinal* nLoc()                      const { return( space_->nLoc() ) ; }

  const Ordinal* bl()                   const { return( space_->bl() ); }
  const Ordinal* bu()                   const { return( space_->bu() ); }

  const Ordinal* sInd() const { return( space_->sInd() ); }
  const Ordinal* eInd() const { return( space_->eInd() ); }

  const Ordinal* sInd(  int fieldType ) const { return( space_->sInd(fieldType)  ); }
  const Ordinal* eInd(  int fieldType ) const { return( space_->eInd(fieldType) ); }

  const Ordinal* sIndB( int fieldType ) const { return( space_->sIndB(fieldType) ); }
  const Ordinal* eIndB( int fieldType ) const { return( space_->eIndB(fieldType) ); }

  const int*     bcL() const { return( space_->getDomain()->getBCLocal()->getBCL() ); }
  const int*     bcU() const { return( space_->getDomain()->getBCLocal()->getBCU() ); }

  const int* rankL() const { return( space_->getProcGrid()->getRankL() ); }
  const int* rankU() const { return( space_->getProcGrid()->getRankU() ); }

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
          sInd(), eInd(),
          ones,
          nLoc(),
          dir+1, vel_dir+1,
          vec_[vel_dir]);
      exchangedState_[vel_dir][dir] = true;
    }
  }
  void exchange() const {
    for( int vel_dir=0; vel_dir<dim(); ++vel_dir )
      for( int dir=0; dir<dim(); ++dir )
        exchange( vel_dir, dir );
  }



}; // end of class VectorField


//template<class Field>
//Teuchos::RCP<Field> create( const Teuchos::RCP<const Space<typename Field::Ordinal,Field::dimension> >& space ) {
//  return( Teuchos::null );
//}
//
//
///// \brief creates a vector field belonging to a \c FieldSpace and two \c IndexSpaces
///// \relates VectorField
//template<>
//Teuchos::RCP< VectorField<double,int,3> > create< VectorField<double,int,3> >( const Teuchos::RCP< const Space<int,3> >& space ) {
//
//  return( Teuchos::RCP<VectorField<double,int,3> > (
//      new VectorField<double,int,3>( space ) ) );
//
//}
//
///// \brief creates a vector field belonging to a \c FieldSpace and two \c IndexSpaces
///// \relates VectorField
//template<>
//Teuchos::RCP< VectorField<double,int,4> > create<VectorField<double,int,4> >(
//    const Teuchos::RCP< const Space<VectorField<double,int,4>::Ordinal,VectorField<double,int,4>::dimensionina> >& space ) {
//
//  return( Teuchos::RCP<VectorField<double,int,4> > (
//      new VectorField<double,int,4>( space ) ) );
//
//}

/// \brief creates a vector field belonging to a \c FieldSpace and two \c IndexSpaces
/// \relates VectorField
template<class S=double, class O=int, int d=3>
Teuchos::RCP< VectorField<S,O,d> > createVectorField( const Teuchos::RCP< const Space<S,O,d> >& space ) {

//  return( create< VectorField<S,O,d> >() );
  return( Teuchos::rcp(
      new VectorField<S,O,d>( space ) ) );

}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_VECTORFIELD_HPP
