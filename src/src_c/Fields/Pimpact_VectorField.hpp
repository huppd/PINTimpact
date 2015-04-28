#pragma once
#ifndef PIMPACT_VECTORFIELD_HPP
#define PIMPACT_VECTORFIELD_HPP

#include <vector>
#include <iostream>
#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_extern_ScalarField.hpp"
#include "Pimpact_extern_VectorField.hpp"

#include "Pimpact_AbstractField.hpp"

#include "Pimpact_ScalarField.hpp"




namespace Pimpact {



/// \brief important basic Vector class
/// it wraps three ScalarFields.
/// \ingroup Field
/// \relates ScalarField
template<class SpaceType>
class VectorField : private AbstractField<SpaceType> {


public:

  typedef SpaceType SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

protected:

  typedef Scalar* ScalarArray;

  typedef VectorField<SpaceT> VF;
  typedef ScalarField<SpaceT> SF;

  ScalarArray vec_;

  const bool owning_;

  Teuchos::Tuple< Teuchos::RCP<SF>, 3 > sFields_;

public:

  VectorField( const Teuchos::RCP< const SpaceT >& space, bool owning=true ):
    AbstractField<SpaceT>( space ),
    owning_(owning) {

    for( int i=0; i<3; ++i )
      sFields_[i] = Teuchos::rcp( new SF( space, false, EField(i) ) );

    if( owning_ ) {

      Ordinal N = getStorageSize()/3;

      vec_ = new Scalar[3*N];

      for( int i=0; i<3; ++i )
        sFields_[i]->setStoragePtr( vec_+i*N );

//#ifdef DEBUG
//		 for( int i=0; i<3*N; ++i )
//			 vec_[i] = 0.;
//#endif // end of #ifdef DEBUG
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
      sFields_[i] = vF.getConstFieldPtr(i)->clone(); // copytype doesnot matter here, because it's not owning

    if( owning_ ) {

      Ordinal n = getStorageSize()/3;

      vec_ = new Scalar[3*n];

      for( int i=0; i<3; ++i )
        sFields_[i]->setStoragePtr( vec_+i*n );

      switch( copyType ) {
      case ShallowCopy:
//#ifdef DEBUG
//			 for( int i=0; i<3*n; ++i )
//				 vec_[i] = 0.;
//#endif // end of #ifdef DEBUG
				initField();
        break;
      case DeepCopy:
        for( int i=0; i<3*n; ++i )
          vec_[i] = vF.vec_[i];
//        for( int i=0; i<space()->dim(); ++i )
//          sFields_[i]->assign( *vF.sFields_[i] );
//        changed();
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

    auto bc = space()->getDomain()->getBCGlobal();

    Ordinal n = 0;
    for( int i=0; i<space()->dim(); ++i )
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
  void abs(const VF& y) {
    for( int i=0; i<space()->dim(); ++i )
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
  void scale(const VF& a) {
    // add test for consistent VectorSpaces in debug mode
    for(int i=0; i<space()->dim(); ++i)
      sFields_[i]->scale( *a.sFields_[i] );
    changed();
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  Scalar dot ( const VF& a, bool global=true ) const {
    Scalar b = 0.;

    for( int i=0; i<space()->dim(); ++i )
      b += sFields_[i]->dot( *a.sFields_[i], false );

    if( global ) this->reduceNorm( comm(), b );

    return( b );

  }


  ///\}
  /// @name Norm method
  ///\{

  /// \brief Compute the norm of each individual vector.
  ///
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2^2\f$, the \c i-th column of \c this.
  /// \attention the two norm is not the real two norm but its square
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = 0.;

    for( int i=0; i<space()->dim(); ++i )
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

    for( int i=0; i<space()->dim(); ++i )
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

    for( int i=0; i<space()->dim(); ++i)
      sFields_[i]->assign( *a.sFields_[i] );
    changed();
  }


  /// \brief Replace the vectors with a random vectors.
  ///
  /// depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is requiered
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


  ///  \brief initializes VectorField with the initial field defined in Fortran
  void initField( EVectorField flowType = ConstFlow, Scalar re=0., Scalar om=0., Scalar px = 0., Scalar sca=1. ) {
    switch( flowType ) {
    case ZeroFlow :
      for( int i=0; i<space()->dim(); ++i )
        sFields_[i]->initField( ConstField );
      break;
    case ConstFlow :
        sFields_[U]->initField( ConstField, re );
        sFields_[V]->initField( ConstField, om );
        sFields_[W]->initField( ConstField, px );
      break;
    case PoiseuilleFlow2D_inX :
      for( int i=0; i<space()->dim(); ++i )
        if( U==i )
          sFields_[i]->initField( Poiseuille2D_inY );
        else
          sFields_[i]->initField( ConstField );
      break;
    case PoiseuilleFlow2D_inY :
      for( int i=0; i<space()->dim(); ++i )
        if( V==i )
          sFields_[i]->initField( Poiseuille2D_inX );
        else
          sFields_[i]->initField( ConstField );
      break;
    case PoiseuilleFlow2D_inZ :
      for( int i=0; i<space()->dim(); ++i )
        if( W==i )
          sFields_[i]->initField( Poiseuille2D_inX );
        else
          sFields_[i]->initField( ConstField );
      break;
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
          space()->getDomain()->getDomainSize()->getSize(1),
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
          space()->getDomain()->getDomainSize()->getSize(0),
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
          space()->getDomain()->getDomainSize()->getSize(1),
          space()->getCoordinatesLocal()->getX(Y,EField::S),
          re, om, px,
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case Pulsatile2D_inYS:
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
          space()->getDomain()->getDomainSize()->getSize(0),
          space()->getCoordinatesLocal()->getX(X,EField::S),
          re, om, px,
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case Streaming2D:
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
          space()->getDomain()->getDomainSize()->getSize(0),
          space()->getCoordinatesLocal()->getX(X,EField::S),
          re,
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case Streaming2DC:
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
          space()->getDomain()->getDomainSize()->getSize(0),
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
          space()->getDomain()->getDomainSize()->getSize(0),
          space()->getCoordinatesLocal()->getX(X,EField::S),
          re,
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case Circle2D:
			sFields_[U]->initField( Grad2D_inY, -1. );
			sFields_[V]->initField( Grad2D_inX,  1. );
			sFields_[W]->initField();
      break;
    case Circle2D_inXZ:
			sFields_[U]->initField( Grad2D_inY, -1. );
			sFields_[V]->initField();
			sFields_[W]->initField( Grad2D_inX,  1. );
      break;
    case RankineVortex2D:
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
          space()->getDomain()->getDomainSize()->getSize(),
          space()->getCoordinatesLocal()->getX(X,EField::S),
          space()->getCoordinatesLocal()->getX(Y,EField::S),
          space()->getCoordinatesLocal()->getX(X,EField::U),
          space()->getCoordinatesLocal()->getX(Y,EField::V),
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case GaussianForcing1D:
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
          space()->getDomain()->getDomainSize()->getSize(0),
          space()->getCoordinatesLocal()->getX(X,EField::U),
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case BoundaryFilter1D:
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
          space()->getDomain()->getDomainSize()->getSize(0),
          space()->getCoordinatesLocal()->getX(X,EField::U),
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case GaussianForcing2D:
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
          space()->getDomain()->getDomainSize()->getSize(),
          space()->getCoordinatesLocal()->getX(X,EField::S),
          space()->getCoordinatesLocal()->getX(Y,EField::S),
          space()->getCoordinatesLocal()->getX(X,EField::U),
          space()->getCoordinatesLocal()->getX(Y,EField::V),
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case BoundaryFilter2D:
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
          space()->getDomain()->getDomainSize()->getSize(),
          space()->getCoordinatesLocal()->getX(X,EField::S),
          space()->getCoordinatesLocal()->getX(Y,EField::S),
          space()->getCoordinatesLocal()->getX(X,EField::U),
          space()->getCoordinatesLocal()->getX(Y,EField::V),
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case VPoint2D:
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
          space()->getDomain()->getDomainSize()->getSize(),
          space()->getCoordinatesLocal()->getX(X,EField::U),
          space()->getCoordinatesLocal()->getX(Y,EField::S),
          re,
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
      break;
    case Disc2D:
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
    case RotationDisc2D:
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
		case SweptHiemenzFlow:
			VF_init_SHBF( 
					space()->rankST(),      
					space()->getProcGrid()->getShift(0),
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
          space()->getCoordinatesGlobal()->getX(X,EField::S),
          space()->getCoordinatesGlobal()->getX(X,EField::U),
          space()->getCoordinatesLocal()->getX(Z,EField::W),
					space()->getInterpolateV2S()->getC( X ),
					space()->getDomain()->getDomainSize()->getRe(),
					0,  // nonDim  
					0., // kappa
					0., // sweep_angle_degrees,  
					0., // sweep_angle,          
					0., // angle_attack,         
          sFields_[U]->getRawPtr(),
          sFields_[V]->getRawPtr(),
          sFields_[W]->getRawPtr() );
			break;
		case Disturbance:
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
					space()->getDomain()->getBCGlobal()->getBCL( Z ),
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
    //    }
    changed();
  }


  /// dirty hack(necessary for MG)
  void level() const {}

  //@}


  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& out=std::cout )  {
    for(int i=0; i<space()->dim(); ++i) {
      sFields_[i]->print( out );
    }
  }


  void write( int count=0, bool restart=false ) const {

		if( 0==space()->rankST() ) {
				std::ofstream xfile;
				std::ostringstream ss;
				ss << std::setw( 5 ) << std::setfill( '0' ) << count;
//        std::string fname = "v_"+ss.str();
				xfile.open( "vel_"+ ss.str() +".xmf", std::ofstream::out );
				xfile<< "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.1\">\n";
				xfile << "\t<Domain>\n";
				xfile << "\t\t<Grid Name=\"3DRectMesh\" GridType=\"Uniform\">\n";
				xfile << "\t\t\t<Topology TopologyType=\"3DRectMesh\" Dimensions=\""<< space()->nGlo(2) << " " << space()->nGlo(1) << " " << space()->nGlo(0) << "\"/>\n";
				xfile << "\t\t\t<Geometry GeometryType=\"VXVYVZ\">\n";
				xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
				xfile << "\t\t\t\t\tDimensions=\""<< space()->nGlo(0) << "\"\n";
				xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
				xfile << "\t\t\t\t\tPrecision=\"8\"\n";
				xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
				xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/VectorX\n";
				xfile << "\t\t\t\t</DataItem>\n";
				xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
				xfile << "\t\t\t\t\tDimensions=\""<< space()->nGlo(1) << "\"\n";
				xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
				xfile << "\t\t\t\t\tPrecision=\"8\"\n";
				xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
				xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/VectorY\n";
				xfile << "\t\t\t\t</DataItem>\n";
				xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
				xfile << "\t\t\t\t\tDimensions=\""<< space()->nGlo(2) << "\"\n";
				xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
				xfile << "\t\t\t\t\tPrecision=\"8\"\n";
				xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
				xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/VectorZ\n";
				xfile << "\t\t\t\t</DataItem>\n";
				xfile << "\t\t\t</Geometry>\n";
				xfile << "\t\t\t<Attribute Name=\"VelX\" AttributeType=\"Scalar\" Center=\"Node\">\n";
				xfile << "\t\t\t\t<DataItem Dimensions=\""<< space()->nGlo(2) << " " << space()->nGlo(1) << " " << space()->nGlo(0) << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
				xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/velX\n";
				xfile << "\t\t\t\t</DataItem>\n";
				xfile << "\t\t\t</Attribute>\n";
				xfile << "\t\t\t<Attribute Name=\"VelY\" AttributeType=\"Scalar\" Center=\"Node\">\n";
				xfile << "\t\t\t\t<DataItem Dimensions=\""<< space()->nGlo(2) << " " << space()->nGlo(1) << " " << space()->nGlo(0) << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
				xfile << "\t\t\t\t\tvelY_"<< ss.str() << ".h5:/velY\n";
				xfile << "\t\t\t\t</DataItem>\n";
				xfile << "\t\t\t</Attribute>\n";
				xfile << "\t\t\t<Attribute Name=\"VelZ\" AttributeType=\"Scalar\" Center=\"Node\">\n";
				xfile << "\t\t\t\t<DataItem Dimensions=\""<< space()->nGlo(2) << " " << space()->nGlo(1) << " " << space()->nGlo(0) << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
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


  Ordinal getStorageSize() const {
    return( sFields_[0]->getStorageSize()*3 );
  }

  void setStoragePtr( Scalar*  array ) {
    Ordinal n = getStorageSize()/3;

    vec_ = array;

    for( int i=0; i<3; ++i )
      sFields_[i]->setStoragePtr( vec_+i*n );
  }

  Scalar* getRawPtr() {
    return( vec_ );
  }

  Scalar* getRawPtr ( int i )       { return( sFields_[i]->getRawPtr() ); }
  Scalar* vecC( int i ) const { return( sFields_[i]->getRawPtr() ); }


  Teuchos::RCP<SF> getFieldPtr( int i ) { return(  sFields_[i] ); }
  SF& getField   ( int i ) { return( *sFields_[i] ); }

  Teuchos::RCP<const SF> getConstFieldPtr( int i ) const { return(  sFields_[i] ); }
  const SF&  getConstField   ( int i ) const { return( *sFields_[i] ); }

  Teuchos::RCP<const SpaceT> space() const { return( AbstractField<SpaceT>::space_ ); }

  const MPI_Comm& comm() const { return(space()->comm()); }

protected:

  void changed( const int& vel_dir, const int& dir ) const {
    getConstFieldPtr( vel_dir )->changed( dir );
  }

public:

  void changed() const {
    for( int vel_dir=0; vel_dir<space()->dim(); ++vel_dir )
      for( int dir=0; dir<space()->dim(); ++dir ) {
        changed( vel_dir, dir );
      }
  }


  bool is_exchanged( const int& vel_dir, const int& dir ) const {
    return( getConstFieldPtr( vel_dir )->is_exchanged( dir ) );
  }

protected:

  bool is_exchanged() const {

    bool all_exchanged = true;
    for( int vel_dir=0; vel_dir<space()->dim(); ++vel_dir )
      for( int dir=0; dir<space()->dim(); ++dir )
        all_exchanged = all_exchanged && is_exchanged(vel_dir,dir);

    return( all_exchanged );

  }

public:

  /// \brief updates ghost layers
  void exchange( const int& vel_dir, const int& dir ) const {

    getConstFieldPtr(vel_dir)->exchange(dir);
  }

  void exchange() const {
    for( int vel_dir=0; vel_dir<space()->dim(); ++vel_dir )
      for( int dir=0; dir<space()->dim(); ++dir )
        exchange( vel_dir, dir );
  }



}; // end of class VectorField



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_VECTORFIELD_HPP
