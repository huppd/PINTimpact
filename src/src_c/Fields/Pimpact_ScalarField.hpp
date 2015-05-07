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

//#ifdef DEBUG
//		 for(int i=0; i<N; ++i)
//			 s_[i] = 0.;
//#endif // end of #ifdef DEBUG
			initField();
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
//#ifdef DEBUG
//			 for(int i=0; i<N; ++i)
//				 s_[i] = 0;
//#endif // end of #ifdef DEBUG
				initField();
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


  ///  \brief initializes VectorField with the initial field defined in Fortran
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
          space()->getDomain()->getDomainSize()->getSize( X ),
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
          space()->getDomain()->getDomainSize()->getSize( Y ),
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
          space()->getDomain()->getDomainSize()->getSize( Z ),
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
          space()->getDomain()->getDomainSize()->getSize( X ),
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
          space()->getDomain()->getDomainSize()->getSize( Y ),
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
          space()->getDomain()->getDomainSize()->getSize( Z ),
          space()->getCoordinatesLocal()->getX( Z, fType_ ),
          s_ );
      break;
		case FPoint :
			Scalar xc[3] =
			{ 
				1.,
				1.,
//				space()->getDomain()->getDomainSize()->getSize( X )/4.,
//				space()->getDomain()->getDomainSize()->getSize( Y )/2.,
				space()->getDomain()->getDomainSize()->getSize( Z )/2. };
			Scalar amp = alpha; //2./space()->getDomain()->getDomainSize()->getRe();
			Scalar sig[3] = { 0.1, 0.1, 0.1 };
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


	void level() const {

		if( EField::S == fType_ ) {
			auto m = getLength();
			auto n = space()->nGlo();

			auto bcl =  space()->getDomain()->getBCGlobal()->getBCL();
			auto bcu =  space()->getDomain()->getBCGlobal()->getBCL();

			if( bcl[0]>0 && bcl[1]>0 ) m -= n[2]-1;
			if( bcl[0]>0 && bcu[1]>0 ) m -= n[2]-1;
			if( bcu[0]>0 && bcl[1]>0 ) m -= n[2]-1;
			if( bcu[0]>0 && bcu[1]>0 ) m -= n[2]-1;
					
			if( bcl[0]>0 && bcl[2]>0 ) m -= n[1]-1;
			if( bcl[0]>0 && bcu[2]>0 ) m -= n[1]-1;
			if( bcu[0]>0 && bcl[2]>0 ) m -= n[1]-1;
			if( bcu[0]>0 && bcu[2]>0 ) m -= n[1]-1;
					
			if( bcl[1]>0 && bcl[2]>0 ) m -= n[0]-1;
			if( bcl[1]>0 && bcu[2]>0 ) m -= n[0]-1;
			if( bcu[1]>0 && bcl[2]>0 ) m -= n[0]-1;
			if( bcu[1]>0 && bcu[2]>0 ) m -= n[0]-1;


			//		set corners to zero, such that level depends only on inner field
			SF_handle_corner(
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->getDomain()->getBCLocal()->getBCL(),
					space()->getDomain()->getBCLocal()->getBCU(),
					s_ );

			SF_level(
					MPI_Comm_c2f( space()->comm() ),
					m,
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->sIndB(fType_),
					space()->eIndB(fType_),
					s_ );

			SF_handle_corner(
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->getDomain()->getBCLocal()->getBCL(),
					space()->getDomain()->getBCLocal()->getBCU(),
					s_ );

			changed();
		}
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
				Teuchos::Tuple<Ordinal,3> N;
				for( int i=0; i<3; ++i ) {
					N[i] = space()->nGlo(i);
          if( space()->getDomain()->getBCGlobal()->getBCL(i)==Pimpact::PeriodicBC )
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
            MPI_Comm_c2f( space()->comm() ),
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
            space()->getDomain()->getDomainSize()->getRe(),
            space()->getDomain()->getDomainSize()->getAlpha2() );

      }
    }
    else {

      int stride[3] = {1,1,1};

      write_hdf_3D(
          space()->rankS(),
          MPI_Comm_c2f( space()->comm() ),
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

  const MPI_Comm& comm() const { return( space()->comm() ); }

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
          space()->getDomain()->getBCLocal()->getBCL(),
          space()->getDomain()->getBCLocal()->getBCU(),
          space()->sInd(EField::S),
          space()->eInd(EField::S),
//				 space()->sIndB(fType_), // should it work
//          space()->eIndB(fType_),
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
