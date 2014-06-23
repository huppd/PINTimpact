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




namespace Pimpact {


/// \brief important basic Vector class
/// vector for a scalar field, e.g.: pressure,
/// \note all indexing is done in Fortran
/// \ingroup Field
template<class S, class O, int dimension=3>
class ScalarField {

  template<class S1,class O1,int dimension1>
  friend class Grad;
  template<class S1,class O1,int dimension1>
  friend class Div;
  template<class S1,class O1>
  friend class DivGradOp;
  template<class S1,class O1,int dimension1>
  friend class MGVDivGradOp;

public:

  typedef S Scalar;
  typedef O Ordinal;

protected:

  typedef Scalar* ScalarArray;
  typedef ScalarField<Scalar,Ordinal,dimension> MV;
  typedef Teuchos::Tuple<bool,3> State;

  Teuchos::RCP<const Space<Ordinal,dimension> > space_;

  ScalarArray s_;

  bool owning_;

  State exchangedState_;

public:

  ScalarField():
    s_(0),
    space_(Teuchos::null),
    owning_(true),
    exchangedState_( Teuchos::tuple(true,true,true) ) {};

  ScalarField( const Teuchos::RCP<const Space<Ordinal,dimension> >& space, bool owning=true ):
    space_(space),
    owning_(owning),
    exchangedState_( Teuchos::tuple(true,true,true) ) {

    if( owning_ ) {
      Ordinal N = 1;
      for(int i=0; i<3; ++i)
        N *= nLoc(i)+bu(i)-bl(i)+1;

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
    exchangedState_( sF.exchangedState_ ) {

    if( owning_ ) {
      Ordinal N = 1;
      for(int i=0; i<3; ++i)
        N *= nLoc(i)+bu(i)-bl(i);

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

  //  Teuchos::RCP<const FieldSpace<Ordinal> > getFieldSpace() const { return( fieldSpace_ ); }
  //  Teuchos::RCP<const IndexSpace<Ordinal> > getIndexSpace() const { return( indexSpace_ ); }

  /// \brief returns the length of Field.
  Ordinal getLength( bool dummy=false ) const {
    Ordinal vl = 1;
    for(int i = 0; i<dim(); ++i)
      vl *= nGlo(i);
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
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
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
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
        s_, y.s_ );
    changed();
  }


  /// \brief Scale each element of the vector with \c alpha.
  void scale( const Scalar& alpha ) {
    SF_scale(
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
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
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
        s_, a.s_ );
    changed();
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  /// \todo add test in debuging mode for testing equality of VectorSpaces
  Scalar dot ( const MV& a, bool global=true ) const {
    Scalar b;
    SF_dot(
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
        s_, a.s_, b);
    if( global ) {
      Scalar b_global;
      MPI_Allreduce( &b, &b_global, 1, MPI_REAL8, MPI_SUM, comm() );
      b = b_global;
    }
    return( b );
  }


  ///@}
  /// @name Norm method
  ///@{


  /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  /// \todo implement OneNorm
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {
    bool twoNorm_yes = false;
    bool infNorm_yes = false;

    switch(type) {
    case Belos::TwoNorm: twoNorm_yes = true; break;
    case Belos::InfNorm: infNorm_yes = true; break;
    case Belos::OneNorm: std::cout << "!!! Warning Belos::OneNorm not implemented \n"; return(0.);
    default: std::cout << "!!! Warning unknown Belos::NormType:\t" << type << "\n"; return(0.);
    }

    Scalar normvec;
    SF_compNorm(
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
        s_,
        infNorm_yes, twoNorm_yes,
        normvec, normvec );

    if( type==Belos::TwoNorm ) {
      if( global ) {
        Scalar normvec_global;
        MPI_Allreduce( &normvec, &normvec_global, 1, MPI_REAL8, MPI_SUM, comm() );
        normvec = normvec_global;
      }
      return( std::sqrt(normvec) );
    }
    else {
      if( global ) {
        Scalar normvec_global;
        MPI_Allreduce( &normvec, &normvec_global, 1, MPI_REAL8, MPI_MAX, comm() );
        normvec = normvec_global;
      }
      return( normvec );
    }
  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  double norm(const MV& weights) const {
    Scalar normvec;
    SF_weightedNorm(
        commf(),
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
        s_, weights.s_,
        normvec );
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
#ifdef DEBUG
    for(int i=0; i<3; ++i) {
      TEST_EQUALITY( nLoc(i), a.Nloc(i) )
			                                TEST_EQUALITY( bu(i), a.bu(i) )
			                                TEST_EQUALITY( bl(i), a.bl(i) )
    }
#endif

    Ordinal N = 1;
    for(int i=0; i<3; ++i)
      N *= nLoc(i)+bu(i)-bl(i);

    for(int i=0; i<N; ++i)
      s_[i] = a.s_[i];

    for( int dir=0; dir<dim(); ++dir )
      exchangedState_[dir] = a.exchangedState_[dir];
  }


  /// \brief Replace the vectors with a random vectors.
  /// depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is requiered
  void random( bool useSeed = false, int seed = 1 ) {
    SF_random(
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
        s_);
    changed();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init( const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    SF_init(
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
        s_, alpha);
    changed();
  }


  //@}

  /// Print the vector.  To be used for debugging only.
  void print( std::ostream& os )  const {
    int rank;
    MPI_Comm_rank(comm(),&rank);
    for(int i=0; i<3; ++i) {
      os << "rank: " << rank << " :dir: " << i << "\n";
      os << "rank: " << rank << " :nGlo: " << nGlo(i) << "\n";
      os << "rank: " << rank << " :nLoc: " << nLoc(i) << "\n";
      os << "rank: " << rank << " :sInd: " << sInd(i) << "\n";
      os << "rank: " << rank << " :eInd: " << eInd(i) << "\n";
      os << "rank: " << rank << " :bl: " << bl(i) << "\n";
      os << "rank: " << rank << " :bu: " << bu(i) << "\n\n";
    }
    //		Ordinal N = 1;
    //		for(int i=0; i<3; ++i)
    //			N *= nLoc(i)+bu(i)-bl(i);
    //		Ordinal Nx = nLoc(0)+bu(0)-bl(0);
    //		Ordinal Ny = nLoc(1)+bu(1)-bl(1);
    //		Ordinal Nz = nLoc(2)+bu(2)-bl(2);

    //		Scalar bla = s_[0];
    //					std::cout << "s_[0]: "<< bla <<"\n" ;
    //		for(int ix=0; ix<Nx; ++ix) {
    //			for(int iy=0; iy<Ny; ++iy) {
    //				for(int iz=0; iz<Nz;++iz) {
    //					std::cout << "rank: " << rank << " " <<
    //							"ind: (" << ix+bl(0) << ", " << iy+bl(1) << ", " << iz+bl(2) <<  ") u(ind): " << s_[iz + Nz*iy + Nz*Ny*ix] << "\n" ;
    //				}
    //				std::cout << "\n";
    //			}
    //			std::cout << "\n";
    //		}
    std::cout << "rank: " << rank << "\n";
    SF_print(
        nLoc(0), nLoc(1), nLoc(2),
        sInd(0), sInd(1), sInd(2),
        eInd(0), eInd(1), eInd(2),
        bl(0),   bl(1),   bl(2),
        bu(0),   bu(1),   bu(2),
        s_ );

  }


  void write( int count=0 ) {
    // exchange?
    SF_write( s_, count );
  }



public:

  const MPI_Fint& commf() const { return( space_->commf() ); }
  MPI_Comm        comm()  const { return( space_->comm() ); }
  const int&      dim()   const { return( space_->dim()   ); }

  Ordinal getStorageSize() const {

    Ordinal N = 1;
    for(int i=0; i<3; ++i)
      N *= nLoc(i)+bu(i)-bl(i)+1;

    return( N );
  }
  void setStoragePtr( Scalar*  array ) {
    Ordinal N = 1;
    for(int i=0; i<3; ++i)
      N *= nLoc(i)+bu(i)-bl(i)+1;

    s_ = array;
  }
  Scalar* getStoragePtr() {
    return( s_ );
  }

protected:

  const Ordinal& nGlo(int i)                 const { return( space_->nGlo()[i] ); }
  const Ordinal& nLoc(int i)                 const { return( space_->nLoc()[i]) ; }

  const Ordinal& sInd(int i)  const { return( space_->sInd()[i] ); }
  const Ordinal& eInd(int i)  const { return( space_->eInd()[i] ); }

  const Ordinal& bl(int i)                   const { return( space_->bl()[i] ); }
  const Ordinal& bu(int i)                   const { return( space_->bu()[i] ); }

  const Ordinal* nLoc()                      const { return( space_->nLoc() ) ; }

  const Ordinal* bl()                   const { return( space_->bl() ); }
  const Ordinal* bu()                   const { return( space_->bu() ); }

  const Ordinal* sInd() const { return( space_->sInd() ); }
  const Ordinal* eInd() const { return( space_->eInd() ); }

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
    if( !exchangedState_[dir] ) {
      F_exchange(
          commf(),
          dir+1, 0,
          1, 1, 1,
          nLoc(0), nLoc(1), nLoc(2),
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
    const Teuchos::RCP<const Space<O,d> >& space=Teuchos::null ) {
//  if( space.is_null() )
//    return( Teuchos::rcp(
//        new ScalarField<S,O,d>( createSpace<O,d>() ) ) );
//  else
    return( Teuchos::rcp(
        new ScalarField<S,O,d>( space ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SCALARFIELD_HPP
