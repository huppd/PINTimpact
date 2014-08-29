#pragma once
#ifndef PIMPACT_TIMEFIELD_HPP
#define PIMPACT_TIMEFIELD_HPP


#include "mpi.h"

#include <vector>
#include <Teuchos_Array.hpp>
#include "Teuchos_RCP.hpp"
#include <Teuchos_Range1D.hpp>

#include <BelosTypes.hpp>


#include "Pimpact_Space.hpp"
#include "Pimpact_Types.hpp"

#include "Pimpact_AbstractField.hpp"



namespace Pimpact {


/// \brief templated class which is the interface to \c Belos and \c NOX
///
/// has multiple \c Field's, where \c Field can be a \c Pimpact:ScalarField, \c
/// Pimpact:VectorField or some combination with \c Pimpact::ModeField or \c
/// Pimpact::CompoundField \note if this is heavily used for many Field's, then
/// the implementation should be improved such that communication is done such
/// that only done once per MV not per Field
///
/// \todo decide if it is better to modifie ghostlayer or exchange eventuell more
/// \todo maybe move functionality in Scalar/VectorField<S,O,4>
/// \ingroup Field
template<class Field>
class TimeField : private AbstractField<typename Field::Scalar,typename Field::Ordinal> {

  template<class Op,bool CNY>
  friend class TimeOpWrap;
  template<class S, class O>
  friend class DtTimeOp;
  template<class S, class O,bool CNY>
  friend class TimeNonlinearJacobian;
//  template<class S,class O>
//  friend void initVectorTimeField();

public:

  typedef typename Field::Scalar Scalar;
  typedef typename Field::Ordinal Ordinal;

  static const int dimension = 4;

  typedef Space<Scalar,Ordinal,4> SpaceT;

public:

  typedef Pimpact::TimeField<Field> MV;
  typedef Scalar* ScalarArray;

  Teuchos::RCP<const SpaceT > space_;

  Teuchos::Array< Teuchos::RCP<Field> > mfs_;

public:

  typedef Teuchos::Array< Teuchos::RCP<Field> > FieldArray;
  typedef typename FieldArray::iterator Iter;

  Iter beginI_;
  Iter endI_;

protected:

  ScalarArray array_;

  bool exchangedState_;

public:

  TimeField( Teuchos::RCP<const SpaceT > space ):
    space_(space),exchangedState_(true) {

    Ordinal nt = space_->nLoc()[3]+space_->bu()[3]-space_->bl()[3];

    mfs_ = Teuchos::Array< Teuchos::RCP<Field> >( nt );

    for( int i=0; i<nt; ++i )
      mfs_[i] = Teuchos::rcp( new Field( space_, false ) );

    Ordinal nx = mfs_[0]->getStorageSize();

    array_ = new Scalar[nx*nt];
    for( int i=0; i<nt; ++i )
      array_[i] = 0.;

    for( int i=0; i<nt; ++i )
      mfs_[i]->setStoragePtr( array_+i*nx );

    beginI_ = mfs_.begin()-space_->bl()[3];
    endI_   = mfs_.end()  -space_->bu()[3];

  }



  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
  TimeField(const TimeField& field, ECopyType copyType=DeepCopy):
    space_(field.space_),exchangedState_(field.exchangedState_) {

    Ordinal nt = space_->nLoc()[3]+space_->bu()[3]-space_->bl()[3];

    mfs_ = Teuchos::Array< Teuchos::RCP<Field> >(nt);

    for( int i=0; i<nt; ++i )
      mfs_[i] = Teuchos::rcp( new Field( space_, false ) );
//      mfs_[i] = field.mfs_[i]->clone(copyType);

    Ordinal nx = mfs_[0]->getStorageSize();

    array_ = new Scalar[nx*nt];

    for( int i=0; i<nt; ++i )
      mfs_[i]->setStoragePtr( array_+i*nx );

    if( DeepCopy==copyType )
      for( int i=0; i<nt; ++i ) {
        mfs_[i]->assign( *(field.mfs_[i]) );
      }
    else
      for( int i=0; i<nt*nx; ++i )
        array_[i] = 0.;

    beginI_ = mfs_.begin()-space_->bl()[3];
    endI_ = mfs_.end()-space_->bu()[3];

    if( ShallowCopy ) exchangedState_ = true;
  }

  ~TimeField() {
    delete[] array_;
  }



  /// \brief Create a new \c TimeField with
  Teuchos::RCP< MV > clone( ECopyType ctype = DeepCopy ) const {
    auto mv_ = Teuchos::rcp( new MV(*this,ctype) );
    return( mv_ );
  }


  /// \brief returns the length of Field.
  ///
  /// \param nox_vec if \c TimeField is used for NOX the Vector length is
  /// considered for all Fields
  Ordinal getLength( bool noxVec=true ) const {
    return( space_->nGlo()[3]*mfs_[0]->getLength(noxVec) );
  }


  /// \brief get number of stored Field's
private:

  int getNumberVecs() const {  return( mfs_.size() ); }

public:

  /// \brief is true
  bool HasConstantStride() const { return( true ); }

  //@}
  /// \name Update methods
  //@{


  /// \brief <tt>mv := alpha*A + beta*B</tt>
  ///
  ///	The Tpetra specialization of this method ignores and completely
  ///	overwrites any NaN or Inf entries in A.  Thus, it does <i>not</i> mean
  ///	the same thing as <tt>mv := 0*mv + alpha*A + beta*B</tt> in IEEE 754
  ///	floating-point arithmetic. (Remember that NaN*0 = NaN.)
  void add( Scalar alpha, const MV& A, Scalar beta, const MV& B ) {
    Iter j = const_cast<MV&>(A).beginI_;
    Iter k = const_cast<MV&>(B).beginI_;
    for( Iter i=beginI_; i<endI_; ++i ) {
      (*i)->add( alpha, **(j++), beta, **(k++) );
    }
    changed();
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void abs(const MV& y) {
    Iter j=const_cast<MV&>(y).beginI_;
    for( Iter i=beginI_; i<endI_; ++i )
      (*i)->abs( **(j++) );
    changed();
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal(const MV& y){
    Iter j=const_cast<MV&>(y).beginI_;
    for( Iter i=beginI_; i<endI_; ++i )
      (*i)->reciprocal( **(j++) );
    changed();
  }


  /// \brief Scale each element of every \c Field by \c gamma.
  ///
  /// Here x represents on \c Field, and we update it as
  /// \f[ x_i = \alpha x_i \quad \mbox{for } i=1,\dots,n \f]
  void scale( const Scalar& alpha ) {
    for( Iter i=beginI_; i<endI_; ++i )
      (*i)->scale(alpha);
    changed();
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot y_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void scale(const MV& y) {
    Iter j=const_cast<MV&>(y).beginI_;
    for( Iter i=beginI_; i<endI_; ++i )
      (*i)->scale( **(j++) );
    changed();
  }



  //@}



  /// \brief Compute the inner product for the \c TimeField considering it as one Vector.
  Scalar dot( const MV& A, bool global=true ) const {

    Scalar b = 0.;

    Iter j = const_cast<MV&>(A).beginI_;
    for( Iter i=beginI_; i<endI_; ++i )
      b+= (*i)->dot( **(j++), false );

    if( global ) this->reduceNorm( comm(), b );

    return( b );
  }



  /// \brief Compute the norm for the \c TimeField as it is considered as one Vector .
  /// \todo implement OneNorm
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = 0.;

    for( Iter i=beginI_; i<endI_; ++i ) {
      switch(type) {
      case Belos::OneNorm:
        normvec += (*i)->norm(type,false);
        break;
      case Belos::TwoNorm:
        normvec += (*i)->norm(type,false);
        break;
      case Belos::InfNorm:
        normvec = std::max( (*i)->norm(type,false), normvec ) ;
        break;
      }
    }

    if( global ) this->reduceNorm( comm(), normvec, type );

    return( normvec );

  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  double norm( const MV& weights, bool global=true ) const {

    double nor=0.;

    Iter j = const_cast<MV&>(weights).beginI_;

    for( Iter i=beginI_; i<endI_; ++i )
      nor+= (*i)->norm( **(j++) );

    if( global ) this->reduceNorm( comm(), nor, Belos::TwoNorm );

    return( nor );
  }


  /// \brief mv := A.
  ///
  /// assign (deep copy) A into mv.
  void assign( const MV& A ) {
    Iter j = const_cast<MV&>(A).beginI_;
    for( Iter i=beginI_; i<endI_; ++i )
      (*i)->assign( **(j++) );
    changed();
  }


  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed = false, int seed = 1) {
    for( Iter i=beginI_; i<endI_; ++i )
      (*i)->random();
    changed();
  }


  /// \brief \f[ *this = \alpha \f]
  void init( Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    for( Iter i=beginI_; i<endI_; ++i )
      (*i)->init(alpha);
    changed();
  }


  /// \param os
  void print( std::ostream& os ) {
    for( Iter i=mfs_.begin(); i<mfs_.end(); ++i )
      (*i)->print( os );
  }



  void write( int count=0 )  {
    //    for( Iter i=mfs_.begin(); i<mfs_.end(); ++i )
    //      (*i)->write(count++ + 2.*space_->shift()[3] );
    for( Iter i=beginI_; i<endI_; ++i )
      (*i)->write(count++ + space_->shift()[3] );
  }


  MPI_Comm comm() const { return( space_->commST() ); }

  Teuchos::RCP<const SpaceT > getSpace() const { return( space_ ); }


public:
  void changed() {
    exchangedState_ = false;
  }

  void exchange()   {

    if( !exchangedState_ ) {
      if( space_->getNProc(3)>=1 ) {
        int transL = beginI_-mfs_.begin();
        //      int transU = mfs_.end()-endI_;

        //    std::cout << "transl: " <<  transl<< "\n";
        //    std::cout << "transu: " <<  transu<< "\n";
        int rankU = space_->getProcGrid()->getRankU(3);
        int rankL = space_->getProcGrid()->getRankL(3);

        MPI_Request reqL;
        //      MPI_Request reqU;

        MPI_Status statusL;
        //      MPI_Status statusU;

        Ordinal lengthL = transL * mfs_[0]->getStorageSize();
        //      Ordinal lengthU = transU * mfs_[0]->getStorageSize();

        Scalar* ghostUR = mfs_[0]->getRawPtr();
        //      Scalar* ghostLR = (*(endI_))->getRawPtr();

        Scalar* ghostUS = ( *(endI_-transL) )->getRawPtr();
        //      Scalar* ghostLS = ;
        //
        if( transL>0 ) MPI_Irecv( ghostUR, lengthL, MPI_REAL8, rankL, 1, comm(), &reqL);
        //      if( transU>0 ) MPI_Irecv( ghostLR, lengthU, MPI_REAL8, rankU, 2, comm(), &reqU);
        //
        if( transL>0 ) MPI_Send ( ghostUS, lengthL, MPI_REAL8, rankU, 1, comm() );
        //      if( transL>0 ) MPI_Send ( ghostLS, lengthU, MPI_REAL8, rankL, 2, comm() );
        //
        mfs_[0]->changed();
        if( transL>0 ) MPI_Wait( &reqL, &statusL );

      }
      else {
        mfs_[0]->assign( **(mfs_.end()-1) );
        mfs_[0]->changed();
      }
    }

    exchangedState_ = true;
  }

  Teuchos::RCP<Field> getFieldPtr( int i ) { return(  mfs_[i] ); }
  Field& getField   ( int i ) { return( *mfs_[i] ); }

  Teuchos::RCP<const Field> getConstFieldPtr( int i ) const { return(  mfs_[i] ); }
  const Field&  getConstField   ( int i ) const { return( *mfs_[i] ); }

}; // end of class TimeField



/// \brief factory for \c TimeField
/// \relates TimeField
template<class Field>
Teuchos::RCP< TimeField<Field> >
createTimeField( const Teuchos::RCP<const Space<typename Field::Scalar,typename Field::Ordinal,4> >& space ) {

  return( Teuchos::rcp( new TimeField<Field>( space ) ) );

}


#include "Pimpact_VectorField.hpp"
#include "cmath"

template<class S,class O>
Teuchos::RCP<TimeField<VectorField<S,O,4> > >
initVectorTimeField(
    Teuchos::RCP<TimeField<VectorField<S,O,4> > > field,
    EFlowType flowType=Zero2DFlow,
    S xm=0.5,
    S ym=0.5,
    S rad=0.1,
    S amp=0.25 ) {
  typedef TimeField<VectorField<S,O,4> > Field;
  typedef typename Field::Iter Iter;

  auto space = field->getSpace();

  auto offset = space->shift()[3];

  O i = -1;
  S pi = 4.*std::atan(1.);

//  for( Iter j = field->beginI_; j<field->endI_; ++j )
  for( Iter j = field->mfs_.begin(); j<field->endI_; ++j )
    switch( flowType ) {
    case Zero2DFlow:
      (*j)->initField( ZeroProf );
      break;
    case Poiseuille_inX:
      (*j)->initField( Poiseuille2D_inX );
      break;
    case Poiseuille_inY:
      (*j)->initField( Poiseuille2D_inY );
      break;
    case Streaming2DFlow: {
      S ampt = std::sin( 2.*pi*((S)i+++(S)offset)/(S)space->nGlo()[3] );
      (*j)->initField( Streaming2D, ampt );
      break;
    }
    case OscilatingDisc2D: {
      S ymt = ym+amp*std::sin( 2.*pi*((S)i+++(S)offset)/(S)space->nGlo()[3] );
      S xmt = xm;
      (*j)->initField( Disc2D, xmt, ymt, rad );
      break;
    }
    case OscilatingDisc2DVel: {
      S yvelt = amp*std::cos( 2.*pi*((S)i+++(S)offset)/(S)space->nGlo()[3] );
      S xvelt = 0;
      (*j)->init( Teuchos::tuple( xvelt, yvelt, 0.) );
      break;
    }
    default:
      (*j)->initField( ZeroProf );
      break;
    }
//  field->changed();
//  field->exchange();

  return( field );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMEFIELD_HPP
