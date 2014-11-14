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
class TimeField : private AbstractField<typename Field::SpaceT> {

  template<class Op,bool CNY>
  friend class TimeOpWrap;
  template<class SpaceTT>
  friend class DtTimeOp;
  template<class SpaceTT, bool CNY>
  friend class TimeNonlinearJacobian;
//  template<class S,class O>
//  friend void initVectorTimeField();

public:

  typedef typename Field::Scalar Scalar;
  typedef typename Field::Ordinal Ordinal;

  static const int dimension = Field::dimension;

  typedef typename Field::SpaceT SpaceT;

  typedef Pimpact::TimeField<Field> MV;

  typedef Scalar* ScalarArray;

  typedef AbstractField<SpaceT> AF;


  Teuchos::Array< Teuchos::RCP<Field> > mfs_;

public:

  typedef Teuchos::Array< Teuchos::RCP<Field> > FieldArray;
  typedef typename FieldArray::iterator Iter;

  Iter sInd_;
  Iter eInd_;

protected:

  ScalarArray array_;

  mutable bool exchangedState_;

public:

  TimeField( Teuchos::RCP<const SpaceT> space ):
    AF( space ),
    exchangedState_(true) {

    Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

    mfs_ = Teuchos::Array< Teuchos::RCP<Field> >( nt );

    for( int i=0; i<nt; ++i )
      mfs_[i] = Teuchos::rcp( new Field( space, false ) );

    Ordinal nx = mfs_[0]->getStorageSize();

    array_ = new Scalar[nx*nt];
    for( int i=0; i<nt; ++i )
      array_[i] = 0.;

    for( int i=0; i<nt; ++i )
      mfs_[i]->setStoragePtr( array_+i*nx );

    sInd_ = mfs_.begin() - space()->bl(3);
    eInd_ = mfs_.end()   - space()->bu(3);

  }



  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ShallowCopy is done but allows also to deepcopy the field
  TimeField(const TimeField& field, ECopyType copyType=DeepCopy):
    AF( field.space() ),
//    space()(field.space()),
    exchangedState_(field.exchangedState_) {

    Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

    mfs_ = Teuchos::Array< Teuchos::RCP<Field> >(nt);

    for( int i=0; i<nt; ++i )
//      mfs_[i] = Teuchos::rcp( new Field( space(), false ) );
      mfs_[i] = field.mfs_[i]->clone(copyType);

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

    sInd_ = mfs_.begin() - space()->bl(3);
    eInd_ = mfs_.end() - space()->bu(3);

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
  /// \param noxVec if \c TimeField is used for NOX the Vector length is
  /// considered for all Fields
  Ordinal getLength( bool noxVec=true ) const {
    return( space()->nGlo()[3]*mfs_[0]->getLength(noxVec) );
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
    Iter j = const_cast<MV&>(A).sInd_;
    Iter k = const_cast<MV&>(B).sInd_;
    for( Iter i=sInd_; i<eInd_; ++i ) {
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
    Iter j=const_cast<MV&>(y).sInd_;
    for( Iter i=sInd_; i<eInd_; ++i )
      (*i)->abs( **(j++) );
    changed();
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f]
  /// \return Reference to this object
  void reciprocal(const MV& y){
    Iter j=const_cast<MV&>(y).sInd_;
    for( Iter i=sInd_; i<eInd_; ++i )
      (*i)->reciprocal( **(j++) );
    changed();
  }


  /// \brief Scale each element of every \c Field by \c gamma.
  ///
  /// Here x represents on \c Field, and we update it as
  /// \f[ x_i = \alpha x_i \quad \mbox{for } i=1,\dots,n \f]
  void scale( const Scalar& alpha ) {
    for( Iter i=sInd_; i<eInd_; ++i )
      (*i)->scale(alpha);
    changed();
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot y_i \quad \mbox{for } i=1,\dots,n \f]
  /// \return Reference to this object
  void scale(const MV& y) {
    Iter j=const_cast<MV&>(y).sInd_;
    for( Iter i=sInd_; i<eInd_; ++i )
      (*i)->scale( **(j++) );
    changed();
  }



  //@}



  /// \brief Compute the inner product for the \c TimeField considering it as one Vector.
  Scalar dot( const MV& A, bool global=true ) const {

    Scalar b = 0.;

    Iter j = const_cast<MV&>(A).sInd_;
    for( Iter i=sInd_; i<eInd_; ++i )
      b+= (*i)->dot( **(j++), false );

    if( global ) this->reduceNorm( comm(), b );

    return( b );
  }



  /// \brief Compute the norm for the \c TimeField as it is considered as one Vector .
  /// \todo implement OneNorm
  Scalar norm(  Belos::NormType type = Belos::TwoNorm, bool global=true ) const {

    Scalar normvec = 0.;

    for( Iter i=sInd_; i<eInd_; ++i ) {
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

    Iter j = const_cast<MV&>(weights).sInd_;

    for( Iter i=sInd_; i<eInd_; ++i )
      nor+= (*i)->norm( **(j++) );

    if( global ) this->reduceNorm( comm(), nor, Belos::TwoNorm );

    return( nor );
  }


  /// \brief mv := A.
  ///
  /// assign (deep copy) A into mv.
  void assign( const MV& A ) {
    Iter j = const_cast<MV&>(A).sInd_;
    for( Iter i=sInd_; i<eInd_; ++i )
      (*i)->assign( **(j++) );
    changed();
  }


  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed = false, int seed = 1) {
    for( Iter i=sInd_; i<eInd_; ++i )
      (*i)->random();
    changed();
  }


  /// \brief \f[ *this = \alpha \f]
  void init( Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() ) {
    for( Iter i=sInd_; i<eInd_; ++i )
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
    //      (*i)->write(count++ + 2.*space()->getShift()[3] );
    for( Iter i=sInd_; i<eInd_; ++i )
      (*i)->write(count++ + space()->getShift()[3] );
  }


  MPI_Comm comm() const { return( space()->commST() ); }

  Teuchos::RCP<const SpaceT> space() const { return( AF::space_ ); }

public:

  void changed() {
    exchangedState_ = false;
  }

  /// \note shoud be constant but weirdly then Iter becomes const iter and can't be converted to int
  void exchange() const {

    if( !exchangedState_ ) {
      if( space()->getNProc(3)>=1 ) {

        int transL = std::abs( space()->bl(3) );
//        int transU = std::abs( space()->bu(3) );

        // std::cout << "transl: " <<  transl<< "\n";
        // std::cout << "transu: " <<  transu<< "\n";

        int rankU = space()->getProcGrid()->getRankU(3);
        int rankL = space()->getProcGrid()->getRankL(3);

        MPI_Request reqL;
//        MPI_Request reqU;

        MPI_Status statusL;
//        MPI_Status statusU;

        Ordinal lengthL = transL * mfs_[0]->getStorageSize();
//        Ordinal lengthU = transU * mfs_[0]->getStorageSize();

        Scalar* ghostUR = mfs_[0]->getRawPtr();
//        Scalar* ghostLR = (*(eInd_))->getRawPtr();

        Scalar* ghostUS = ( *(eInd_-transL) )->getRawPtr();
        //      Scalar* ghostLS = ;

        if( transL>0 ) MPI_Irecv( ghostUR, lengthL, MPI_REAL8, rankL, 1, comm(), &reqL);
        //      if( transU>0 ) MPI_Irecv( ghostLR, lengthU, MPI_REAL8, rankU, 2, comm(), &reqU);
        //
        if( transL>0 ) MPI_Send ( ghostUS, lengthL, MPI_REAL8, rankU, 1, comm() );
        //      if( transL>0 ) MPI_Send ( ghostLS, lengthU, MPI_REAL8, rankL, 2, comm() );

        // depends on if field from sender was exchanged, so to be sure
        mfs_[0]->changed();
        if( transL>0 ) MPI_Wait( &reqL, &statusL );
//        if( transU>0 ) MPI_Wait( &reqU, &statusU );

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
/// \todo substitue
/// \deprecated
template<class FieldT, class SpaceT>
Teuchos::RCP< TimeField<FieldT> >
createTimeField( const Teuchos::RCP<const SpaceT>& space ) {

  return( Teuchos::rcp( new TimeField<FieldT>( space ) ) );

}


#include "Pimpact_VectorField.hpp"
#include "cmath"

template<class SpaceT>
Teuchos::RCP<TimeField<VectorField<SpaceT> > >
initVectorTimeField(
    Teuchos::RCP<TimeField<VectorField<SpaceT> > > field,
    EFlowType flowType=Zero2DFlow,
    typename SpaceT::Scalar xm=0.5,
    typename SpaceT::Scalar ym=0.5,
    typename SpaceT::Scalar rad=0.1,
    typename SpaceT::Scalar amp=0.25 ) {

  typedef typename SpaceT::Scalar S;
  typedef typename SpaceT::Ordinal O;

  typedef TimeField<VectorField<SpaceT> > Field;
  typedef typename Field::Iter Iter;

  auto space = field->space();

  auto offset = space->getShift()[3];

  O i = -1;
  S pi = 4.*std::atan(1.);

//  for( Iter j = field->sInd_; j<field->eInd_; ++j )
  for( Iter j = field->mfs_.begin(); j<field->mfs_.end(); ++j )
    switch( flowType ) {
    case Zero2DFlow:
      (*j)->initField( ZeroFlow );
      break;
    case Poiseuille_inX:
      (*j)->initField( PoiseuilleFlow2D_inX );
      break;
    case Poiseuille_inY:
      (*j)->initField( PoiseuilleFlow2D_inY );
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
      (*j)->initField( ZeroFlow );
      break;
    }

  return( field );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMEFIELD_HPP
