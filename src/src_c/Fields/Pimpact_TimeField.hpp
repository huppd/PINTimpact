#pragma once
#ifndef PIMPACT_TIMEFIELD_HPP
#define PIMPACT_TIMEFIELD_HPP


#include <cmath>
#include <vector>

#include "mpi.h"

#include <Teuchos_Array.hpp>
#include "Teuchos_RCP.hpp"
#include <Teuchos_Range1D.hpp>

#include <BelosTypes.hpp>

#include "Pimpact_AbstractField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \brief templated class which is the interface to \c Belos and \c NOX
///
/// has multiple \c Field's, where \c Field can be a \c Pimpact:ScalarField, \c
/// Pimpact:VectorField or some combination with \c Pimpact::ModeField or \c
/// Pimpact::CompoundField \note if this is heavily used for many Field's, then
/// the implementation should be improved such that communication is done such
/// that only done once per FieldT not per Field
///
/// \todo rm RCP from vector
/// \ingroup Field
template<class Field>
class TimeField : private AbstractField<typename Field::GridT> {

  template<class Op, bool CNY>
  friend class TimeOpWrap;
  template<class GridTT>
  friend class DtTimeOp;

public:

  using GridT = typename Field::GridT;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using ScalarArray = ST*;

  using AF = AbstractField<GridT>;

  using FieldT = Pimpact::TimeField<Field>;

  const Owning owning_;

  Teuchos::Array<Teuchos::RCP<Field> > mfs_;

  ScalarArray s_;

  mutable bool exchangedState_;

  void allocate() {

    OT n = getStorageSize();
    setStoragePtr(new ST[n]);
    std::uninitialized_fill_n(s_, n , 0.);
  }

public:

  constexpr OT getStorageSize() {

    OT nt = grid()->nLoc(3) + grid()->bu(3) - grid()->bl(3);
    OT nx = at(0).getStorageSize();
    return nx*nt;
  }

  void setStoragePtr(ST* array) {
    OT nt = grid()->nLoc(3) + grid()->bu(3) - grid()->bl(3);
    OT nx = at(0).getStorageSize();
    for(int i=0; i<nt; ++i)
      at(i).setStoragePtr(s_+i*nx);
  }

  TimeField(Teuchos::RCP<const GridT> grid, const Owning owning=Owning::Y):
    AF(grid),
    owning_(owning),
    exchangedState_(true) {

    OT nt = grid()->nLoc(3) + grid()->bu(3) - grid()->bl(3);

    mfs_ = Teuchos::Array<Teuchos::RCP<Field> >(nt);

    for(int i=0; i<nt; ++i)
      mfs_[i] = Teuchos::rcp(new Field(grid, Owning::N));

    if(owning_==Owning::Y) allocate();
  }



  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param field
  /// \param copyType by default a ECopy::Shallow is done but allows also to deepcopy the field
  TimeField(const TimeField& field, const ECopy copyType=ECopy::Deep):
    AF(field.grid()),
    owning_(field.owning_),
    exchangedState_(field.exchangedState_) {

    OT nt = grid()->nLoc(3) + grid()->bu(3) - grid()->bl(3);

    mfs_ = Teuchos::Array<Teuchos::RCP<Field> >(nt);

    for(int i=0; i<nt; ++i)
      mfs_[i] = Teuchos::rcp(new Field(grid(), Owning::N));

    if(owning_==Owning::Y) {
      allocate();
      switch(copyType) {
        case ECopy::Shallow:
          break;
        case ECopy::Deep:
          *this = field;
          break;
      }
    }
  }


  ~TimeField() {
    if(owning_==Owning::Y) delete[] s_;
  }


  /// \brief Create a new \c TimeField with
  Teuchos::RCP<FieldT> clone(const ECopy ctype = ECopy::Deep) const {

    Teuchos::RCP<FieldT> mv = Teuchos::rcp(new FieldT(grid()));

    switch(ctype) {
      case ECopy::Shallow:
        break;
      case ECopy::Deep:
        *mv = *this;
        break;
    }

    return mv;
  }


  /// \brief returns the length of Field.
  constexpr OT getLength() {
    return grid()->nGlo(3)*at(0).getLength();
  }

public:

  /// \}
  /// \name Update methods
  /// \{


  /// \brief <tt>mv := alpha*a + beta*b</tt>
  void add(ST alpha, const FieldT& a, ST beta, const FieldT& b, const B wb=B::Y) {

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).add(alpha, a(i), beta, b(i), wb);
    changed();
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void abs(const FieldT& y) {

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).abs(y(i));
    changed();
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1, \dots, n  \f]
  /// \return Reference to this object
  void reciprocal(const FieldT& y) {

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).reciprocal(y(i));
    changed();
  }


  /// \brief Scale each element of every \c Field by \c gamma.
  ///
  /// Here x represents on \c Field, and we update it as
  /// \f[ x_i = \alpha x_i \quad \mbox{for } i=1, \dots, n \f]
  void scale(const ST alpha, const B wB=B::Y) {

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).scale(alpha, wB);
    changed();
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot y_i \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void scale(const FieldT& y) {
    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).scale(y(i));
    changed();
  }

  /// \}


  /// \brief Compute the inner product for the \c TimeField considering it as one Vector.
  constexpr ST dotLoc(const FieldT& A) const {

    ST b = 0.;

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      b+= at(i).dotLoc(A(i));

    return b;
  }

  /// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
  constexpr ST dot(const FieldT& y) const {

    return this->reduce(comm(), dotLoc(y));
  }


  /// \brief Compute the norm for the \c TimeField as it is considered as one Vector .
  constexpr ST normLoc(const ENorm type=ENorm::Two, const B bcYes=B::Y) const {

    ST normvec = 0.;

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      normvec =
        ((ENorm::Inf==type)?
          std::max(at(i).normLoc(type, bcYes), normvec) :
          (normvec + at(i).normLoc(type, bcYes)));

    return normvec;
  }

/// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm/
  constexpr ST norm(const ENorm type=ENorm::Two, const B bcYes=B::Y) const {

    ST normvec = this->reduce(comm(), normLoc(type, bcYes),
        (ENorm::Inf==type)?MPI_MAX:MPI_SUM);

    normvec = ((ENorm::Two==type||ENorm::L2==type) ?
                std::sqrt(normvec) :
                normvec);

    return normvec;
  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST normLoc(const FieldT& weights) const {

    ST nor = Teuchos::ScalarTraits<ST>::zero();

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      nor+= at(i).norm(weights(i));

    return nor;
  }

  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST norm(const FieldT& weights) const {
    return std::sqrt(this->reduce(comm(), normLoc(weights)));
  }


  /// \brief *this := a.
  ///
  /// assign (deep copy) A into mv.
  TimeField& operator=(const TimeField& a) {

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i) = a(i);
    changed();

    return *this;
  }


  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed=false, const B bcYes=B::Y, int seed=1) {

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).random(useSeed, bcYes, seed);
    changed();
  }


  /// \brief \f[ *this = \alpha \f]
  void init(const ST alpha = Teuchos::ScalarTraits<ST>::zero(), const B wB=B::Y) {

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).init(alpha, wB);
    changed();
  }

  void extrapolateBC(const Belos::ETrans trans=Belos::NOTRANS) {

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).extrapolateBC(trans);
    changed();
  }

  void level() const {

    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).level();
    changed();
  }


  /// \param os
  void print(std::ostream& os=std::cout) const {
    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).print(os);
  }



  void write(int count=0, const bool restart=false) const {
    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).write(count++ + grid()->getShift(3), restart);
  }

  void read(int count=0) {
    for(OT i=grid()->si(F::S, 3); i<=grid()->ei(F::S, 3); ++i)
      at(i).read(count++ + grid()->getShift(3));
  }

  constexpr const MPI_Comm& comm() const {
    return grid()->getProcGrid()->getCommWorld();
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return AF::grid_;
  }

public:

  void changed() const {
    exchangedState_ = false;
  }

  /// \note shoud be constant but weirdly then Iter becomes const iter and can't be converted to int
  void exchange() const {

    if(!exchangedState_) {
      if(grid()->np(3)>1) {

        int transL = std::abs(grid()->bl(3));
        int transU = std::abs(grid()->bu(3));

        // std::cout <<"transL: " << transL <<"\n";
        // std::cout <<"transU: " << transU <<"\n";

        int rankU = grid()->getProcGrid()->getRankU(3);
        int rankL = grid()->getProcGrid()->getRankL(3);

        MPI_Request reqL;
        MPI_Request reqU;

        MPI_Status statusL;
        MPI_Status statusU;

        OT lengthL = transL * at(0).getStorageSize();
        OT lengthU = transU * at(0).getStorageSize();

        ST* ghostUR = at(grid()->si(F::S, 3)-transU).getRawPtr();
        ST* ghostLR = at(grid()->ei(F::S, 3)  +transL).getRawPtr();

        ST* ghostUS = at(grid()->ei(F::S, 3)).getRawPtr();
        ST* ghostLS = at(grid()->si(F::S, 3)).getRawPtr();

        if(transL>0) MPI_Irecv(ghostUR, lengthL, MPI_REAL8, rankL, 1, comm(), &reqL);
        if(transU>0) MPI_Irecv(ghostLR, lengthU, MPI_REAL8, rankU, 2, comm(), &reqU);

        if(transL>0) MPI_Send (ghostUS, lengthL, MPI_REAL8, rankU, 1, comm());
        if(transU>0) MPI_Send (ghostLS, lengthU, MPI_REAL8, rankL, 2, comm());

        if(transL>0) MPI_Wait(&reqL, &statusL);
        if(transU>0) MPI_Wait(&reqU, &statusU);

        // depends on if field from sender was exchanged, so to be sure
        at(0).changed();
        at(grid()->ei(F::S, 3)).changed();

      } else {
        if(std::abs(grid()->bl(3))>0) {
          *mfs_[ grid()->si(F::S, 3)-1 ] = at(grid()->ei(F::S, 3));
          at(grid()->si(F::S, 3)-1).changed();
        }
        if(std::abs(grid()->bu(3))>0) {
          *mfs_[ grid()->ei(F::S, 3)+1 ] = at(grid()->si(F::S, 3));
          at(grid()->ei(F::S, 3)+1).changed();
        }
      }
    }
    exchangedState_ = true;
  }


protected:
  Field& at(const int i) {
    return *mfs_[i];
  }
  constexpr const Field& at(const int i) {
    return *mfs_[i];
  }

public:

  Field& operator()(const int i) {
    return at(i);
  }
  constexpr const Field& operator()(const int i) {
    return at(i);
  }

  constexpr ScalarArray getRawPtr() {
    return s_;
  }

  constexpr const ST* getConstRawPtr() const {
    return s_;
  }

}; // end of class TimeField



/// \brief factory for \c TimeField
/// \relates TimeField
/// \deprecated
template<class FieldT, class GridT>
Teuchos::RCP<TimeField<FieldT> >
createTimeField(const Teuchos::RCP<const GridT>& grid) {

  return Teuchos::rcp(new TimeField<FieldT>(grid));
}



/// \deprecated
template<class GridT>
void
initVectorTimeField(
  TimeField<VectorField<GridT> >& field,
  EFlowType flowType=Zero2DFlow,
  typename GridT::Scalar xm=0.5,
  typename GridT::Scalar ym=0.5,
  typename GridT::Scalar rad=0.1,
  typename GridT::Scalar amp=0.25) {

  using S = typename GridT::Scalar;
  //using O = typename GridT::Ordinal;

  Teuchos::RCP<const GridT> grid = field.grid();

  S pi = 4.*std::atan(1.);

  S nt = grid->nGlo(3);
  S offset = grid->getShift(3) - grid->si(F::S, 3);

  bool notImplemented = true;
  TEUCHOS_TEST_FOR_EXCEPT(notImplemented);
  //for(O i=grid->si(F::S, 3); i<=grid->ei(F::S, 3); ++i)
  //switch(flowType) {
  //case Zero2DFlow:
  //field(i)->initField(ZeroFlow);
  //break;
  //case Const2DFlow:
  //field(i)->initField(ConstFlow, xm, ym, rad);
  //break;
  //case Poiseuille_inX:
  //field(i)->initField(PoiseuilleFlow2D_inX);
  //break;
  //case Poiseuille_inY:
  //field(i)->initField(PoiseuilleFlow2D_inY);
  //break;
  //case Streaming2DFlow: {
  //S ampt = std::sin(2.*pi*((F::S)i+offset)/nt);
  //field(i)->initField(Streaming2D, ampt);
  //break;
  //}
  //case OscilatingDisc2D: {
  ////			std::cout <<"\ti: " <<i <<"\tt: " <<2.*pi*((F::S)i+offset)/nt <<"\tt: " <<grid->getCoordinatesLocal()->getX(F::S, ECoord::T)[i] <<"\n";
  //S ymt = ym+amp*std::sin(grid->getCoordinatesLocal()->getX(F::S, ECoord::T)[i]);
  //S xmt = xm;
  //field(i)->initField(Disc2D, xmt, ymt, rad);
  //break;
  //}
  //case OscilatingDisc2DVel: {
  //S yvelt = amp*std::cos(2.*pi*((F::S)i+offset)/nt);
  //S xvelt = 0;
  //field(i)->init(Teuchos::tuple(xvelt, yvelt, 0.));
  //break;
  //}
  //case ConstVel_inX:{
  //field(i)->init(Teuchos::tuple(-2*xm*std::cos(grid->getCoordinatesLocal()->getX(F::S, ECoord::T)[i]), 0., 0.)); // here xm = p
  //break;
  //}
  //case Pulsatile_inX: {
  ////field(i)->initField(Pulsatile2D_inX, xm, grid->getCoordinatesLocal()->getX(F::S, ECoord::T)[i], ym, rad); // the arguments are (xmt, i, ymt, rad) --> (re, t, px, alpha)
  //break;
  //}
  //default:
  //field(i)->initField(ZeroFlow);
  //break;
  //}

  //field.changed();
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TIMEFIELD_HPP
