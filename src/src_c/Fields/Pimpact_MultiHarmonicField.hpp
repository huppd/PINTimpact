/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_MULTIHARMONICFIELD_HPP
#define PIMPACT_MULTIHARMONICFIELD_HPP


#include <vector>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "BelosTypes.hpp"


#include "Pimpact_AbstractField.hpp"
#include "Pimpact_ModeField.hpp"




namespace Pimpact {


/// \brief important basic Vector class.
///
/// vector for wrapping many fields into one multiharmonic field
/// \ingroup Field
/// \todo rm std::vector<RCP<...>>
template<class IFT>
class MultiHarmonicField : private AbstractField<typename IFT::GridT> {

public:

  using GridT = typename IFT::GridT;
  using InnerFieldT = IFT;

  enum class Global : bool { Y=true, N=false };

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using ScalarArray =  ST*;

  using FieldT = MultiHarmonicField<IFT>;

  using AF = AbstractField<GridT>;

  const Owning owning_;

  const Global global_;

  IFT field0_;

  std::vector<Teuchos::RCP<ModeField<IFT> > > fields_;

  ScalarArray s_;

  mutable bool exchangedState_;


  void allocate() {

    OT n = getStorageSize();
    setStoragePtr(new ST[n]);
    std::uninitialized_fill_n(s_, n , 0.);
  }

public:

  constexpr OT getStorageSize() {

    return (global_==Global::Y)?
      ((1 + 2*grid()->nGlo(3))*field0_.getStorageSize()):
      ((1 + 2*(grid()->ei(F::U, 3) - std::max(grid()->si(F::U, 3), 1) +
                  1))*field0_.getStorageSize());
  }

  constexpr ST* getRawPtr() {
    return s_;
  }

  void setStoragePtr(ST* array) {

    s_ = array;

    OT nx = field0_.getStorageSize();

    field0_.setStoragePtr(s_);

    if(global_==Global::Y) {

      for(OT i=0; i<grid()->nGlo(3); ++i)
        fields_[i]->setStoragePtr(s_ + nx + 2*nx*i);

    } else {

      for(OT i=0; i<=grid()->ei(F::U, 3) - std::max(grid()->si(F::U, 3), 1); ++i)
        fields_[i]->setStoragePtr(s_ + nx + 2*nx*i);
    }
  }


  MultiHarmonicField(const Teuchos::RCP<const GridT>& grid,
      const Owning owning=Owning::Y):
    MultiHarmonicField(grid, static_cast<Global>(grid->np(3)==1), owning) {};

  MultiHarmonicField(const Teuchos::RCP<const GridT>& grid, const Global global,
      const Owning owning=Owning::Y):
    AF(grid),
    owning_(owning),
    global_(global),
    field0_(grid, Owning::N),
    fields_(grid->nGlo(3)),
    exchangedState_(true) {

    assert(4 == GridT::dimension);
    assert(true == grid()->getStencilWidths()->spectralT());

    if(global_==Global::Y) {
      for(OT i=0; i<grid->nGlo(3); ++i)
        fields_[i] = Teuchos::rcp(new ModeField<IFT>(grid, Owning::N));
    } else {
      for(OT i=0; i<grid()->ei(F::U, 3) - std::max(grid()->si(F::U, 3), 1) + 1; ++i)
        fields_[i] = Teuchos::rcp(new ModeField<IFT>(grid, Owning::N));
    }

    if(owning_==Owning::Y) allocate();
  };


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ECopy::Shallow is done but allows also to deepcopy the field
  MultiHarmonicField(const MultiHarmonicField& vF, const ECopy copyType=ECopy::Deep):
    AF(vF.grid()),
    owning_(vF.owning_),
    global_(static_cast<Global>(vF.grid()->np(3)==1)),
    field0_(vF.field0_, copyType),
    fields_(vF.grid()->nGlo(3)),
    exchangedState_(vF.exchangedState_) {

    if(global_==Global::Y) {
      for(OT i=1; i<=grid()->nGlo(3); ++i)
        fields_[i-1] = Teuchos::rcp(new ModeField<IFT>(vF.getField(i), copyType));
    } else {
      for(OT i=0; i<grid()->ei(F::U, 3) - std::max(grid()->si(F::U, 3), 1) + 1; ++i)
        fields_[i] = Teuchos::rcp(new ModeField<IFT>(vF.getField(i+std::max(grid()->si(F::U, 3), 1)), copyType));
    }

    if(owning_==Owning::Y) {
      allocate();
      switch(copyType) {
        case ECopy::Shallow:
          break;
        case ECopy::Deep:
          *this = vF;
          break;
      }
    }
  };


  ~MultiHarmonicField() {
    if(owning_==Owning::Y) delete[] s_;
  }

  Teuchos::RCP<FieldT> clone(const ECopy ctype=ECopy::Deep) const {
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


  /// \name Attribute methods
  /// \{

protected:

  constexpr OT index(const OT i) {
    return i - 1 + ((global_==Global::Y||0==grid()->si(F::U, 3))?
             0: (-grid()->si(F::U, 3)+1));
  };

public:

  constexpr const Global& global() {
    return global_;
  }

  IFT& get0Field() {
    return field0_;
  }

  constexpr const IFT& get0Field() {
    return field0_;
  }

  ModeField<IFT>& getField(const OT i) {
    return *fields_[index(i)];
  }

  constexpr const ModeField<IFT>& getField(const OT i) {
    return *fields_[index(i)];
  }

  IFT& getCField(const OT i) {
    return fields_[index(i)]->getCField();
  }

  constexpr const IFT& getCField(const OT i) {
    return fields_[index(i)]->getCField();
  }

  IFT& getSField(const OT i) {
    return fields_[index(i)]->getSField();
  }

  constexpr const IFT& getSField(const OT i) {
    return fields_[index(i)]->getSField();
  }

  constexpr const Teuchos::RCP<const GridT>& grid() {
    return AF::grid_;
  }

  constexpr const MPI_Comm& comm() {
    return grid()->getProcGrid()->getCommWorld();
  }


  /// \brief returns the length of Field.
  ///
  /// the vector length is with regard to the inner points
  constexpr OT getLength() {

    OT len = 0;

    //len += get0Field().getLength();
    //len += 2*grid()->nGlo(3)*get0Field().getLength();
    len = 2*(1 + grid()->nGlo(3))*get0Field().getLength();

    return len;
  }


  /// \}
  /// \name Update methods
  /// \{

  /// \brief Replace \c this with \f$\alpha a + \beta b\f$.
  /// \todo add test for consistent VectorGrids in debug mode
  void add(const ST alpha, const FieldT& a, const ST beta, const FieldT& b, const B wb=B::Y) {

    if(0==grid()->si(F::U, 3))
      field0_.add(alpha, a.get0Field(), beta, b.get0Field(), wb);

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).add(alpha, a.getField(i), beta, b.getField(i), wb);

    changed();
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void abs(const FieldT& y) {

    if(0==grid()->si(F::U, 3))
      field0_.abs(y.get0Field());

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).abs(y.getField(i));

    changed();
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1, \dots, n  \f]
  /// \return Reference to this object
  void reciprocal(const FieldT& y) {

    if(0==grid()->si(F::U, 3))
      field0_.reciprocal(y.get0Field());

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).reciprocal(y.getField(i));

    changed();
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale(const ST alpha, const B wB=B::Y) {

    if(0==grid()->si(F::U, 3))
      field0_.scale(alpha, wB);

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).scale(alpha, wB);

    changed();
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void scale(const FieldT& a) {

    if(0==grid()->si(F::U, 3))
      field0_.scale(a.get0Field());

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).scale(a.getField(i));

    changed();
  }



  /// \}
  /// \name Norm method and SP
  /// \{

  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  constexpr ST dotLoc(const FieldT& a) {

    ST b = 0.;

    if(0==grid()->si(F::U, 3))
      b += 2.*get0Field().dotLoc(a.get0Field());
    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      b += getField(i).dotLoc(a.getField(i));

    return b;
  }

  /// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
  constexpr ST dot(const FieldT& y) {

    return this->reduce(comm(), dotLoc(y));
  }

  /// \brief Compute the norm of Field.
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  constexpr ST normLoc(ENorm type=ENorm::Two) {

    ST normvec = 0.;

    if(0==grid()->si(F::U, 3))
      normvec =
        (ENorm::Inf==type)?
        std::max(get0Field().normLoc(type), normvec):
        (2.*get0Field().normLoc(type));

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      normvec =
        (ENorm::Inf==type)?
        std::max(getField(i).normLoc(type), normvec):
        (normvec+getField(i).normLoc(type));

    return normvec;
  }


  /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  constexpr ST norm(ENorm type=ENorm::Two) {

    ST normvec = this->reduce(comm(), normLoc(type),
        (ENorm::Inf==type)?MPI_MAX:MPI_SUM);

    normvec = (ENorm::Two==type||ENorm::L2==type) ?
      std::sqrt(normvec) :
      normvec;

    return normvec;
  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST normLoc(const FieldT& weights) {

    ST normvec= Teuchos::ScalarTraits<ST>::zero();

    if(0==grid()->si(F::U, 3))
      normvec += get0Field().normLoc(weights.get0Field());

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      normvec += getField(i).normLoc(weights.getField(i));

    return normvec;
  }


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST norm(const FieldT& weights) {
    return std::sqrt(this->reduce(comm(), normLoc(weights)));
  }


  /// \}
  /// \name Initialization methods
  /// \{


  /// \brief mv := a
  MultiHarmonicField& operator=(const MultiHarmonicField& a) {

    if(global_==Global::Y and a.global_==Global::Y and a.exchangedState_) {
      field0_ = a.get0Field();

      for(OT i=1; i<=grid()->nGlo(3); ++i)
        getField(i) = a.getField(i);
      exchangedState_ = a.exchangedState_;
    }
    else {

      if(0==grid()->si(F::U, 3))
        field0_ = a.get0Field();

      for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
        getField(i) = a.getField(i);

      changed();
    }

    return *this;
  }


  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed=false, const B bcYes=B::Y, int seed=1) {

    if(0==grid()->si(F::U, 3))
      field0_.random(useSeed, bcYes, seed);

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).random(useSeed, bcYes, seed);

    changed();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init(const ST alpha = Teuchos::ScalarTraits<ST>::zero(), const B wB=B::Y) {

    if(0==grid()->si(F::U, 3))
      field0_.init(alpha, wB);

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).init(alpha, wB);

    changed();
  }


  ///  \brief initializes including boundaries to zero
  void initField(Teuchos::ParameterList& para, const Add add=Add::N) {

    if(0==grid()->si(F::U, 3))
      field0_.initField(para.sublist("0 mode"), add);

    if(grid()->si(F::U, 3)<=1 && 1<=grid()->ei(F::U, 3)) {
      getCField(1).initField(para.sublist("cos mode"), add);
      getSField(1).initField(para.sublist("sin mode"), add);
    }
    changed();
  }


  void extrapolateBC(const Belos::ETrans trans=Belos::NOTRANS) {

    if(0==grid()->si(F::U, 3))
      field0_.extrapolateBC(trans);

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).extrapolateBC(trans);

    changed();
  }

  void level() const {

    if(0==grid()->si(F::U, 3))
      field0_.level();

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).level();

    changed();
  }

  /// \}

  /// Print the vector.  To be used for debugging only.
  void print(std::ostream& os=std::cout) const {

    if(0==grid()->si(F::U, 3))
      get0Field().print(os);

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).print(os);
  }


  void writeEvol(const int count=0) const {

    exchange();

    if(grid()->getProcGrid()->getIB(3)==1) {
      ST pi = 4.*std::atan(1.);
      OT nf = grid()->nGlo(3);
      OT nt = 4*nf;
      Teuchos::RCP<IFT> temp = get0Field().clone(Pimpact::ECopy::Shallow);
      for(OT i=0; i<nt;  ++i) {
        *temp = get0Field();
        for(OT j=1; j<=nf; ++j) {
          temp->add(
              1., *temp,
              std::sin(2.*pi*i*(static_cast<ST>(j))/nt), getSField(j));
          temp->add(
              std::cos(2.*pi*i*(static_cast<ST>(j))/nt), getCField(j),
              1., *temp);
          temp->write(count+i);
        }
      }
    }
  }


  void write(const int count=0, const bool restart=false) const {

    if(0==grid()->si(F::U, 3))
      get0Field().write(count, restart);

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i) {
      getCField(i).write(count+2*i-1, restart);
      getSField(i).write(count+2*i, restart);
    }
  }


  void read(const int count=0, OT nf_restart=-1) {

    if(0==grid()->si(F::U, 3))
      get0Field().read(count);

    if(-1==nf_restart)
      nf_restart = grid()->ei(F::U, 3);
    else if(nf_restart>grid()->ei(F::U, 3)) 
      nf_restart = grid()->ei(F::U, 3);

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=nf_restart; ++i) {
      getCField(i).read(count+2*i-1);
      getSField(i).read(count+2*i);
    }
    changed();
  }


  /// \name comunication methods.
  /// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
  ///
  /// \{

  void changed() const {
    exchangedState_ = false;

    if(0==grid()->si(F::U, 3))
      field0_.changed();

    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
      getField(i).changed();
  }


  void exchange() const {

    assert(global_==Global::Y);

    // check if exchange is necessary
    if(exchangedState_==false) {

      // exchange spatial
      if(0==grid()->si(F::U, 3))
        field0_.exchange();

      for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
        getField(i).exchange();

      // mpi stuff
      int nx = get0Field().getStorageSize();

      // --- sendcount ---
      int sendcount = 0;
      if(0==grid()->si(F::U, 3))
        sendcount += nx;
      for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i)
        sendcount += 2*nx;

      // --- recvcount, displacement ---
      int np = grid()->getProcGrid()->getNP(3);

      int* recvcounts = new int[np];
      int* displs     = new int[np];

      OT nfl  = (grid()->getGridSizeLocal()->get(3)+1)/np;
      OT rem  = (grid()->getGridSizeLocal()->get(3)+1)%np;

      //MPI_Barrier(grid()->getProcGrid()->getCommBar(3));

      for(int rank=0; rank<np; ++rank) {
        if(0==rank)
          recvcounts[rank] = nx + (nfl-1)*2*nx + (rank<rem?2*nx:0);
        else
          recvcounts[rank] = nfl*2*nx + (rank<rem?2*nx:0);
        if(0==rank)
          displs[rank] = 0;
        else
          displs[rank] = displs[rank-1] + recvcounts[rank-1];
//				if(0==grid()->rankST())
//					std::cout << "\trank: " << rank
//						<< "\trecevcount: " << recvcounts[rank]
//						<< "\tdispl: " << displs[rank] << "\n";
      }

      // exchange modes
      MPI_Allgatherv(
          MPI_IN_PLACE,                           // sendbuf
          sendcount,                              // sendcount
          MPI_REAL8,                              // sendtype
          s_,                                     // recvbuf
          recvcounts,                             // recvcounts
          displs,                                 // displs
          MPI_REAL8,                              // recvtype
          grid()->getProcGrid()->getCommBar(3)); // comm
         

      delete[] recvcounts;
      delete[] displs;

      // set non owning spatial block as exchanged
      if(0!=grid()->si(F::U, 3))
        get0Field().setExchanged();

      for(OT i=1; i<grid()->si(F::U, 3); ++i)
        getField(i).setExchanged();

      for(OT i=grid()->ei(F::U, 3)+1; i<=grid()->nGlo(3); ++i)
        getField(i).setExchanged();
    }
  }

  /// \}

}; // end of class MultiHarmonicField



/// \brief creates a multi-harmonic scalar field.
///
/// \relates MultiHarmonicField
/// \param grid scalar Vector Grid to which returned vector belongs
template<class FieldT>
Teuchos::RCP<MultiHarmonicField<FieldT > > createMultiHarmonic(
  const Teuchos::RCP<const typename FieldT::GridT >& grid) {

  return create<MultiHarmonicField<FieldT> >(grid);
}



/// \brief creates a multi-harmonic scalar field.
///
/// \relates MultiHarmonicField
/// \param grid scalar Vector Grid to which returned vector belongs
/// \param global
template<class FieldT>
Teuchos::RCP<MultiHarmonicField<FieldT> > createMultiHarmonic(
  const Teuchos::RCP<const typename FieldT::GridT >& grid, const typename MultiHarmonicField<FieldT>::Global global) {

  return Teuchos::rcp(new  MultiHarmonicField<FieldT>(grid, global));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIHARMONICFIELD_HPP
