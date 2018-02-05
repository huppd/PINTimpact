#pragma once
#ifndef PIMPACT_MULTIFIELD_HPP
#define PIMPACT_MULTIFIELD_HPP


#include <vector>

#include "Teuchos_Array.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_AbstractField.hpp"
#include "Pimpact_Grid.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief templated class which is the interface to \c Belos and \c NOX
///
/// has multiple \c Field's, where \c Field can be a \c Pimpact:ScalarField, \c
/// Pimpact:VectorField or some combination with \c Pimpact::ModeField or \c
/// Pimpact::CompoundField
///
/// \note for better documentation, look at the equivalent documentation in the \c Belos::...
/// \note continous memory is not necessarily desirable because also views are used
/// \ingroup Field
template<class IFT>
class MultiField : private AbstractField<typename IFT::GridT> {

public:

  using InnerFieldT = IFT;

  using GridT = typename InnerFieldT::GridT;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using FieldT = Pimpact::MultiField<InnerFieldT>;

  using AF = AbstractField<GridT>;

  const Owning owning_;

  Teuchos::Array<Teuchos::RCP<InnerFieldT> > mfs_;

  std::vector<ST*> s_;

  void allocate() {

    OT n = getStorageSize();
    setStoragePtr(new ST[n]);
    std::uninitialized_fill_n(s_[0], n , 0.);
  }

public:

  constexpr OT getStorageSize() {

    OT n = getNumberVecs();
    OT nx = mfs_[0]->getStorageSize();
    return nx*n;
  }

  void setStoragePtr(ST* array) {

    OT n = getNumberVecs();
    OT nx = mfs_[0]->getStorageSize();
    for(int i=0; i<n; ++i) {
      s_[i] = array + i*nx;
      mfs_[i]->setStoragePtr(s_[i]);
    }
  }


  /// \brief cheap constructor from one InnerFieldT.
  ///
  /// creates simple wrapper from field(no coppying).
  MultiField(const Teuchos::RCP<InnerFieldT>& field):
    AF(field->grid()),
    owning_(Owning::N),
    mfs_(1, Teuchos::null),
    s_(1) {

      mfs_[0] = field;
      s_[0] = field->getRawPtr();
    }


  /// \note dangerous to use, if not Owning s_ has to be set
  MultiField(const FieldT& mv, const ECopy ctype):
    AF(mv.grid()),
    owning_(mv.owning_),
    mfs_(mv.getNumberVecs()),
    s_(mv.getNumberVecs()){

    for(int i=0; i<getNumberVecs(); ++i)
      mfs_[i] = Teuchos::rcp(new InnerFieldT(grid(), Owning::N));

    if(owning_==Owning::Y) {
      allocate();
      switch(ctype) {
        case ECopy::Shallow:
          break;
        case ECopy::Deep:
          *this = mv;
          break;
      }
    }
  }

  /// \brief  constructor, creates \c numvecs  empty fields
  ///
  /// \param grid
  /// \param numvecs
  /// \return
  /// \note danger if Owning::N, the s_ is no set 
  MultiField(const Teuchos::RCP<const GridT>& grid, const Owning owning):
    AF(grid),
    owning_(owning),
    mfs_(1),
    s_(1) {

    for(int i=0; i<1; ++i)
      mfs_[i] = Teuchos::rcp(new InnerFieldT(grid, Owning::N));

    if(owning_==Owning::Y)
      allocate();
  }

  /// \brief  constructor, creates \c numvecs  empty fields
  ///
  /// \param grid
  /// \param numvecs
  /// \return
  /// \note danger if Owning::N, the s_ is no set 
  MultiField(const Teuchos::RCP<const GridT>& grid, const int numvecs=1,
      const Owning owning=Owning::Y):
    AF(grid),
    owning_(owning),
    mfs_(numvecs),
    s_(numvecs) {

    for(int i=0; i<numvecs; ++i)
      mfs_[i] = Teuchos::rcp(new InnerFieldT(grid, Owning::N));

    if(owning_==Owning::Y)
      allocate();
  }


  ~MultiField() {
    if(owning_==Owning::Y) delete[] s_[0];
  }

  /// \brief Create a new \c MultiField with \c numvecs columns.
  ///
  /// The returned Pimpact::MultiField has the same (distribution over one or
  /// more parallel processes) Its entries are not initialized and have undefined
  /// values.
  Teuchos::RCP<MultiField> clone(const int numvecs) const {

    Teuchos::RCP<MultiField> mv = Teuchos::rcp(new MultiField(grid(), numvecs));

    return mv;
  }


  /// \brief Create a new \c MultiField with \c numvecs columns.
  ///
  /// The returned Pimpact::MultiField has the same (distribution over one or
  /// more parallel processes) Its entries are not initialized and have undefined
  /// values.
  Teuchos::RCP<MultiField> clone(const ECopy ctype = ECopy::Deep) const {

    Teuchos::RCP<MultiField> mv = Teuchos::rcp(new MultiField(grid(), getNumberVecs()));

    switch(ctype) {
    case ECopy::Shallow:
      break;
    case ECopy::Deep:
      *mv = *this;
      break;
    }

    return mv;
  }


  /// \return deep copy of \c MultiField
  /// \note not needed
  Teuchos::RCP<MultiField> CloneCopy() const {
    return clone();
  }


  /// \brief deep copy of \c index fields, the new \c MultiFields stores \c index.size() vector.
  ///
  /// \param index
  /// \return
  Teuchos::RCP<MultiField> CloneCopy(const std::vector<int>& index) const {

    Teuchos::RCP<MultiField> mv_ = Teuchos::rcp(new MultiField(grid(), index.size()));

    for(unsigned int i=0; i<index.size(); ++i) {
      (*mv_->mfs_[i]) = (*mfs_[index[i]]);
    }
    return mv_;
  }


  /// deep copy of \c index fields, the new \c MultiFields stores \c index.size()
  /// \param index here index means an interval
  /// \return
  Teuchos::RCP<MultiField> CloneCopy(const Teuchos::Range1D& range) const {

    Teuchos::RCP<MultiField> mv_ = Teuchos::rcp(new MultiField(grid(), range.size()));

    int j = 0;
    for(int i=range.lbound(); i<=range.ubound(); ++i) {
      (*mv_->mfs_[j]) = (*mfs_[i]);
      ++j;
    }

    return mv_;
  }


  /// \param index
  /// \return nonConst View
  Teuchos::RCP<MultiField> CloneViewNonConst(const std::vector<int>& index) {

    Teuchos::RCP<MultiField> mv_ =
      Teuchos::rcp(new MultiField(grid(), index.size(), Owning::N));

    for(unsigned int i=0; i<index.size(); ++i) {
      mv_->mfs_[i] = mfs_[index[i]];
      mv_->s_[i] = s_[index[i]];
    }

    return mv_;
  }


  Teuchos::RCP<MultiField> CloneViewNonConst(const Teuchos::Range1D& range) {

    Teuchos::RCP<MultiField> mv_ =
      Teuchos::rcp(new MultiField(grid(), range.size(), Owning::N));

    int j=0;
    for(int i=range.lbound(); i<=range.ubound(); ++i) {
      mv_->mfs_[j] = mfs_[i];
      mv_->s_[j] = s_[i];
      ++j;
    }

    return mv_;
  }


  Teuchos::RCP<const MultiField> CloneView(const std::vector<int>& index) const {

    Teuchos::RCP<MultiField> mv_ =
      Teuchos::rcp(new MultiField(grid(), index.size(), Owning::N));

    for(unsigned int i=0; i<index.size(); ++i) {
      mv_->mfs_[i] = mfs_[index[i]];
      mv_->s_[i] = s_[index[i]];
    }

    return mv_;
  }


  Teuchos::RCP<const MultiField> CloneView(const Teuchos::Range1D& range) const {

    Teuchos::RCP<MultiField> mv_ =
      Teuchos::rcp(new MultiField(grid(), range.size(), Owning::N));

    int j=0;
    for(int i=range.lbound(); i<=range.ubound(); ++i) {
      mv_->mfs_[j] = mfs_[i];
      mv_->s_[j] = s_[i];
      ++j;
    }

    return mv_;
  }


  /// \brief returns the length of MultiField.
  constexpr OT getLength() {
    return mfs_[0]->getLength();
  }


  /// \brief get number of stored Field's
  constexpr int getNumberVecs() {
    return mfs_.size();
  }


  /// \brief is true
  constexpr bool HasConstantStride() {
    return true;
  }


  /// \}
  /// \name Update methods
  /// \{

  /// \brief addes new field at end
  void push_back(const Teuchos::RCP<InnerFieldT>& field=Teuchos::null) {
    if(Teuchos::is_null(field)) {
      mfs_.push_back(mfs_.back()->clone(ECopy::Shallow));
      s_.push_back(mfs_.back()->getRawPtr());
    }
    else {
      mfs_.push_back(field);
      s_.push_back(mfs_.back()->getRawPtr());
    }
  }



  /// \brief \f[ *this = \alpha a b + \beta *this \f]
  ///
  /// \param alpha scalar
  /// \param a Vector
  /// \param b Matrix
  /// \param beta
  void TimesMatAdd(const ST alpha, const FieldT& a,
                   const Teuchos::SerialDenseMatrix<int, ST>& b,
                   const ST beta) {

    int m1 = a.getNumberVecs(); ///<is assumed to be equal to number vecs of this and ncolumns and nrows of b
    int m2 = getNumberVecs();   ///<is assumed to be equal to number vecs of this and ncolumns and nrows of b

    OT nx = mfs_[0]->getStorageSize();

    for(int i=0; i<m2; ++i) {
      for(OT k=0; k<nx; ++k) {
        ST temp = 0.;
        for(int j=0; j<m1; ++j)
          temp += alpha*b(j, i)*a.s_[j][k] ;
        s_[i][k] = s_[i][k]*beta + temp;
      }
      mfs_[i]->changed();
    }
  }



  /// \brief <tt>mv := alpha*a + beta*b</tt>
  void add(ST alpha, const FieldT& a, ST beta, const FieldT& b, const B wb=B::Y) {

    for(int i=0; i<getNumberVecs(); ++i)
      mfs_[i]->add(alpha, *a.mfs_[i], beta, *b.mfs_[i], wb);
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void abs(const FieldT& y) {

    for(int i=0; i<getNumberVecs(); ++i)
      mfs_[i]->abs(*y.mfs_[i]);
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1, \dots, n  \f]
  /// \return Reference to this object
  void reciprocal(const FieldT& y) {

    for(int i=0; i<getNumberVecs(); ++i)
      mfs_[i]->reciprocal(*y.mfs_[i]);
  }


  /// \brief Scale each element of every \c Field by \c gamma.
  ///
  /// Here x represents on \c Field, and we update it as
  /// \f[ x_i = \alpha x_i \quad \mbox{for } i=1, \dots, n \f]
  void scale(const ST alpha) {

    for(int i=0; i<getNumberVecs(); ++i)
      mfs_[i]->scale(alpha);
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void scale(const FieldT& a) {

    for(int i=0; i<getNumberVecs(); ++i)
      mfs_[i]->scale(*a.mfs_[i]);
  }


  /// \brief Scale each element of a \c Field by a \c gamma.
  ///
  /// Here x_j represents the j'th field, and we update it as
  /// \f[ x_j[i] = \alpha_j x_j[i] \quad \mbox{for } i=1, \dots, n \f]
  void scale(const std::vector<ST>& alphas) {

    assert(alphas.size()==getNumberVecs());

    for(unsigned int i=0; i<alphas.size(); ++i)
      mfs_[i]->scale(alphas[i]);
  }


  /// \f[ C_{ij} = \alpha A_i^T \c this_j \f]
  /// \param alpha
  /// \param A
  /// \param C
  /// \note copying "twice" only necessary if there is a stride in C
  void Trans(const ST alpha, const FieldT& A,
             Teuchos::SerialDenseMatrix<int, ST>& C) const {

    const int n1 = getNumberVecs();
    const int n2 = A.getNumberVecs();

    Teuchos::SerialDenseMatrix<int, ST> Cloc(n1, n2, true);
    Teuchos::SerialDenseMatrix<int, ST> Cglo(n1, n2, true);

    assert(n1==C.numRows());
    assert(n2==C.numCols());

    for(int i=0; i<n1; ++i)
      for(int j=0; j<n2; ++j)
        Cloc(i, j) = alpha*mfs_[i]->dotLoc(*A.mfs_[j]);

    MPI_Allreduce(Cloc.values(), Cglo.values(), n1*n2, MPI_REAL8, MPI_SUM, comm());

    for(int i=0; i<n1; ++i)
      for(int j=0; j<n2; ++j)
        C(i, j) = Cglo(i, j);
  }


  /// \}

  /// For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
  void dot(const FieldT& A, std::vector<ST>& dots) const {

    const int n = getNumberVecs();
    ST* temp = new ST[n];

    for(int i=0; i<n; ++i)
      temp[i] = A.mfs_[i]->dotLoc(*mfs_[i]);

    MPI_Allreduce(temp, dots.data(), n, MPI_REAL8, MPI_SUM, comm());
    delete[] temp;
  }


  /// \brief Compute the inner product for the \c MultiField considering it as one Vector.
  constexpr ST dotLoc(const FieldT& A) {
    int n = getNumberVecs();

    ST b = 0.;

    for(int i=0; i<n; ++i)
      b+= mfs_[i]->dotLoc(*A.mfs_[i]);

    return b;
  }

  /// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this,
  /// i.e.\f$b = y^H this\f$.
  constexpr ST dot(const FieldT& y) {

    return this->reduce(comm(), dotLoc(y));
  }

  /// \brief Compute the norm of each individual vector.
  /// Upon return, \c normvec[i] holds the value of \f$||this_i||_2\f$, the \c i-th column of \c this.
  void norm(std::vector<typename Teuchos::ScalarTraits<ST>::magnitudeType> &normvec,
      const ENorm type=ENorm::Two) const {

    const int n = getNumberVecs();
    ST* temp = new ST[n];

    for(int i=0; i<n; ++i)
      temp[i] = mfs_[i]->normLoc(type);

    switch(type) {
      case ENorm::One: {
        MPI_Allreduce(temp, normvec.data(), n, MPI_REAL8, MPI_SUM, comm());
        break;
      }
      case ENorm::Two: {
        MPI_Allreduce(temp, normvec.data(), n, MPI_REAL8, MPI_SUM, comm());
        for(int i=0; i<n; ++i)
          normvec[i] = std::sqrt(normvec[i]);
        break;
      }
      case ENorm::Inf: {
        MPI_Allreduce(temp, normvec.data(), n, MPI_REAL8, MPI_MAX, comm());
        break;
      }
      case ENorm::L2: {
        MPI_Allreduce(temp, normvec.data(), n, MPI_REAL8, MPI_SUM, comm());
        for(int i=0; i<n; ++i)
          normvec[i] = std::sqrt(normvec[i]);
        break;
      }
    }
    delete[] temp;
  }


  /// \brief Compute the norm for the \c MultiField as it is considered as one Vector .
  constexpr ST normLoc(ENorm type = ENorm::Two) {

    ST normvec = 0.;

    for(int i=0; i<getNumberVecs(); ++i)
      normvec =
        (ENorm::Inf==type)?
        (std::max(mfs_[i]->normLoc(type), normvec)):
        (normvec+mfs_[i]->normLoc(type));

    return normvec;
  }

/// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  constexpr ST norm(ENorm type = ENorm::Two) {

    ST normvec = this->reduce(
                   comm(),
                   normLoc(type),
                   (ENorm::Inf==type)?MPI_MAX:MPI_SUM);

    normvec =
      (ENorm::Two==type || ENorm::L2==type) ?
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

    ST nor = Teuchos::ScalarTraits<ST>::zero();

    for(int i=0; i<getNumberVecs(); ++i)
      nor += mfs_[i]->normLoc(*weights.mfs_[i]);

    return nor;
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


  /// \test this
  /// \param A
  /// \param index
  void SetBlock(const FieldT& A, const std::vector<int>& index) {

    const int n = index.size();
    for(int i=0; i<n; ++i)
      *mfs_[index[i]] = *A.mfs_[i];
  }

  /// \param A
  /// \param index
  void SetBlock(const FieldT& A, const Teuchos::Range1D& index) {

    for(int i=index.lbound(); i<=index.ubound(); ++i)
      *mfs_[i] = *A.mfs_[i-index.lbound()];
  }

  /// \brief *this := a.
  ///
  /// assign (deep copy) a into *this.
  MultiField& operator=(const MultiField& a) {

    const int n = getNumberVecs();
    for(int i=0; i<n; ++i)
      *mfs_[i] = *a.mfs_[i];

    return *this;
  }

  /// \brief Replace the vectors with a random vectors.
  void random(const bool useSeed = false, const int seed = 1) {

    const int n = getNumberVecs();
    for(int i=0; i<n; ++i)
      mfs_[i]->random();
  }

  /// \brief \f[ *this = \alpha \f]
  void init(const ST alpha = Teuchos::ScalarTraits<ST>::zero(), const B wB=B::Y) {

    const int n = getNumberVecs();
    for(int i=0; i<n; ++i)
      mfs_[i]->init(alpha, wB);
  }

  void extrapolateBC(const Belos::ETrans trans=Belos::NOTRANS) {
    
    const int n = getNumberVecs();
    for(int i=0; i<n; ++i)
      mfs_[i]->extrapolateBC(trans);
  }

  void level() const {

    const int n = getNumberVecs();
    for(int i=0; i<n; ++i)
      mfs_[i]->level();
  }

  /// \param out out stream
  void print(std::ostream& out=std::cout) const {
    const int n = getNumberVecs();
    for(int i=0; i<n; ++i)
      mfs_[i]->print(out);
  }


  void write(const int count=0, const bool restart=false) const {
    const int n = getNumberVecs();
    for(int i=0; i<n; ++i)
      mfs_[i]->write(count + i, restart);
  }


  void read(const int count=0) {
    const int n = getNumberVecs();
    for(int i=0; i<n; ++i)
      mfs_[i]->read(count + i);
  }


  /// \name Attribute methods
  /// @{

  InnerFieldT& getField(const int i) {
    return *mfs_[i];
  }

  constexpr const InnerFieldT& getField(const int i) {
    return *mfs_[i];
  }

  constexpr const Teuchos::RCP<const GridT>& grid() {
    return AF::grid_;
  }

  constexpr const MPI_Comm& comm() {
    return mfs_[0]->comm();
  }

  ///  @}

}; // end of class MultiField



/// \brief factory for \c MultiField.
///
/// simple wrapper.
/// \relates MultiField
///
/// \todo make field reference + const cast(think of giving to flavours
template<class FieldT>
Teuchos::RCP<MultiField<FieldT> > wrapMultiField(const Teuchos::RCP<FieldT>& field) {

  return Teuchos::rcp(new MultiField<FieldT>(field));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIFIELD_HPP
