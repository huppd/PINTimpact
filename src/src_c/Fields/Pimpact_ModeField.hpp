/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_MODEFIELD_HPP
#define PIMPACT_MODEFIELD_HPP


#include <vector>
#include <iostream>

#include "mpi.h"

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_AbstractField.hpp"




namespace Pimpact {



/// \brief important basic Vector class
/// vector for wrapping 2 fields into one mode
/// \ingroup Field
template<class IFT>
class ModeField : private AbstractField<typename IFT::GridT> {

public:

 using GridT = typename IFT::GridT;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using ScalarArray = ST*;

  using AF = AbstractField<GridT>;

  const Owning owning_;

  IFT fieldc_;
  IFT fields_;

  ScalarArray s_;

private:

  void allocate() {
    OT n = getStorageSize();
    setStoragePtr(new ST[n]);
    std::uninitialized_fill_n(s_, n , 0.);
  }

public:

  constexpr OT getStorageSize() {
    return fieldc_.getStorageSize() + fields_.getStorageSize();
  }

  constexpr ST* getRawPtr() {
    return s_;
  }

  void setStoragePtr(ST* array) {
    s_ = array;
    fieldc_.setStoragePtr(s_);
    fields_.setStoragePtr(s_ + fieldc_.getStorageSize());
  }


  ModeField(const Teuchos::RCP<const GridT>& grid, const Owning owning=Owning::Y):
    AF(grid),
    owning_(owning),
    fieldc_(grid, Owning::N),
    fields_(grid, Owning::N) {

    if(owning_==Owning::Y) allocate();
  };


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ECopy::Shallow is done but allows also to deepcopy the field
  ModeField(const ModeField& vF, const ECopy copyType=ECopy::Deep):
    AF(vF.grid()),
    owning_(vF.owning_),
    fieldc_(vF.fieldc_, copyType),
    fields_(vF.fields_, copyType) {

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

  ~ModeField() {
    if(owning_==Owning::Y) delete[] s_;
  }

  Teuchos::RCP<ModeField> clone(const ECopy cType=ECopy::Deep) const {

    Teuchos::RCP<ModeField> mv = Teuchos::rcp(new ModeField(grid()));

    switch(cType) {
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

  IFT& getCField() {
    return fieldc_;
  }
  IFT& getSField() {
    return fields_;
  }

  constexpr const IFT& getCField() {
    return fieldc_;
  }
  constexpr const IFT& getSField() {
    return fields_;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() {
    return AF::grid_;
  }

  constexpr const MPI_Comm& comm() {
    return fieldc_.comm();
  }

  /// \brief returns the length of Field.
  constexpr OT getLength() {
    return fieldc_.getLength() + fields_.getLength();
  }



  /// \}
  /// \name Update methods
  /// \{

  /// \brief Replace \c this with \f$\alpha a + \beta b\f$.
  void add(const ST alpha, const ModeField& a, const ST beta, const ModeField& b,
      const B wb=B::Y) {

    fieldc_.add(alpha, a.fieldc_, beta, b.fieldc_, wb);
    fields_.add(alpha, a.fields_, beta, b.fields_, wb);
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void abs(const ModeField& y) {

    fieldc_.abs(y.fieldc_);
    fields_.abs(y.fields_);
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1, \dots, n  \f]
  /// \return Reference to this object
  void reciprocal(const ModeField& y) {

    fieldc_.reciprocal(y.fieldc_);
    fields_.reciprocal(y.fields_);
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale(const ST alpha, const B wB=B::Y) {

    fieldc_.scale(alpha, wB);
    fields_.scale(alpha, wB);
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void scale(const ModeField& a) {

    fieldc_.scale(a.fieldc_);
    fields_.scale(a.fields_);
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  constexpr ST dotLoc(const ModeField& a) {

    ST b=0.;

    b = fieldc_.dotLoc(a.fieldc_) + fields_.dotLoc(a.fields_);

    return b;
  }


  /// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
  constexpr ST dot(const ModeField& y) {

    return this->reduce(comm(), dotLoc(y));
  }

  /// \}
  /// \name Norm method
  /// \{

  constexpr ST normLoc(ENorm type=ENorm::Two) {

    ST normvec =
      (ENorm::Inf==type)?
      std::fmax(fieldc_.normLoc(type), fields_.normLoc(type)):
      (fieldc_.normLoc(type) + fields_.normLoc(type));

    return normvec;
  }

  /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  constexpr ST norm(const ENorm type=ENorm::Two) {

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
  constexpr ST normLoc(const ModeField& weights) {
    return fieldc_.normLoc(weights.fieldc_) + fields_.normLoc(weights.fields_);
  }

  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST norm(const ModeField& weights) {

    return std::sqrt(this->reduce(comm(), normLoc(weights)));
  }


  /// \}
  /// \name Initialization methods
  /// \{

  /// \brief *this := a
  /// Assign (deep copy) A into mv.
  ModeField& operator=(const ModeField& a) {

    fieldc_ = a.fieldc_;
    fields_ = a.fields_;

    return *this;
  }

  /// \brief Replace the vectors with a random vectors.
  void random(bool useSeed=false, const B bcYes=B::Y, int seed=1) {

    fieldc_.random(useSeed, bcYes, seed);
    fields_.random(useSeed, bcYes, seed);
  }

  /// \brief Replace each element of the vector  with \c alpha.
  void init(const ST alpha = Teuchos::ScalarTraits<ST>::zero(), const B wB=B::Y) {

    fieldc_.init(alpha, wB);
    fields_.init(alpha, wB);
  }

  void extrapolateBC(const Belos::ETrans trans=Belos::NOTRANS) {

    fieldc_.extrapolateBC(trans);
    fields_.extrapolateBC(trans);
  }

  void level() const {

    fieldc_.level();
    fields_.level();
  }

  void changed() const {

    fieldc_.changed();
    fields_.changed();
  }

  /// \}
  /// Print the vector.  To be used for debugging only.
  void print(std::ostream& out=std::cout) const {

    fieldc_.print(out);
    fields_.print(out);
  }

  void write(const int count=0, const bool restart=false) const {

    fieldc_.write(count, restart);
    fields_.write(count+1, restart);
  }

  void read(const int count=0) {

    fieldc_.read(count);
    fields_.read(count+1);
  }

  /// \name comunication methods.
  /// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
  ///
  /// \{

  void exchange() const {

    fieldc_.exchange();
    fields_.exchange();
  }

  void setExchanged() const {

    fieldc_.setExchanged();
    fields_.setExchanged();
  }


  /// \}

}; // end of class ModeField


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MODEFIELD_HPP
