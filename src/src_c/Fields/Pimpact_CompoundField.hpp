#pragma once
#ifndef PIMPACT_COMPOUNDFIELD_HPP
#define PIMPACT_COMPOUNDFIELD_HPP


#include <vector>
#include <iostream>

#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_AbstractField.hpp"




namespace Pimpact {


/// \brief important basic Vector class
/// vector for wrapping 2 fields into one mode
/// \ingroup Field
/// \todo use attributes methods in vectorspace functions???
/// \todo rm RCP, therefore make move constructor
template<class VField, class SField>
class CompoundField : private AbstractField<typename VField::GridT> {

public:

  using GridT = typename VField::GridT;

protected:

  using ST = typename GridT::Scalar;
  using OT =typename GridT::Ordinal;

  using ScalarArray =  ST*;

  using AF = AbstractField<GridT>;

  const Owning owning_;

  VField vfield_;
  SField sfield_;

  ScalarArray s_;

private:

  void allocate() {
    OT n = getStorageSize();
    setStoragePtr(new ST[n]);
    std::uninitialized_fill_n(s_, n , 0.);
  }

public:

  constexpr OT getStorageSize() {
    return vfield_.getStorageSize() + sfield_.getStorageSize();
  }

  constexpr ST* getRawPtr() {
    return s_;
  }

  void setStoragePtr(ST* array) {
    s_ = array;
    vfield_.setStoragePtr(s_);
    sfield_.setStoragePtr(s_ + vfield_.getStorageSize());
  }


  CompoundField(const Teuchos::RCP<const GridT>& grid, const Owning owning=Owning::Y):
    AF(grid),
    owning_(owning),
    vfield_(grid, Owning::N),
    sfield_(grid, Owning::N) {

    if(owning_==Owning::Y) allocate();
  };


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param field
  /// \param copyType by default a ECopy::Shallow is done but allows also to deepcopy the field
  CompoundField(const CompoundField& field, const ECopy copyType=ECopy::Deep):
    AF(field.grid()),
    owning_(field.owning_),
    vfield_(field.vfield_, copyType),
    sfield_(field.sfield_, copyType) {

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
    };


  ~CompoundField() {
    if(owning_==Owning::Y) delete[] s_;
  }
  
  Teuchos::RCP<CompoundField> clone(const ECopy ctype=ECopy::Deep) const {

    Teuchos::RCP<CompoundField> mv = Teuchos::rcp(new CompoundField(grid()));

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

  VField& getVField() {
    return vfield_;
  }
  SField& getSField() {
    return sfield_;
  }

  constexpr const VField& getVField() {
    return vfield_;
  }
  constexpr const SField& getSField() {
    return sfield_;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() {
    return AF::grid_;
  }

  constexpr const MPI_Comm& comm() {
    return vfield_.comm();
  }


  /// \brief get Vect length
  /// shoud be the same as 2*vfield_->getVecLength()
  /// the vector length is withregard to the inner points
  /// \return vector length
  /// \brief returns the length of Field.
  constexpr OT getLength() {
    return vfield_.getLength() + sfield_.getLength();
  }



  /// \}
  /// \name Update methods
  /// \{

  /// \brief Replace \c this with \f$\alpha a + \beta b\f$.
  void add(const ST alpha, const CompoundField& a, const ST beta, const CompoundField& b,
      const B wb=B::Y) {

    // add test for consistent Grids in debug mode
    vfield_.add(alpha, a.vfield_, beta, b.vfield_, wb);
    sfield_.add(alpha, a.sfield_, beta, b.sfield_, wb);
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void abs(const CompoundField& y) {

    vfield_.abs(y.vfield_);
    sfield_.abs(y.sfield_);
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1, \dots, n  \f]
  /// \return Reference to this object
  void reciprocal(const CompoundField& y) {

    vfield_.reciprocal(y.vfield_);
    sfield_.reciprocal(y.sfield_);
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale(const ST alpha) {

    vfield_.scale(alpha);
    sfield_.scale(alpha);
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void scale(const CompoundField& a) {

    vfield_.scale(a.vfield_);
    sfield_.scale(a.sfield_);
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  constexpr ST dotLoc(const CompoundField& a) {

    ST b = 0.;

    b = vfield_.dotLoc(a.vfield_) + sfield_.dotLoc(a.sfield_);

    return b;
  }


  /// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this,
  /// i.e.\f$b = y^H this\f$.
  constexpr ST dot(const CompoundField& y) {

    return this->reduce(comm(), dotLoc(y));
  }


  /// \}
  /// \name Norm method
  /// \{

  /// \brief Compute the norm of the field.
  constexpr ST normLoc(const ENorm type=ENorm::Two) {

    return (ENorm::Inf==type)?
      std::max(vfield_.normLoc(type), sfield_.normLoc(type)):
      (vfield_.normLoc(type) + sfield_.normLoc(type));
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
  constexpr ST normLoc(const CompoundField& weights) {

    return vfield_.normLoc(weights.vfield_) + sfield_.normLoc(weights.sfield_);
  }

  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST norm(const CompoundField& weights) {

    return std::sqrt(this->reduce(comm(), normLoc(weights)));
  }


  /// \}
  /// \name Initialization methods
  /// \{

  /// \brief *this := a
  /// Assign (deep copy) a into mv.
  CompoundField& operator=(const CompoundField& a) {

    vfield_ = a.vfield_;
    sfield_ = a.sfield_;

    return *this;
  }

  /// \brief Replace the vectors with a random vectors.
  void random(const bool useSeed = false, const int seed=1) {
    vfield_.random();
    sfield_.random();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  void init(const ST alpha=Teuchos::ScalarTraits<ST>::zero(), const B wB=B::Y) {
    vfield_.init(alpha, wB);
    sfield_.init(alpha, wB);
  }


  void extrapolateBC(const Belos::ETrans trans=Belos::NOTRANS) {
    vfield_.extrapolateBC(trans);
    sfield_.extrapolateBC(trans);
  }

  void level() const {
    vfield_.level();
    sfield_.level();
  }

  void changed() {
    vfield_.changed();
    sfield_.changed();
  }


  /// \}

  /// Print the vector.  To be used for debugging only.
  void print(std::ostream& out=std::cout) const {
    vfield_.print(out);
    sfield_.print(out);
  }


  void write(const int count=0, const bool restart=false) const {
    vfield_.write(count, restart);
    sfield_.write(count, restart);
  }


  void read(const int count=0) {
    vfield_.read(count);
    sfield_.read(count);
  }


}; // end of class CompoundField


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COMPOUNDFIELD_HPP
