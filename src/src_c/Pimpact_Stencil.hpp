#pragma once
#ifndef PIMPACT_STENCIL_HPP
#define PIMPACT_STENCIL_HPP


#include <iostream>

#include "Teuchos_ScalarTraits.hpp"



namespace Pimpact {


/// \brief Array
///
/// \tparam Scalar
/// \tparam Ordinal
/// \tparam ss start index
template<class Scalar, class Ordinal, int ss >
class Array {

protected:

  using ScalarArray = Scalar*;

  ScalarArray c_;

  Ordinal nn_;

  constexpr Ordinal size() {
    return nn_-ss+1;
  }

public:

  Array() : c_(nullptr), nn_(ss) {}

  Array(const Ordinal nn) : nn_(nn) {

    assert(ss<=nn_);

    c_ = new Scalar[ size() ];

    std::fill_n(c_, size(), Teuchos::ScalarTraits<Scalar>::zero());
  }

  Array(const Array& that) : Array(that.nn_) {

    std::copy_n(that.c_, size(), c_);
  }

  friend void swap(Array& one, Array& two) {
    using std::swap;
    swap(one.nn_, two.nn_);
    swap(one.c_, two.c_);
  }

  Array(Array&& that) : Array() {
    swap(*this, that);
  }

  Array& operator=(Array that) {
    swap(*this, that);
    return *this;
  }

  ~Array() {
    delete[] c_;
  }

  constexpr ScalarArray get() const {
    return c_;
  }

  constexpr const Ordinal NN() {
    return nn_;
  }

  Scalar& operator[](const Ordinal index) {
    assert(index>=ss);
    assert(index<=nn_);
    return c_[ index-ss ];
  };

  constexpr Scalar operator[](const Ordinal index) {
    assert(index>=ss);
    assert(index<=nn_);
    return c_[ index-ss ];
  };


  void print(std::ostream& out=std::cout) const {

    out << std::scientific;
    out << std::setw(3) << "i" << std::setw(14) << "x" << "\n";
    for(Ordinal i=ss; i<=nn_; ++i) {
      out << std::setw(3) << i << std::setw(14) << (*this)[i] << "\n";
    }
  }

}; // end of class Array



/// \brief Stencil
///
/// \tparam Scalar
/// \tparam Ordinal
/// \tparam ss start index
/// \tparam lb lower stencil bound
/// \tparam ub_ upper stencil bound
template<class Scalar, class Ordinal, int ss, int lb, int ub>
class Stencil {

protected:

  using ScalarArray = Scalar*;

  ScalarArray c_;

  Ordinal nn_;
  static const int w_ = ub-lb+1;

  constexpr Ordinal size() {
    return (nn_ - ss + 1)*w_;
  }

public:

  Stencil() : c_(nullptr), nn_(ss) {}

  Stencil(const Ordinal nn) : nn_(nn) {

    static_assert(lb<=ub, "Stencil width cannot be negative");
    assert(ss<=nn_);

    c_ = new Scalar[ size() ];

    std::fill_n(c_, size(), Teuchos::ScalarTraits<Scalar>::zero());
  }

  Stencil(const Stencil& that) : Stencil(that.nn_) {

    std::copy_n(that.c_, size(), c_);
  }

  friend void swap(Stencil& one, Stencil& two) {
    using std::swap;
    swap(one.nn_, two.nn_);
    swap(one.c_, two.c_);
  }

  Stencil(Stencil&& that) : Stencil() {
    swap(*this, that);
  }

  Stencil& operator=(Stencil that) {
    swap(*this, that);
    return *this;
  }

  ~Stencil() {
    delete[] c_;
  }

  constexpr ScalarArray get() const {
    return c_;
  }

  Scalar& operator()(const Ordinal index, const int offset) {

    assert(offset>=lb);
    assert(offset<=ub);
    assert(index>=ss);
    assert(index<=nn_);

    return c_[ offset-lb + (index-ss)*w_ ];
  };

  constexpr Scalar operator()(const Ordinal index, const int offset) {

    assert(offset>=lb);
    assert(offset<=ub);
    assert(index>=ss);
    assert(index<=nn_);

    return c_[ offset-lb + (index-ss)*w_ ];
  };


  static constexpr int bl() {
    return lb;
  }

  static constexpr int bu() {
    return ub;
  }

  constexpr const Ordinal NN() {
    return nn_;
  }

  void print(std::ostream& out=std::cout) const {
    out << std::setw(8) << "bl: ";
    out << std::scientific;
    out << std::setprecision(3);
    for(int k=lb; k<=ub; ++k)
      out << std::setw(12) << k ;
    out << " :bu\n";

    for(int bla=0; bla<12*w_+5; bla++)
      out << "-";
    out << "\n";

    for(Ordinal i=ss; i<=nn_; ++i) {
      out << "i: " << std::setw(3) << i << " (";
      for(int ii=lb; ii<=ub; ++ii)
        out << std::setw(12) << (*this)(i, ii);
      out << ")\n";
    }
  }

}; // end of class Stencil


} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_STENCIL_HPP
