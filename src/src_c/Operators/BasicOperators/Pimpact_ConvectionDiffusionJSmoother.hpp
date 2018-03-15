#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP
#define PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP

#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_DiffusionOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"



namespace Pimpact {


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
/// \note todo merge with SORSmoother or make interface
template<class OperatorT>
class ConvectionDiffusionJSmoother {

public:

  using GridT = typename OperatorT::GridT;

  using FluxFieldT = ScalarField<GridT>[3];

  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using SW = typename GridT::SW;

  ST omega_;
  int nIter_;

  const Teuchos::RCP<const OperatorT> op_;

  constexpr ST getHC(const ECoord dir, const F ftype, OT i, OT ii) {
    return op_->getHelmOp()->getC(dir, ftype, i, ii);
  }

public:

  /// \brief constructor
  ///
  /// These options include the following:
  /// - "omega" - damping parameter
  /// - "numIters" - an \c int specifying the maximum number of iterations the
  ConvectionDiffusionJSmoother(
    const Teuchos::RCP<const OperatorT>& op,
    Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList()):
    omega_(pl->get<ST>("omega", 0.5)),
    nIter_(pl->get("numIters", 10)),
    op_(op) {}



  void apply(const FluxFieldT& wind, const DomainFieldT& x, RangeFieldT& y, ST mul, ST mulI, ST mulC, ST mulL) const {
    std::cout << "not implmented\n";
  }


protected:

  void applyStep(const FluxFieldT& wind, const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y) const {

    const F f = y.getType();

    x.exchange();

    //applyBC(b, x, y);

    if(3==GridT::sdim) {
      for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
        for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST diag =
              op_->getMulI()
              + op_->getMulC() * op_->getConvSOp()->innerDiag3D(
                wind[0](i, j, k),
                wind[1](i, j, k),
                wind[2](i, j, k), f, i, j, k)
              - op_->getMulL() * op_->getHelmOp()->innerDiag3D(f, i, j, k) ;
            assert(diag!=0);
            y(i, j, k) = x(i, j, k) + omega_*(b(i, j, k)
                                           - op_->getMulI() * x(i, j, k)
                                           - op_->getMulC() * op_->getConvSOp()->innerStenc3D(
                                             wind[0](i, j, k),
                                             wind[1](i, j, k),
                                             wind[2](i, j, k), x, i, j, k)
                                           + op_->getMulL() * op_->getHelmOp()->innerStenc3D(x, f, i, j, k)) / diag;
          }
    } else {

      for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
        for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST diag =
              op_->getMulI()
              + op_->getMulC() * op_->getConvSOp()->innerDiag2D(
                wind[0](i, j, k),
                wind[1](i, j, k), f, i, j, k)
              - op_->getMulL() * op_->getHelmOp()->innerDiag2D(f, i, j, k) ;
            assert(diag!=0);
            y(i, j, k) = x(i, j, k) + omega_*(b(i, j, k)
                                           - op_->getMulI() * x(i, j, k)
                                           - op_->getMulC() * op_->getConvSOp()->innerStenc2D(
                                             wind[0](i, j, k),
                                             wind[1](i, j, k), x, i, j, k)
                                           + op_->getMulL() * op_->getHelmOp()->innerStenc2D(x, f, i, j, k)) / diag;
          }
    }

    applyBC(b, x, y);
    y.changed();
  }

  /// \brief implements smoothing for Dirichlet boundary conditions as identity
  /// in tangential/velocity direction or interpolation in wand normal
  /// direction
  void applyBC(const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y) const {

    assert(b.getType()==y.getType());
    assert(x.getType()==y.getType());

    const F f = y.getType();

    const ST omegaBC = omega_;

    // U-field
    if(F::U==f) {

      // tangential direction: Y
      if(0<grid()->bcl(Y)) {
        OT j = grid()->si(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST temp = 0.;
            for(OT jj=0; jj<=SW::BU(Y); ++jj)
              temp += getHC(Y, f, j, jj)*x(i, j+jj, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }
      if(0<grid()->bcu(Y)) {
        OT j = grid()->ei(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST temp = 0.;
            for(OT jj=SW::BL(Y); jj<=0; ++jj)
              temp += getHC(Y, f, j, jj)*x(i, j+jj, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }

      // tangential direction: Z
      if(0<grid()->bcl(Z)) {
        OT k = grid()->si(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST temp = 0.;
            for(OT kk=0; kk<=SW::BU(Z); ++kk)
              temp += getHC(Z, f, k, kk)*x(i, j, k+kk);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }
      if(0<grid()->bcu(Z)) {
        OT k = grid()->ei(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST temp = 0.;
            for(OT kk=SW::BL(Z); kk<=0; ++kk)
              temp += getHC(Z, f, k, kk)*x(i, j, k+kk);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }

      // normal direction: X
      if(0<grid()->bcl(X)) {
        OT i = grid()->si(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            ST temp = 0.;
            for(OT ii=0; ii<=SW::BU(X); ++ii)
              temp += getHC(X, f, i, ii)*x(i+ii, j, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }
      if(0<grid()->bcu(X)) {
        OT i = grid()->ei(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            ST temp = 0.;
            for(OT ii=SW::BL(X); ii<=0; ++ii)
              temp += getHC(X, f, i, ii)*x(i+ii, j, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }
    }

    // V-field
    if(F::V==f) {

      // tangential direction: X
      if(0<grid()->bcl(X)) {
        OT i = grid()->si(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j) {
            ST temp = 0.;
            for(OT ii=0; ii<=SW::BU(X); ++ii)
              temp += getHC(X, f, i, ii)*x(i+ii, j, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }
      if(0<grid()->bcu(X)) {
        OT i = grid()->ei(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j) {
            ST temp = 0.;
            for(OT ii=SW::BL(X); ii<=0; ++ii)
              temp += getHC(X, f, i, ii)*x(i+ii, j, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }

      // tangential direction: Z
      if(0<grid()->bcl(Z)) {
        OT k = grid()->si(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT kk=0; kk<=SW::BU(Z); ++kk)
              temp += getHC(Z, f, k, kk)*x(i, j, k+kk);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }
      if(0<grid()->bcu(Z)) {
        OT k = grid()->ei(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT kk=SW::BL(Z); kk<=0; ++kk)
              temp += getHC(Z, f, k, kk)*x(i, j, k+kk);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }

      // normal direction: Y
      if(0<grid()->bcl(Y)) {
        OT j = grid()->si(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT jj=0; jj<=SW::BU(Y); ++jj)
              temp += getHC(Y, f, j, jj)*x(i, j+jj, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }
      if(0<grid()->bcu(Y)) {
        OT j = grid()->ei(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT jj=SW::DL(Y); jj<=SW::DU(Y); ++jj)
              temp += getHC(Y, f, j, jj)*x(i, j+jj, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }
    }

    // W-field
    if(F::W==f) {

      // tangential direction: X
      if(0<grid()->bcl(X)) {
        OT i = grid()->si(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            ST temp = 0.;
            for(OT ii=0; ii<=SW::BU(X); ++ii)
              temp += getHC(X, f, i, ii)*x(i+ii, j, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }
      if(0<grid()->bcu(X)) {
        OT i = grid()->ei(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            ST temp = 0.;
            for(OT ii=SW::BL(X); ii<=0; ++ii)
              temp += getHC(X, f, i, ii)*x(i+ii, j, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }

      // tangential direction: Y
      if(0<grid()->bcl(Y)) {
        OT j = grid()->si(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT jj=0; jj<=SW::BU(Y); ++jj)
              temp += getHC(Y, f, j, jj)*x(i, j+jj, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }
      if(0<grid()->bcu(Y)) {
        OT j = grid()->ei(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT jj=SW::BL(Y); jj<=0; ++jj)
              temp += getHC(Y, f, j, jj)*x(i, j+jj, k);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }

      // normal direction: Z
      if(0<grid()->bcl(Z)) {
        OT k = grid()->si(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT kk=0; kk<=SW::BU(Z); ++kk)
              temp += getHC(Z, f, k, kk)*x(i, j, k+kk);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }
      if(0<grid()->bcu(Z)) {
        OT k = grid()->ei(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT kk=SW::BL(Z); kk<=0; ++kk)
              temp += getHC(Z, f, k, kk)*x(i, j, k+kk);
            y(i, j, k) = x(i, j, k) + omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }
    }
  }

public:

  void apply(const FluxFieldT& wind, const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N) const {

    const F m = y.getType();

    DomainFieldT temp(grid(), Owning::Y, m);

    assert(y.getType() == x.getType());

    for(int i =0; i<GridT::sdim; ++i)
      assert(wind[i].getType() == x.getType());

    for(int vel_dir=0; vel_dir<GridT::sdim; ++vel_dir)
      wind[vel_dir].exchange();

    for(int i=0; i<nIter_; ++i) {
      applyStep(wind, x, y, temp);
      applyStep(wind, x, temp, y);
    }
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  void print(std::ostream& out=std::cout) const {
    out << "--- " << getLabel() << " ---\n";
    op_->print();
  }


  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const std::string getLabel() const {
    return "ConvectionDiffusionJSmoother ";
  };

}; // end of class ConvectionDiffusionJSmoother


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP
