/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_MODESMOOTHER_HPP
#define PIMPACT_MODESMOOTHER_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_ModeField.hpp"




namespace Pimpact {



/// \ingroup ModeOperator
template<class OpT>
class ModeSmoother {

public:

  using GridT = typename OpT::GridT;

  using DomainFieldT = typename OpT::DomainFieldT;
  using RangeFieldT  = typename OpT::RangeFieldT;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using SW = typename GridT::SW;

  using InnerOpT = typename OpT::InnerOpT;

  ST mulI_;
  ST mulC_;
  ST mulL_;

  ST omega_;

  int numIter_;

  int type_;

  Teuchos::RCP<OpT> op_;
  Teuchos::RCP<InnerOpT> innerOp_;

public:

  ModeSmoother(
    const Teuchos::RCP<OpT>& op,
    const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList()):
    mulI_(0.),
    mulC_(1.),
    mulL_(1./op->grid()->getDomainSize()->getRe()),
    omega_(pl->get<ST>("omega", 1.)),
    numIter_(pl->get<int>("numIters", 1)),
    type_(pl->get<int>("type", -1)),
    op_(op),
    innerOp_(op->getInnerOpPtr()) {
      //std::cout << "type: " << type_ << "\n";
      //std::cout << "omega: " << omega_ << "\n";
    };


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    switch(type_) {
      case 2: {
        applyGSplaine(x, y);
        break;
      }
      case 3: {
        applyGSall(x, y);
        break;
      }
      case 4: {
        applyGSSH(x, y);
        break;
      }
      default: {
        applyJ(x, y);
        break;
      }
    }
  }


  void applyJ(const DomainFieldT& x, RangeFieldT& y) const {

    const B wnB = B::N;

    RangeFieldT temp(grid());

    for(F m=F::U; m<GridT::sdim; ++m) {
      for(int iter=0; iter<numIter_; ++iter) {

        y.getCField()(m).exchange();
        y.getSField()(m).exchange();

        for(OT k=grid()->si(m, Z, wnB); k<=grid()->ei(m, Z, wnB); ++k)
          for(OT j=grid()->si(m, Y, wnB); j<=grid()->ei(m, Y, wnB); ++j)
            for(OT i=grid()->si(m, X, wnB); i<=grid()->ei(m, X, wnB); ++i) {
              if(3==GridT::sdim) {

                ST diag = mulC_*innerOp_->getSOp()->getConvSOp()->innerDiag3D(
                    innerOp_->getConvField(m)[0](i, j, k),
                    innerOp_->getConvField(m)[1](i, j, k),
                    innerOp_->getConvField(m)[2](i, j, k), m, i, j, k)
                  - mulL_ * innerOp_->getSOp()->getHelmOp()->innerDiag3D(m, i, j, k) ;

                assert(diag!=0);

                temp.getCField()(m)(i, j, k) = y.getCField()(m)(i, j, k) + omega_*(
                    x.getCField()(m)(i, j, k)
                    - mulI_*y.getSField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc3D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      innerOp_->getConvField(m)[2](i, j, k),
                      y.getCField()(m), i, j, k)
                    +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc3D(
                      y.getCField()(m), m, i, j, k))/diag;
                temp.getSField()(m)(i, j, k) = y.getSField()(m)(i, j, k) + omega_*(
                    x.getSField()(m)(i, j, k)
                    + mulI_*y.getCField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc3D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      innerOp_->getConvField(m)[2](i, j, k),
                      y.getSField()(m), i, j, k)
                    +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc3D(
                      y.getSField()(m), m, i, j, k))/diag;
              }
              else {
                ST diag = mulC_*innerOp_->getSOp()->getConvSOp()->innerDiag2D(
                    innerOp_->getConvField(m)[0](i, j, k),
                    innerOp_->getConvField(m)[1](i, j, k), m, i, j, k)
                  - mulL_ * innerOp_->getSOp()->getHelmOp()->innerDiag2D(m, i, j, k) ;
                assert(diag!=0);

                temp.getCField()(m)(i, j, k) = y.getCField()(m)(i, j, k) + omega_*(
                    x.getCField()(m)(i, j, k)
                    - mulI_*y.getSField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc2D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      y.getCField()(m), i, j, k)
                    +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc2D(
                      y.getCField()(m), m, i, j, k))/diag;
                temp.getSField()(m)(i, j, k) = y.getSField()(m)(i, j, k) + omega_*(
                    x.getSField()(m)(i, j, k)
                    + mulI_*y.getCField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc2D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      y.getSField()(m), i, j, k)
                    + mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc2D(
                      y.getSField()(m), m, i, j, k))/diag;
              }
            }

        applyBCJ(x.getCField()(m), y.getCField()(m), temp.getCField()(m));
        applyBCJ(x.getSField()(m), y.getSField()(m), temp.getSField()(m));

        temp.getCField()(m).changed();
        temp.getSField()(m).changed();

        temp.getCField()(m).exchange();
        temp.getSField()(m).exchange();

        for(OT k=grid()->si(m, Z, wnB); k<=grid()->ei(m, Z, wnB); ++k)
          for(OT j=grid()->si(m, Y, wnB); j<=grid()->ei(m, Y, wnB); ++j)
            for(OT i=grid()->si(m, X, wnB); i<=grid()->ei(m, X, wnB); ++i) {
              if(3==GridT::sdim) {

                ST diag = mulC_*innerOp_->getSOp()->getConvSOp()->innerDiag3D(
                    innerOp_->getConvField(m)[0](i, j, k),
                    innerOp_->getConvField(m)[1](i, j, k),
                    innerOp_->getConvField(m)[2](i, j, k), m, i, j, k)
                  - mulL_ * innerOp_->getSOp()->getHelmOp()->innerDiag3D(m, i, j, k) ;

                assert(diag!=0);

                y.getCField()(m)(i, j, k) = temp.getCField()(m)(i, j, k) + omega_ * (
                    x.getCField()(m)(i, j, k)
                    - mulI_*temp.getSField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc3D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      innerOp_->getConvField(m)[2](i, j, k),
                      temp.getCField()(m), i, j, k)
                    +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc3D(
                      temp.getCField()(m), m, i, j, k))/diag;
                y.getSField()(m)(i, j, k) = temp.getSField()(m)(i, j, k) + omega_*(
                    x.getSField()(m)(i, j, k)
                    + mulI_*temp.getCField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc3D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      innerOp_->getConvField(m)[2](i, j, k),
                      temp.getSField()(m), i, j, k)
                    + mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc3D(
                      temp.getSField()(m), m, i, j, k))/diag;
              }
              else {
                ST diag = mulC_*innerOp_->getSOp()->getConvSOp()->innerDiag2D(
                    innerOp_->getConvField(m)[0](i, j, k),
                    innerOp_->getConvField(m)[1](i, j, k), m, i, j, k)
                  - mulL_ * innerOp_->getSOp()->getHelmOp()->innerDiag2D(m, i, j, k) ;
                assert(diag!=0);

                y.getCField()(m)(i, j, k) = temp.getCField()(m)(i, j, k) + omega_*(
                    x.getCField()(m)(i, j, k)
                    - mulI_*temp.getSField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc2D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      temp.getCField()(m), i, j, k)
                    + mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc2D(
                      temp.getCField()(m), m, i, j, k))/diag;

                y.getSField()(m)(i, j, k) = temp.getSField()(m)(i, j, k) + omega_*(
                    x.getSField()(m)(i, j, k)
                    + mulI_*temp.getCField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc2D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      temp.getSField()(m), i, j, k)
                    + mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc2D(
                      temp.getSField()(m), m, i, j, k))/diag;
              }
            }

        applyBCJ(x.getCField()(m), temp.getCField()(m), y.getCField()(m));
        applyBCJ(x.getSField()(m), temp.getSField()(m), y.getSField()(m));

        y.getCField()(m).changed();
        y.getSField()(m).changed();
      }
    }
  }


  void applyGSplaine(const DomainFieldT& x, RangeFieldT& y) const {

    Teuchos::Tuple<bool, 4> downwinding = Teuchos::tuple(true, true, true, true);
    for(int iter=0; iter<numIter_; ++iter) {
      applyGS(x, y, downwinding) ;
    }
  }


  void applyGSall(const DomainFieldT& x, RangeFieldT& y) const {

    for(int iter=0; iter<numIter_; ++iter) {

      Teuchos::Tuple<bool, 4> downwinding = Teuchos::tuple(true, true, true, true);

      for(int i=0; i<2; ++i) {
        downwinding[0] = i;
        for(int j=0; j<2; ++j) {
          downwinding[1] = j;
          for(int k=0; k<2; ++k) {
            downwinding[2] = k;
            for(int l=0; l<2; ++l) {
              downwinding[3] = l;

              //if(0==grid()->rankST())
                //std::cout << downwinding << "\n";

              applyGS(x, y, downwinding) ;
            }
          }
        }
      }
    }
    //if(0==grid()->rankST())
      //std::cout << "\n";
  }


  void applyGSSH(const DomainFieldT& x, RangeFieldT& y) const {

    for(int iter=0; iter<numIter_; ++iter) {

      Teuchos::Tuple<bool, 4> downwinding = Teuchos::tuple(false, true, true, true);

      for(int i=0; i<2; ++i) {
        downwinding[0] = i;
        for(int k=0; k<2; ++k) {
          downwinding[2] = k;
          for(int l=0; l<2; ++l) {
            downwinding[3] = l;

            applyGS(x, y, downwinding) ;
          }
        }
      }
    }
  }

  void assignField(const DomainFieldT& mv) {};

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op_->grid();
  };

  //Teuchos::RCP<OpT> getOperator() const {
  //return modeOp_;
  //};

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {

    if(para->name()!="Linear Solver") {
      mulI_ = para->get<ST>("mulI");
      mulC_ = para->get<ST>("mulC");
      mulL_ = para->get<ST>("mulL");
    }
  }

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "ModeSmoother";
  };

  void print(std::ostream& out=std::cout) const {
    out << getLabel() << ":\n";
    op_->print(out);
  }

private:

  /// \brief implements smoothing for Dirichlet boundary conditions as identity
  /// in tangential/velocity direction or interpolation in wand normal
  /// direction
  void applyBCJ(const ScalarField<GridT>& b, const ScalarField<GridT>& x,
      ScalarField<GridT>& y) const {

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


  /// \brief implements smoothing for Dirichlet boundary conditions as identity
  /// in tangential/velocity direction or interpolation in wand normal
  /// direction
  void applyBCGS(const ScalarField<GridT>& b, ScalarField<GridT>& y) const {

    assert(b.getType()==y.getType());

    const F f = y.getType();

    const ST omegaBC = 1.;

    // U-field
    if(F::U==f) {

      // tangential direction: Y
      if(0<grid()->bcl(Y)) {
        OT j = grid()->si(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST temp = 0.;
            for(OT jj=0; jj<=SW::BU(Y); ++jj)
              temp += getHC(Y, f, j, jj)*y(i, j+jj, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }
      if(0<grid()->bcu(Y)) {
        OT j = grid()->ei(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST temp = 0.;
            for(OT jj=SW::BL(Y); jj<=0; ++jj)
              temp += getHC(Y, f, j, jj)*y(i, j+jj, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }

      // tangential direction: Z
      if(0<grid()->bcl(Z)) {
        OT k = grid()->si(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST temp = 0.;
            for(OT kk=0; kk<=SW::BU(Z); ++kk)
              temp += getHC(Z, f, k, kk)*y(i, j, k+kk);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }
      if(0<grid()->bcu(Z)) {
        OT k = grid()->ei(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::N); i<=grid()->ei(f, X, B::N); ++i) {
            ST temp = 0.;
            for(OT kk=SW::BL(Z); kk<=0; ++kk)
              temp += getHC(Z, f, k, kk)*y(i, j, k+kk);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }

      // normal direction: X
      if(0<grid()->bcl(X)) {
        OT i = grid()->si(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            ST temp = 0.;
            for(OT ii=0; ii<=SW::BU(X); ++ii)
              temp += getHC(X, f, i, ii)*y(i+ii, j, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }
      if(0<grid()->bcu(X)) {
        OT i = grid()->ei(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            ST temp = 0.;
            for(OT ii=SW::BL(X); ii<=0; ++ii)
              temp += getHC(X, f, i, ii)*y(i+ii, j, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
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
              temp += getHC(X, f, i, ii)*y(i+ii, j, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }
      if(0<grid()->bcu(X)) {
        OT i = grid()->ei(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j) {
            ST temp = 0.;
            for(OT ii=SW::BL(X); ii<=0; ++ii)
              temp += getHC(X, f, i, ii)*y(i+ii, j, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }

      // tangential direction: Z
      if(0<grid()->bcl(Z)) {
        OT k = grid()->si(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT kk=0; kk<=SW::BU(Z); ++kk)
              temp += getHC(Z, f, k, kk)*y(i, j, k+kk);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }
      if(0<grid()->bcu(Z)) {
        OT k = grid()->ei(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::N); j<=grid()->ei(f, Y, B::N); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT kk=SW::BL(Z); kk<=0; ++kk)
              temp += getHC(Z, f, k, kk)*y(i, j, k+kk);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }

      // normal direction: Y
      if(0<grid()->bcl(Y)) {
        OT j = grid()->si(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT jj=0; jj<=SW::BU(Y); ++jj)
              temp += getHC(Y, f, j, jj)*y(i, j+jj, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }
      if(0<grid()->bcu(Y)) {
        OT j = grid()->ei(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::Y); k<=grid()->ei(f, Z, B::Y); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT jj=SW::DL(Y); jj<=SW::DU(Y); ++jj)
              temp += getHC(Y, f, j, jj)*y(i, j+jj, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
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
              temp += getHC(X, f, i, ii)*y(i+ii, j, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }
      if(0<grid()->bcu(X)) {
        OT i = grid()->ei(f, X, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j) {
            ST temp = 0.;
            for(OT ii=SW::BL(X); ii<=0; ++ii)
              temp += getHC(X, f, i, ii)*y(i+ii, j, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(X, f, i, 0);
          }
      }

      // tangential direction: Y
      if(0<grid()->bcl(Y)) {
        OT j = grid()->si(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT jj=0; jj<=SW::BU(Y); ++jj)
              temp += getHC(Y, f, j, jj)*y(i, j+jj, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }
      if(0<grid()->bcu(Y)) {
        OT j = grid()->ei(f, Y, B::Y);
        for(OT k=grid()->si(f, Z, B::N); k<=grid()->ei(f, Z, B::N); ++k)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT jj=SW::BL(Y); jj<=0; ++jj)
              temp += getHC(Y, f, j, jj)*y(i, j+jj, k);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Y, f, j, 0);
          }
      }

      // normal direction: Z
      if(0<grid()->bcl(Z)) {
        OT k = grid()->si(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT kk=0; kk<=SW::BU(Z); ++kk)
              temp += getHC(Z, f, k, kk)*y(i, j, k+kk);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }
      if(0<grid()->bcu(Z)) {
        OT k = grid()->ei(f, Z, B::Y);
        for(OT j=grid()->si(f, Y, B::Y); j<=grid()->ei(f, Y, B::Y); ++j)
          for(OT i=grid()->si(f, X, B::Y); i<=grid()->ei(f, X, B::Y); ++i) {
            ST temp = 0.;
            for(OT kk=SW::BL(Z); kk<=0; ++kk)
              temp += getHC(Z, f, k, kk)*y(i, j, k+kk);
            y(i, j, k) += omegaBC*(b(i, j, k) - temp)/getHC(Z, f, k, 0);
          }
      }
    }
  }

  constexpr ST getHC(const ECoord dir, const F ftype, const OT i, const OT ii) {
    return innerOp_->getSOp()->getHelmOp()->getC(dir, ftype, i, ii);
  }


protected:

  void applyGS(
      const DomainFieldT& x, RangeFieldT& y, const Teuchos::Tuple<bool, 4>& downwinding) const {

    Teuchos::Tuple<OT, 3> ss;
    Teuchos::Tuple<OT, 3> nn;
    Teuchos::Tuple<OT, 3> ii;

    for(F m=F::U; m<GridT::sdim; ++m) {

      for(int i=0; i<3; ++i) {
        ss[i] = (downwinding[i]>0)?(grid()->si(m, i, B::N) ):(grid()->ei(m, i, B::N) );
        nn[i] = (downwinding[i]>0)?(grid()->ei(m, i, B::N)+1):(grid()->si(m, i, B::N)-1);
        ii[i] = (downwinding[i]>0)?1:-1;
      }


      y.getCField()(m).exchange();
      y.getSField()(m).exchange();

      for(OT k=ss[Z]; k!=nn[Z]; k+=ii[Z])
        for(OT j=ss[Y]; j!=nn[Y]; j+=ii[Y])
          for(OT i=ss[X]; i!=nn[X]; i+=ii[X]) {
            if(3==GridT::sdim) {

              ST diag = mulC_*innerOp_->getSOp()->getConvSOp()->innerDiag3D(
                  innerOp_->getConvField(m)[0](i, j, k),
                  innerOp_->getConvField(m)[1](i, j, k),
                  innerOp_->getConvField(m)[2](i, j, k), m, i, j, k)
                - mulL_ * innerOp_->getSOp()->getHelmOp()->innerDiag3D(m, i, j, k) ;

              assert(diag!=0);

              if(downwinding[4]) {
                y.getCField()(m)(i, j, k) += (x.getCField()(m)(i, j, k)
                    - mulI_*y.getSField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc3D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      innerOp_->getConvField(m)[2](i, j, k),
                      y.getCField()(m), i, j, k)
                    +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc3D(
                      y.getCField()(m), m, i, j, k))/diag;
              y.getSField()(m)(i, j, k) += (x.getSField()(m)(i, j, k)
                  + mulI_*y.getCField()(m)(i, j, k)
                  - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc3D(
                    innerOp_->getConvField(m)[0](i, j, k),
                    innerOp_->getConvField(m)[1](i, j, k),
                    innerOp_->getConvField(m)[2](i, j, k),
                    y.getSField()(m), i, j, k)
                  +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc3D(
                    y.getSField()(m), m, i, j, k))/diag;
              }
              else {
                y.getSField()(m)(i, j, k) += (x.getSField()(m)(i, j, k)
                    + mulI_*y.getCField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc3D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      innerOp_->getConvField(m)[2](i, j, k),
                      y.getSField()(m), i, j, k)
                    +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc3D(
                      y.getSField()(m), m, i, j, k))/diag;
                y.getCField()(m)(i, j, k) += (x.getCField()(m)(i, j, k)
                    - mulI_*y.getSField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc3D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      innerOp_->getConvField(m)[2](i, j, k),
                      y.getCField()(m), i, j, k)
                    +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc3D(
                      y.getCField()(m), m, i, j, k))/diag;
              }
            }
            else {
              ST diag = mulC_*innerOp_->getSOp()->getConvSOp()->innerDiag2D(
                  innerOp_->getConvField(m)[0](i, j, k),
                  innerOp_->getConvField(m)[1](i, j, k), m, i, j, k)
                - mulL_ * innerOp_->getSOp()->getHelmOp()->innerDiag2D(m, i, j, k) ;
              assert(diag!=0);

              if(downwinding[4]) {
                y.getCField()(m)(i, j, k) += (x.getCField()(m)(i, j, k)
                    - mulI_*y.getSField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc2D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      y.getCField()(m), i, j, k)
                    +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc2D(
                      y.getCField()(m), m, i, j, k))/diag;
                y.getSField()(m)(i, j, k) += (x.getSField()(m)(i, j, k)
                    + mulI_*y.getCField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc2D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      y.getSField()(m), i, j, k)
                    + mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc2D(
                      y.getSField()(m), m, i, j, k))/diag;
              }
              else{
                y.getSField()(m)(i, j, k) += (x.getSField()(m)(i, j, k)
                    + mulI_*y.getCField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc2D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      y.getSField()(m), i, j, k)
                    + mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc2D(
                      y.getSField()(m), m, i, j, k))/diag;
                y.getCField()(m)(i, j, k) += (x.getCField()(m)(i, j, k)
                    - mulI_*y.getSField()(m)(i, j, k)
                    - mulC_*innerOp_->getSOp()->getConvSOp()->innerStenc2D(
                      innerOp_->getConvField(m)[0](i, j, k),
                      innerOp_->getConvField(m)[1](i, j, k),
                      y.getCField()(m), i, j, k)
                    +mulL_*innerOp_->getSOp()->getHelmOp()->innerStenc2D(
                      y.getCField()(m), m, i, j, k))/diag;
              }
            } // end of if(dimm...
          } // end of loop

      applyBCGS(x.getCField()(m), y.getCField()(m));
      applyBCGS(x.getSField()(m), y.getSField()(m));

      y.getCField()(m).changed();
      y.getSField()(m).changed();
    } // end of F looop
  } // end of applyGS 

}; // end of class ModeSmoother


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MODESMOOTHER_HPP
