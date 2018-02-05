#pragma once
#ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
#define PIMPACT_DIVGRADO2JSMOOTHER_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_TeuchosTransfer.hpp"




namespace Pimpact {



/// \brief \f$\omega\f$-Jacobian smoother for second Order DivGradOp.
///
///
/// \relates DivGradO2Op
/// \ingroup BaseOperator
/// \note todo instead of hardcode 2nd Order it would be pretty to use new grid with StencilWidth<3, 2>
template<class OperatorT>
class DivGradO2JSmoother {

public:

  using GridT = typename OperatorT::GridT;

  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using SolverT = TeuchosSolver<OperatorT>;

  ST omega_;
  int nIter_;

  bool levelYes_;

  bool jacobi_;

  const Teuchos::RCP<const OperatorT> op_;

public:

  /// \brief constructor
  ///
  /// \param[in] op pointer to operator that is smoothed
  /// \param[in] pl  Parameter list of options for the multi grid solver.
  ///   These are the options accepted by the solver manager:
  ///   - "omega" - a \c ST damping factor. Default: for 2D 0.8 for 3D 6./7.  /
  ///   - "numIters" - a \c int number of smoothing steps. Default: 4  /
  ///   - "BC smoothing" - a \c int type of BC smoothing 0, 0: Jacbobian, else: direct. Default: 0 /
  ///   - "level" - a \c bool if pressure is leveled. Default: false  /
  ///   - "Jacobi" - a \c bool if Jacbonian is used or GausSeidel. Default: true
  DivGradO2JSmoother(
      const Teuchos::RCP<const OperatorT>& op,
      const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList()):
    omega_(pl->get<ST>("omega", (2==GridT::sdim)?0.8:6./7.)),
    nIter_(pl->get<int>("numIters", 2)),
    levelYes_(pl->get<bool>("level", false)),
    jacobi_(pl->get<bool>("Jacobi", true)),
    op_(op) {

      if(!jacobi_)
        omega_ = 1.;
    }


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}(x - A y_k) \f]
  void apply(const DomainFieldT& b, RangeFieldT& y, const Add add=Add::N) const {
    if(jacobi_)
      applyJ(b, y, add);
    else
      applyGS(b, y, add);
  }

protected:

  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}(x - A y_k) \f]
  void applyJ(const DomainFieldT& b, RangeFieldT& y, const Add add=Add::N) const {

    DomainFieldT temp(grid());

    for(int i=0; i<nIter_; ++i) {

      y.exchange();

      if(3==GridT::sdim)
        for(OT k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
          for(OT j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
            for(OT i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
              temp(i, j, k) = innerStenc3D(b, y, i, j, k);
            }
      else
        for(OT k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
          for(OT j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
            for(OT i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
              temp(i, j, k) = innerStenc2D(b, y, i, j, k);
            }

      temp.changed();
      temp.exchange();

      if(3==GridT::sdim)
        for(OT k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
          for(OT j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
            for(OT i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
              y(i, j, k) = innerStenc3D(b, temp, i, j, k);
            }
      else
        for(OT k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
          for(OT j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
            for(OT i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
              y(i, j, k) = innerStenc2D(b, temp, i, j, k);
            }

      y.changed();
    }
    if(levelYes_)
      y.level();
  }


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}(x - A y_k) \f]
  void applyGS(const DomainFieldT& b, RangeFieldT& y, const Add add=Add::N) const {

    for(int i=0; i<nIter_; ++i) {

      y.exchange();

      if(3==GridT::sdim)
        for(OT k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
          for(OT j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
            for(OT i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
              y(i, j, k) = innerStenc3D(b, y, i, j, k);
            }
      else
        for(OT k=grid()->si(F::S, Z); k<=grid()->ei(F::S, Z); ++k)
          for(OT j=grid()->si(F::S, Y); j<=grid()->ei(F::S, Y); ++j)
            for(OT i=grid()->si(F::S, X); i<=grid()->ei(F::S, X); ++i) {
              y(i, j, k) = innerStenc2D(b, y, i, j, k);
            }

      y.changed();
    }
    if(levelYes_)
      y.level();
  }

public:

  void assignField(const DomainFieldT& mv) {};

  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  void print(std::ostream& out=std::cout) const {
    out <<"--- " <<getLabel() <<" ---\n";
    out <<"\t omega: " <<omega_ <<"\n";
    out <<"\t numIter: " <<nIter_ <<"\n";
    op_->print(out);
  }

  const std::string getLabel() const {
    return "DivGradO2JSmoother";
  };

protected:

  constexpr ST innerStenc3D(const DomainFieldT& b, const DomainFieldT& x,
                             const OT i, const OT j, const OT k) const {

    const bool bcX = (grid()->getBCLocal()->getBCL(X) > 0 && i==grid()->si(F::S, X)) ||
                     (grid()->getBCLocal()->getBCU(X) > 0 && i==grid()->ei(F::S, X)) ;
    const bool bcY = (grid()->getBCLocal()->getBCL(Y) > 0 && j==grid()->si(F::S, Y)) ||
                     (grid()->getBCLocal()->getBCU(Y) > 0 && j==grid()->ei(F::S, Y)) ;
    const bool bcZ = (grid()->getBCLocal()->getBCL(Z) > 0 && k==grid()->si(F::S, Z)) ||
                     (grid()->getBCLocal()->getBCU(Z) > 0 && k==grid()->ei(F::S, Z)) ;

    //const ST eps = 0.1;
    const ST eps = 1./static_cast<ST>(OperatorT::epsI);

    const ST epsX = (bcY||bcZ)?eps:1.;
    const ST epsY = (bcX||bcZ)?eps:1.;
    const ST epsZ = (bcX||bcY)?eps:1.;

    ST diag = epsX*getC(X, i, 0) + epsY*getC(Y, j, 0) + epsZ*getC(Z, k, 0);

    diag = (diag==0.)?1.:diag;
    return (1.-omega_)*x(i, j, k) +
      omega_/ diag *(
          b(i, j, k) -
          epsX*getC(X, i, -1)*x(i-1, j  , k) - epsX*getC(X, i, 1)*x(i+1, j  , k) -
          epsY*getC(Y, j, -1)*x(i  , j-1, k) - epsY*getC(Y, j, 1)*x(i  , j+1, k) -
          epsZ*getC(Z, k, -1)*x(i  , j  , k-1) - epsZ*getC(Z, k, 1)*x(i  , j  , k+1));
  }

  constexpr ST innerStenc2D(const DomainFieldT& b, const DomainFieldT& x, const OT i,
      const OT j, const OT k) const {

    const bool bcX = (grid()->getBCLocal()->getBCL(X) > 0 && i==grid()->si(F::S, X)) ||
                     (grid()->getBCLocal()->getBCU(X) > 0 && i==grid()->ei(F::S, X)) ;
    const bool bcY = (grid()->getBCLocal()->getBCL(Y) > 0 && j==grid()->si(F::S, Y)) ||
                     (grid()->getBCLocal()->getBCU(Y) > 0 && j==grid()->ei(F::S, Y)) ;

    //const ST eps = 0.1;
    const ST eps = 1./static_cast<ST>(OperatorT::epsI);

    const ST epsX = bcY?eps:1.;
    const ST epsY = bcX?eps:1.;

    ST diag = (epsX*getC(X, i, 0) + epsY*getC(Y, j, 0));
    diag = (diag!=0.)?diag:1.;
    return (1.-omega_)*x(i, j, k) +
      omega_/diag*(b(i, j, k) -
          epsX*getC(X, i, -1)*x(i-1, j  , k) - epsX*getC(X, i, 1)*x(i+1, j  , k) -
          epsY*getC(Y, j, -1)*x(i  , j-1, k) - epsY*getC(Y, j, 1)*x(i  , j+1, k));
  }

  constexpr const ST* getC(const ECoord dir) const  {
    return op_->getC(dir);
  }

  constexpr const ST* getC(const int dir) const  {
    return op_->getC(dir);
  }

  constexpr const ST getC(const ECoord dir, OT i, OT off) const  {
    return op_->getC(dir, i, off);
  }

  constexpr const ST getC(const int dir, OT i, OT off) const  {
    return op_->getC(dir, i, off);
  }

}; // end of class DivGradO2JSmoother



template<template<class> class SmootherT, class OperatorT>
Teuchos::RCP<SmootherT<OperatorT> >
create(
  const Teuchos::RCP<OperatorT>& op,
  const Teuchos::RCP<Teuchos::ParameterList>& pl) {

  return Teuchos::rcp(new SmootherT<OperatorT>(op, pl));
}


/// \todo move somewhere better
template<class SmootherT, class OperatorT>
Teuchos::RCP<SmootherT >
create(
  const Teuchos::RCP<OperatorT>& op,
  const Teuchos::RCP<Teuchos::ParameterList>& pl) {

  return Teuchos::rcp(new SmootherT(op, pl));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
