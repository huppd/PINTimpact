#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONSOP_HPP
#define PIMPACT_CONVECTIONDIFFUSIONSOP_HPP

#include "Teuchos_ScalarTraits.hpp"

#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_DiffusionOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
template<class ST>
class ConvectionDiffusionSOp {

public:

  using GridT = ST;

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  using FluxFieldT = ScalarField<GridT>[3];
  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

protected:

  Teuchos::RCP<const ConvectionSOp<GridT> > convSOp_;
  Teuchos::RCP<const DiffusionOp<GridT> > helmOp_;

  Scalar mulI_;
  Scalar mulC_;
  Scalar mulL_;

public:

  ConvectionDiffusionSOp(const Teuchos::RCP<const GridT>& grid ):
    convSOp_(create<ConvectionSOp>(grid)),
    helmOp_(create<DiffusionOp>(grid)),
    mulI_(0.),
    mulC_(1.),
    mulL_(1./grid()->getDomainSize()->getRe())	{};


  void assignField(const RangeFieldT& mv) {};


  /// \f[ y =   (wind\cdot\nabla) x - \frac{1}{Re} \Delta x \f]
  void apply(const FluxFieldT& wind, const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N) const {

    //std::cout << "mulI: " << mulI_ << "\n";
    //std::cout << "mulC: " << mulC_ << "\n";
    //std::cout << "mulL: " << mulL_ << "\n";
    apply(wind, x, y, mulI_, mulC_, mulL_, add);
  }


  /// \f[ z = mul z + mulI y + mulC(x\cdot\nabla)y - mulL \Delta y \f]
  void apply(const FluxFieldT& wind, const DomainFieldT& y, RangeFieldT& z,
              const Scalar mulI, const Scalar mulC, const Scalar mulL, const Add add=Add::N) const {

    assert(z.getType() == y.getType());
    for(int i=0; i<GridT::sdim; ++i)
      assert(wind[i].getType() == y.getType());

    for(int vel_dir=0; vel_dir<GridT::sdim; ++vel_dir)
      wind[vel_dir].exchange();


    const F m = y.getType();

    const B wB = ((Add::N==add) ? B::Y : B::N);
    const B wnB = B::N;

    y.exchange();


    //std::cout << "mulI: " << mulI << "\tmulC: " << mulC << "\tmulL: " << mulL << "\tsize: " << wind[0].getLength() << "\tnorm: " << wind[0].norm() << "\n";

    if(3==GridT::sdim) {
      for(Ordinal k=grid()->si(m, Z, wnB); k<=grid()->ei(m, Z, wnB); ++k)
        for(Ordinal j=grid()->si(m, Y, wnB); j<=grid()->ei(m, Y, wnB); ++j)
          for(Ordinal i=grid()->si(m, X, wnB); i<=grid()->ei(m, X, wnB); ++i) {
            if(Add::N==add) z(i, j, k) = 0.;
            z(i, j, k) +=
              + mulI * y(i, j, k)
              + mulC * convSOp_->innerStenc3D(
                wind[0](i, j, k),
                wind[1](i, j, k),
                wind[2](i, j, k), y, i, j, k)
              - mulL * helmOp_->innerStenc3D(y, m, i, j, k);
          }
    } else {
      for(Ordinal k=grid()->si(m, Z, wnB); k<=grid()->ei(m, Z, wnB); ++k)
        for(Ordinal j=grid()->si(m, Y, wnB); j<=grid()->ei(m, Y, wnB); ++j)
          for(Ordinal i=grid()->si(m, X, wnB); i<=grid()->ei(m, X, wnB); ++i) {
            if(Add::N==add) z(i, j, k) = 0.;
            z(i, j, k) +=
              + mulI*y(i, j, k)
              + mulC*convSOp_->innerStenc2D(wind[0](i, j, k), wind[1](i, j, k), y, i, j, k)
              - mulL * helmOp_->innerStenc2D(y, m, i, j, k);
          }
    }

    if(B::Y==wB) helmOp_->applyBC(y, z);

    z.changed();
  }

  void computeResidual(const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res) const {
    apply(x, res);
    res.add(1., b, -1., res);
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return helmOp_->grid();
  };

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {

    if(para->name()!="Linear Solver") {
      mulI_ = para->get<Scalar>("mulI");
      mulC_ = para->get<Scalar>("mulC");
      mulL_ = para->get<Scalar>("mulL");
    }
  }

  constexpr const Teuchos::RCP<const ConvectionSOp<GridT> >& getConvSOp() const {
    return convSOp_;
  }
  constexpr const Teuchos::RCP<const DiffusionOp<GridT> >& getHelmOp() const {
    return helmOp_;
  }

  void print(std::ostream& out=std::cout) const {
    out << "--- " << getLabel() << " ---\n";
    out << "mulI: " << mulI_ << "\n";
    out << "mulC: " << mulC_ << "\n";
    out << "mulL: " << mulL_ << "\n";
    convSOp_->print(out);
    helmOp_->print(out);
  }

  constexpr const Scalar getMulI() const {
    return mulI_;
  }
  constexpr const Scalar getMulC() const {
    return mulC_;
  }
  constexpr const Scalar getMulL() const {
    return mulL_;
  }

  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const std::string getLabel() const {
    return "ConvectionDiffusion";
  };


}; // end of class ConvectionDiffusionSOp


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONSOP_HPP
