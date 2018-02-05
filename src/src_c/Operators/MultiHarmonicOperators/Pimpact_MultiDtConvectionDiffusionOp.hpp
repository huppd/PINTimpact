#pragma once
#ifndef PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP
#define PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP


#include "Pimpact_ConvectionField.hpp"
#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_NonlinearVWrap.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class SpT>
class MultiDtConvectionDiffusionOp {

public:

  using GridT = SpT;

  using DomainFieldT = MultiHarmonicField<VectorField<GridT> >;
  using RangeFieldT = MultiHarmonicField<VectorField<GridT> >;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  Teuchos::RCP<NonlinearWrap<ConvectionDiffusionSOp<GridT> > > op_;

  Teuchos::RCP<ConvectionField<GridT> > wind0_;
  Teuchos::Array<Teuchos::RCP<ConvectionField<GridT> > > windc_;
  Teuchos::Array<Teuchos::RCP<ConvectionField<GridT> > > winds_;

  using FieldTensorT = typename ConvectionField<GridT>::FieldTensor;

public:

  MultiDtConvectionDiffusionOp(const Teuchos::RCP<const GridT>& grid):
    op_(create<NonlinearWrap>(create<ConvectionDiffusionSOp<GridT> >(grid))),
    wind0_(create<ConvectionField>(grid)),
    windc_(grid->nGlo(3)),
    winds_(grid->nGlo(3)) {

    for(OT i=0; i<grid->nGlo(3); ++i) {
      windc_[i] = create<ConvectionField>(grid);
      winds_[i] = create<ConvectionField>(grid);
    }
  };


  void assignField(const DomainFieldT& y_ref) {

    Teuchos::RCP<const DomainFieldT> y;

    if(y_ref.global()==DomainFieldT::Global::Y)
      y = Teuchos::rcpFromRef(y_ref);
    else {
      Teuchos::RCP<DomainFieldT> temp =
        Teuchos::rcp(new DomainFieldT(grid(), DomainFieldT::Global::Y));
      *temp = y_ref;
      y = temp;
      //std::cout <<"assign op: y->global(): " <<y->global() <<"\n";
    }

    y->exchange();

    wind0_->assignField(y->get0Field());

    for(OT i=1; i<=grid()->nGlo(3); ++i) {
      windc_[i-1]->assignField(y->getCField(i));
      winds_[i-1]->assignField(y->getSField(i));
    }
  };


  void apply(const DomainFieldT& y_ref, RangeFieldT& z, bool init_yes=true) const {

    Teuchos::RCP<const DomainFieldT> y;
    if(y_ref.global()==DomainFieldT::Global::Y)
      y = Teuchos::rcpFromRef(y_ref);
    else {
      Teuchos::RCP<DomainFieldT> temp =
        Teuchos::rcp(new DomainFieldT(grid(), DomainFieldT::Global::Y));
      *temp = y_ref; // needed because of const
      y = temp;
    }

    y->exchange();

    OT Nf = grid()->nGlo(3);
    ST iRe = 1./op_->grid()->getDomainSize()->getRe();
    ST a2 = op_->grid()->getDomainSize()->getAlpha2()*iRe;

    ST mulI;

    // computing zero mode of z
    if(0==grid()->si(F::U, 3)) {

      op_->apply(get0Wind(), y->get0Field(), z.get0Field(), 0., 1., iRe, Add::N);

      for(OT i=1; i<=Nf; ++i) {
        op_->apply(getCWind(i), y->getCField(i), z.get0Field(), 0., 0.5, 0., Add::Y);
        op_->apply(getSWind(i), y->getSField(i), z.get0Field(), 0., 0.5, 0., Add::Y);
      }
    }

    // computing cos mode of z
    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i) {

      op_->apply(get0Wind(), y->getCField(i), z.getCField(i), 0., 1., iRe, Add::N);
      op_->apply(getCWind(i), y->get0Field(), z.getCField(i), 0., 1., 0.,  Add::Y);

      for(OT k=1; k+i<=Nf; ++k) { // thats fine

        mulI = (k==i)?(a2*i):0;

        op_->apply(getCWind(k+i), y->getCField(k), z.getCField(i),   0., 0.5, 0., Add::Y);
        op_->apply(getCWind(k), y->getCField(k+i), z.getCField(i),   0., 0.5, 0., Add::Y);
        op_->apply(getSWind(k+i), y->getSField(k), z.getCField(i), mulI, 0.5, 0., Add::Y);
        op_->apply(getSWind(k), y->getSField(k+i), z.getCField(i),   0., 0.5, 0., Add::Y);
      }
    }

    // computing sin mode of y
    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i) {

      op_->apply(get0Wind(),  y->getSField(i), z.getSField(i), 0., 1., iRe, Add::N);
      op_->apply(getSWind(i), y->get0Field(),  z.getSField(i), 0., 1., 0. , Add::Y);

      for(OT k=1; k+i<=Nf; ++k) { // that is fine

        mulI = (k==i)?(a2*i):0;

        op_->apply(getCWind(k+i), y->getSField(k), z.getSField(i),    0., -0.5, 0., Add::Y);
        op_->apply(getCWind(k), y->getSField(k+i), z.getSField(i),    0.,  0.5, 0., Add::Y);
        op_->apply(getSWind(k+i), y->getCField(k), z.getSField(i), -mulI,  0.5, 0., Add::Y);
        op_->apply(getSWind(k), y->getCField(k+i), z.getSField(i),    0., -0.5, 0., Add::Y);
      }
    }

    // rest of time
    for(OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i) {
      if(Nf/2+1<=i && i<=Nf) {
        mulI = a2*i;
        z.getCField(i).add(1., z.getCField(i),  mulI, y->getSField(i), B::N);
        z.getSField(i).add(1., z.getSField(i), -mulI, y->getCField(i), B::N);
      }
    }

    // strange terms
    OT i;
    for(OT k=1; k<=Nf; ++k) {
      for(OT l=1; l<=Nf; ++l) { // that is fine
        i = k+l;
        if(i<=Nf) { // do something here
          if(std::max(grid()->si(F::U, 3), 1)<=i && i<=grid()->ei(F::U, 3)) {
            op_->apply(getCWind(k), y->getCField(l), z.getCField(i), 0.,  0.5, 0., Add::Y);
            op_->apply(getSWind(k), y->getSField(l), z.getCField(i), 0., -0.5, 0., Add::Y);

            op_->apply(getCWind(k), y->getSField(l), z.getSField(i), 0.,  0.5, 0., Add::Y);
            op_->apply(getSWind(k), y->getCField(l), z.getSField(i), 0.,  0.5, 0., Add::Y);
          }
        }
      }
    }

    z.changed();
  }


  /// computes residual of higher modes >N_f, and returns the according residual of the
  /// modes that are bigger than the cutoff
  std::pair<ST, OT> compRefRes(const DomainFieldT& y_ref, ST cutoff=1.e-6) const {

    Teuchos::RCP<const DomainFieldT> y;
    if(y_ref.global()==DomainFieldT::Global::Y)
      y = Teuchos::rcpFromRef(y_ref);
    else {
      Teuchos::RCP<DomainFieldT> temp =
        Teuchos::rcp(new DomainFieldT(grid(), DomainFieldT::Global::Y));
      *temp = y_ref; // needed because of const
      y = temp;
    }

    y->exchange();

    OT Nf = grid()->nGlo(3);

    MultiField<ModeField<typename DomainFieldT::InnerFieldT> > z(grid(), Nf);


    // strange terms
    for(OT k=1; k<=Nf; ++k) {
      for(OT l=1; l<=Nf; ++l) { // that is fine
        OT i = k+l;
        if(Nf<i && i<=2*Nf) { // do something here
          op_->apply(getCWind(k), y->getCField(l), z.getField(i-Nf-1).getCField(), 0.,  0.5, 0., Add::Y);
          op_->apply(getSWind(k), y->getSField(l), z.getField(i-Nf-1).getCField(), 0., -0.5, 0., Add::Y);

          op_->apply(getCWind(k), y->getSField(l), z.getField(i-Nf-1).getSField(), 0.,  0.5, 0., Add::Y);
          op_->apply(getSWind(k), y->getCField(l), z.getField(i-Nf-1).getSField(), 0.,  0.5, 0., Add::Y);
        }
      }
    }

    std::vector<ST> norms(Nf);
    z.norm(norms, ENorm::L2);

    OT NfR = Nf-1;
    while(NfR>=0 && norms[NfR]<cutoff)
      --NfR;

    ST res = 0.;
    for(OT i=0; i<=NfR; ++i)
      res += std::pow(norms[i], 2);

    return std::make_pair<ST, OT>(std::sqrt(res), NfR+1);
  }


  void applyBC(const DomainFieldT& x, RangeFieldT& y) const {

    if(0==grid()->si(F::U, 3))
      op_->getSOp()->getHelmholtzOp()->applyBC(x.get0Field(), y.get0Field());

    for(typename GridT::OT i=std::max(grid()->si(F::U, 3), 1); i<=grid()->ei(F::U, 3); ++i) {
      op_->getSOp()->getHelmholtzOp()->applyBC(x.getCField(i), y.getCField(i));
      op_->getSOp()->getHelmholtzOp()->applyBC(x.getSField(i), y.getSField(i));
    }
  }


  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "MHDtConvectionDiffusion";
  };

  void print(std::ostream& out=std::cout) const {
    out << getLabel() <<":\n";
    op_->print(out);
  }

protected:

  constexpr const FieldTensorT& get0Wind() const {
    return wind0_->get();
  }

  constexpr const FieldTensorT& getCWind(const OT i) const {
    return windc_[i-1]->get();
  }
  constexpr const FieldTensorT& getSWind(const OT i) const {
    return winds_[i-1]->get();
  }

}; // end of class MultiDtConvectionDiffusionOp



/// \relates MultiDtConvectionDiffusionOp
template<class GridT>
Teuchos::RCP<MultiDtConvectionDiffusionOp<GridT> >
createMultiDtConvectionDiffusionOp(const Teuchos::RCP<const GridT>& grid) {

  return Teuchos::rcp(new MultiDtConvectionDiffusionOp<GridT>(grid));
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP
