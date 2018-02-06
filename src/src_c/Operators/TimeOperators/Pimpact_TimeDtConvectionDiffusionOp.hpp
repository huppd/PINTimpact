#pragma once
#ifndef PIMPACT_TIMEDTCONVECTIONDIFFUSIONOP_HPP
#define PIMPACT_TIMEDTCONVECTIONDIFFUSIONOP_HPP


#include "Pimpact_ConvectionField.hpp"
#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_NonlinearVWrap.hpp"
#include "Pimpact_TimeField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \ingroup TimeHarmonicOperator
/// \ingroup NonliearOperator
/// \deprecated
template<class SpT, int meth=1 >
class TimeDtConvectionDiffusionOp {

public:

  //static const int method = meth;

  using GridT = SpT;

  using DomainFieldT = TimeField<VectorField<GridT> >;
  using RangeFieldT = TimeField<VectorField<GridT> >;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  Teuchos::RCP<NonlinearWrap<ConvectionDiffusionSOp<GridT> > > op_;

  Teuchos::Array<Teuchos::RCP<ConvectionField<GridT> > > wind_;

public:


  TimeDtConvectionDiffusionOp(const Teuchos::RCP<const GridT>& grid):
    op_(create<NonlinearWrap>(create<ConvectionDiffusionSOp<GridT> >(grid))),
    wind_(grid()->nLoc(3) + grid()->bu(3) - grid()->bl(3)) {

    OT nt = grid()->nLoc(3) + grid()->bu(3) - grid()->bl(3);

    for(OT i=0; i<nt; ++i)
      wind_[i] = create<ConvectionField>(grid);
  };

  void assignField(const DomainFieldT& mv) {

    mv.exchange();

    OT nt = grid()->nLoc(3) + grid()->bu(3) - grid()->bl(3);

    for(OT i=0; i<nt; ++i) {
      wind_[i]->assignField(mv(i));
    }
  };


  void apply(const DomainFieldT& y, RangeFieldT& z, bool init_yes=true) const {

    OT sInd = grid()->si(F::S, 3);
    OT eInd = grid()->ei(F::S, 3);

    ST iRe = 1./grid()->getDomainSize()->getRe();
    ST a2 = grid()->getDomainSize()->getAlpha2()*iRe;

    ST pi = 4.*std::atan(1.);

    ST mulI = a2*(static_cast<ST>(grid()->nGlo(3)))/2./pi;


    y.exchange();

    switch(meth) {
    case 0 : {
      for(OT i=sInd; i<=eInd; ++i) {
        op_->apply(wind_[i]->get(), y(i), z(i), mulI, 1., iRe, Add::N);
        z(i).add(1., z(i), -mulI, y(i-1), B::N);
      }
      break;
    }
    case 1: {
      for(OT i=sInd; i<=eInd; ++i) { // explicit looping
        op_->apply(wind_[i  ]->get(), y(i), z(i),  mulI, 0.5, iRe*0.5, Add::N);
        op_->apply(wind_[i-1]->get(), y(i-1), z(i), -mulI, 0.5, iRe*0.5, Add::Y);
      }
      break;
    }
    case 2: {
      for(OT i=sInd; i<=eInd; ++i) {
        op_->apply(wind_[i]->get(), y(i), z(i), 0., 1., iRe, Add::N);
        z(i).add(1., z(i), -mulI/2., y(i-1), B::N);
        z(i).add(1., z(i),  mulI/2., y(i+1), B::N);
      }
      break;
    }
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
    return "TimeDtConvectionDiffusionOp";
  };

  void print(std::ostream& out=std::cout) const {
    out << getLabel() << ":\n";
    op_->print(out);
  }

}; // end of class TimeDtConvectionDiffusionOp



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_TIMEDTCONVECTIONDIFFUSIONOP_HPP
