#pragma once
#ifndef PIMPACT_TIMEDTCONVECTIONDIFFUSIONPREC_HPP
#define PIMPACT_TIMEDTCONVECTIONDIFFUSIONPREC_HPP


#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_NonlinearVWrap.hpp"
#include "Pimpact_TimeField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \ingroup TimeHarmonicOperator
/// \ingroup NonliearOperator
/// \deprecated
template<class ConvDiffInv, int meth=1>
class TimeDtConvectionDiffusionPrec {

public:

  using SpaceT = typename ConvDiffInv::SpaceT;

  using DomainFieldT = TimeField<VectorField<SpaceT> >;
  using RangeFieldT = TimeField<VectorField<SpaceT> >;

protected:

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  Teuchos::RCP<ConvDiffInv> op_;

  //Teuchos::Array<Teuchos::RCP<ConvectionField<SpaceT> > > wind_;
  TimeField<VectorField<SpaceT> > wind_;

public:


  TimeDtConvectionDiffusionPrec(const Teuchos::RCP<ConvDiffInv>& op):
    op_(op),
    wind_(op->space()) {};

  void assignField(const DomainFieldT& mv) {

    wind_ = mv;
    //mv.exchange();

    //OT nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

    //for(OT i=0; i<nt; ++i) {
    //wind_[i]->assignField(mv(i));
    //}
  };


  void apply(const DomainFieldT& x, RangeFieldT& y, bool init_yes=true) const {

    OT sInd = space()->si(F::S, 3);
    OT eInd = space()->ei(F::S, 3);

    ST iRe = 1./space()->getDomainSize()->getRe();
    ST a2 = space()->getDomainSize()->getAlpha2()*iRe;

    ST pi = 4.*std::atan(1.);

    ST mulI = a2*(static_cast<ST>(space()->nGlo(3)))/2./pi;

    auto para = Teuchos::parameterList();
    switch(meth) {
    case 0 : {
      para->set<ST>("mulI", mulI);
      para->set<ST>("mulC", 1.);
      para->set<ST>("mulL", iRe);
      break;
    }
    case 1 : {
      para->set<ST>("mulI", mulI);
      para->set<ST>("mulC", 0.5);
      para->set<ST>("mulL", iRe/2.);
      break;
    }
    case 2 : {
      para->set<ST>("mulI", 0);
      para->set<ST>("mulC", 1.);
      para->set<ST>("mulL", iRe);
      break;
    }
    }
    op_->setParameter(para);


    for(OT i=sInd; i<=eInd; ++i) {
      op_->assignField(wind_(i));
      op_->apply(x(i), y(i));
    }
    y.changed();
  }



  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op_->space();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {
    op_->setParameter(para);
  }

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "TimeDtConvectionDiffusionPrec";
  };

  void print(std::ostream& out=std::cout) const {
    out <<getLabel() <<":\n";
    op_->print(out);
  }

}; // end of class TimeDtConvectionDiffusionPrec


/// \relates MultiHarmonicDiagOp
template<class OpT>
Teuchos::RCP<TimeDtConvectionDiffusionPrec<OpT> >
createTimeDtConvectionDiffusionPrec(const Teuchos::RCP<OpT>& zeroOp) {

  return Teuchos::rcp(new TimeDtConvectionDiffusionPrec<OpT>(zeroOp));
}


} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_TIMEDTCONVECTIONDIFFUSIONPREC_HPP
