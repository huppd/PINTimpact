/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_NONLINEARVOP_HPP
#define PIMPACT_NONLINEARVOP_HPP


#include "Pimpact_ConvectionField.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_NonlinearVWrap.hpp"




namespace Pimpact {



/// \brief Convection Operator for Velocity fields
/// \ingroup BaseOperator
/// \ingroup NonlinearOperator
template<class CSOPT>
class NonlinearOp {

public:

  using ConvSOpT = CSOPT;

  using GridT = typename ConvSOpT::GridT;

  using DomainFieldT = VectorField<GridT>;
  using RangeFieldT = VectorField<GridT>;

protected:

  using Scalar = typename GridT::Scalar;

  Teuchos::RCP<NonlinearWrap<ConvSOpT> > convVWrap_;

  Teuchos::RCP<ConvectionField<GridT> > convField_;

public:

  NonlinearOp(const Teuchos::RCP<const GridT>& grid):
    convVWrap_(create<NonlinearWrap<ConvSOpT> >(create<ConvSOpT>(grid))),
    convField_(create<ConvectionField>(grid)) {};

  template<class ConvSOpTT >
  NonlinearOp(const Teuchos::RCP<ConvSOpTT>& op):
    convVWrap_(create<NonlinearWrap<ConvSOpT> >(create<ConvSOpT>(op->grid()))),
    convField_(op->getConvField()) {}

  void assignField(const DomainFieldT& mv) const {
    convField_->assignField(mv);
  };

  /// \note Operator's wind has to be assigned correctly
  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    convVWrap_->apply(convField_->get(), x, y);
  }


  void computeResidual(const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res) const {
    apply(x, res);
    res.add(1., b, -1., res);
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return convVWrap_->grid();
  }

  Teuchos::RCP<ConvectionField<GridT> >
  getConvField() const {
    return convField_;
  }

  const ScalarField<GridT>* getConvField(F f) const {
    return convField_->get()[static_cast<int>(f)];
  }

  Teuchos::RCP<NonlinearWrap<ConvSOpT> >
  getOp() const {
    return convVWrap_;
  }

  Teuchos::RCP<const ConvSOpT>
  getSOp() const {
    return convVWrap_->getSOp();
  }


  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    convVWrap_->setParameter(para);
  }


  bool hasApplyTranspose() const {
    return false;
  }

  void print(std::ostream& out=std::cout) const {
    convVWrap_->print(out);
  }

  const std::string getLabel() const {
    return convVWrap_->getSOp()->getLabel() + "VOp";
  };

}; // end of class NonlinearOp



/// \relates NonlinearOp
template<class GridT>
Teuchos::RCP<NonlinearOp<NonlinearWrap<ConvectionSOp<GridT> > > >
createNonlinearOp(
    const Teuchos::RCP<const GridT>& grid) {

  Teuchos::RCP<NonlinearOp<NonlinearWrap<ConvectionSOp<GridT> > > > sop =
    Pimpact::create<Pimpact::ConvectionSOp>(grid) ;

  return Teuchos::rcp(new NonlinearOp<ConvectionSOp<GridT> >(sop));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_NONLINEARVOP_HPP
