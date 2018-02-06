#pragma once
#ifndef PIMPACT_DIVGRADO2INV_HPP
#define PIMPACT_DIVGRADO2INV_HPP


#include "Teuchos_LAPACK.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_TeuchosTransfer.hpp"




namespace Pimpact {




/// \brief inverse for second Order DivGradOp.
///
/// \todo make the same for ConvectionDiffusionSOp, or better more general
/// \relates DivGradO2Op
/// \ingroup BaseOperator
template<class OperatorT>
class DivGradO2Inv {

public:

  using GridT = typename OperatorT::GridT;

  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

protected:

  bool levelYes_;

  const TeuchosSolver<OperatorT> solver_;

public:

  /// \brief constructor
  ///
  /// \param[in] op pointer to operator that is smoothed
  /// \param[in] pl  Parameter list of options for the multi grid solver.
  ///   These are the options accepted by the solver manager: none
  DivGradO2Inv(const Teuchos::RCP<const OperatorT>& op,
                const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList()):
    levelYes_(pl->get<bool>("level", false)),
    solver_(op) { }


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}(x - N y_k) \f]
  void apply(const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N) const {

    solver_.apply(x, y);

    if(levelYes_)
      y.level();
  }

  void assignField(const DomainFieldT& mv) {};

  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return solver_.getOperator()->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  void print(std::ostream& out=std::cout) const {
    out << "--- " << getLabel() << " ---\n";
    solver_.getOperator()->print(out);
    //out << "\n" << *A_ << "\n";
  }

  const std::string getLabel() const {
    return "DivGradO2Inv";
  };

}; // end of class DivGradO2Inv



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DIVGRADO2INV_HPP
