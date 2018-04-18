/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_MULTIGRID_HPP
#define PIMPACT_MULTIGRID_HPP


#include "Pimpact_MGFields.hpp"
#include "Pimpact_MGOperators.hpp"
#include "Pimpact_MGTransfers.hpp"
#include "Pimpact_MGSmoothers.hpp"
#include "Pimpact_MGGrids.hpp"




/// \defgroup MG MultiGrid
///
/// Multi Grid



namespace Pimpact {




/// \brief basic multi grid calls
///
/// \tparam MGGridsT grid hierarchy
/// \tparam FieldT field type
/// \tparam TransT transfer operator type
/// \tparam RestrT restriction operator type
/// \tparam InterT interpolation operator type
/// \tparam FOperatorT high order operator type
/// \tparam COperatorT low order operator type
/// \tparam SmootherT smoother type
/// \tparam CGST coarse grid solver type
///
/// \todo constructor with FOperator
/// \ingroup MG
template<
  class MGGridsT,
  template<class> class FieldT,
  template<class, class> class TransT,
  template<class> class RestrT,
  template<class> class InterT,
  template<class> class FOperatorT,
  template<class> class COperatorT,
  template<class> class SmootherT,
  template<class> class CGST >
class MultiGrid {

  using FGridT = typename MGGridsT::FGridT;
  using CGridT = typename MGGridsT::CGridT;


  using MGTransfersT = MGTransfers<MGGridsT, TransT, RestrT, InterT>;

  using MGFieldsT = MGFields<MGGridsT, FieldT>;

  using MGOperatorsT = MGOperators<MGGridsT, FOperatorT, COperatorT>;

  using MGSmoothersT = MGSmoothers<MGOperatorsT, SmootherT>;

public:

  using GridT =  FGridT;

  using DomainFieldT = FieldT<FGridT>;
  using RangeFieldT = FieldT<FGridT>;

  using CGridSolverT = CGST<COperatorT<CGridT> >;

protected:

  bool defectCorrection_;
  int numCycles_;
  int cycleType_;
  int numGrids_;

  const Teuchos::RCP<const MGGridsT> mgGrids_;

  Teuchos::RCP<const MGTransfersT> mgTrans_;

  Teuchos::RCP<const MGOperatorsT> mgOps_;

  Teuchos::RCP<const MGSmoothersT> mgSms_;

  Teuchos::RCP<CGridSolverT> cGridSolver_;

public:

  /// \brief constructor
  ///
  /// \param mgGrids
  /// \param pl  Parameter list of options for the multi grid solver.
  ///   These are the options accepted by the solver manager:
  ///   - "defect correction" - a \c bool specifying if defect correction is
  ///     applied before cycling. Default: true  /
  ///   - "init zero" - a \c bool specifying if defect correction is
  ///     applied before cylcing. Default: false  /
  ///   - "numCycles" - a \c int number of cycles. Default:
  ///     FGridT::dimNC - CGridT::dimNC+1)  /
  ///   - "Smoother" - a \c sublist for smoothers
  ///   - "Coarse Grid Solver" - a \c sublist for coarse grid solver
  /// \todo FGridOperator should be passed here
  MultiGrid(
    const Teuchos::RCP<const MGGridsT>& mgGrids,
    const Teuchos::RCP<FOperatorT<typename MGGridsT::FGridT> >& fOperator,
    const Teuchos::RCP<Teuchos::ParameterList>& pl):
    defectCorrection_(pl->get<bool>("defect correction", true)),
    numCycles_(pl->get<int>("numCycles", FGridT::dimNC-CGridT::dimNC+1)),
    cycleType_(pl->get<int>("cycle type", 0)),
    numGrids_(pl->get<int>("numGrids", -1)),
    mgGrids_(mgGrids),
    mgTrans_(createMGTransfers<TransT, RestrT, InterT>(mgGrids)),
    mgOps_(  createMGOperators<FOperatorT, COperatorT>(mgGrids, fOperator)),
    mgSms_(  createMGSmoothers<SmootherT>(mgOps_, Teuchos::sublist(pl, "Smoother"))),
    //cGridSolver_(mgGrids_->participating(-1)?create<CGridSolverT>(mgOps_->get(-1),
          //Teuchos::sublist(pl, "Coarse Grid Solver")):Teuchos::null) {
    cGridSolver_(mgGrids_->participating(numGrids_)?create<CGridSolverT>(mgOps_->get(numGrids_),
          Teuchos::sublist(pl, "Coarse Grid Solver")):Teuchos::null) {

      //print();
      //assert(!(numGrids_==-1 || numGrids_<mgGrids_->getNGrids()));
      if(numGrids_==-1) {
        numGrids_ = mgGrids_->getNGrids()-1;
        //pl->set<int>("numGrids", numGrids_); 
      }
    }



  void apply(const DomainFieldT& rhs, RangeFieldT& y) const {

    if(0==cycleType_)
      applyVcycle(rhs, y);
    else
      applyVFull(rhs, y);
  }


  /// \brief solves \f$ L y = x \f$
  /// defect correction\f$ \hat{L}u_{k+1} = f-L u_k +\hat{L}u_k \f$ and V-cylce for solving with \f$\hat{L}\f$
  /// \note todo extract smooth/restrict/interpolate method???
  /// \note todo template cycle method
  void applyVcycle(const DomainFieldT& rhs, RangeFieldT& y) const {

    MGFieldsT x(   mgGrids_);
    MGFieldsT temp(mgGrids_);
    MGFieldsT b(   mgGrids_);

    // === no defect correction
    if(!defectCorrection_) {
      mgTrans_->getTransferOp()->apply(y, x.get(0));
      mgTrans_->getTransferOp()->apply(rhs, b.get(0));
    }

    for(int j=0; j<numCycles_; ++j) {

      if(defectCorrection_) {
        // defect correction rhs \hat{f}= b = x - L y
        mgOps_->get()->computeResidual(rhs, y, b.get());

        // transfer init y and \hat{f} to coarsest coarse
        mgTrans_->getTransferOp()->apply(y, x.get(0));
        mgTrans_->getTransferOp()->apply(b.get(), b.get(0));

        // residual temp = \hat(L) y
        mgOps_->get(0)->apply(x.get(0), temp.get(0));
        // b = x - L y + \hat{L} y
        b.get(0).add(1., b.get(0), 1, temp.get(0));
      }

      int i;
      for(i=0; i<=numGrids_-1; ++i) {
        if(i>0) x.get(i).init(0.); // necessary? for DivGradOp yes

        if(mgGrids_->participating(i)) {
          mgSms_->get(i)->apply(b.get(i), x.get(i));
          mgOps_->get(i)->computeResidual(b.get(i), x.get(i), temp.get(i));
          mgTrans_->getRestrictionOp(i)->apply(temp.get(i), b.get(i+1));
        }
      }

      // coarse grid solution
      i = numGrids_;
      //i = -1;
      if(mgGrids_->participating(i)) {
        x.get(i).init(0.);
        try {
          cGridSolver_->apply(b.get(i), x.get(i));
        } catch(std::logic_error& e) {
          std::cout << "error in MG on coarse grid:\n";
          cGridSolver_->print();
          b.get(i).print();
          b.get(i).write(111);
          x.get(i).write(222);
          throw(e);
        }
        //x.get(i).level();
      }

      for(i=numGrids_-1; i>=0; --i) {
        // interpolate/correct/smooth
        if(mgGrids_->participating(i)) {
          mgTrans_->getInterpolationOp(i)->apply(x.get(i+1), temp.get(i));
          x.get(i).add(1., temp.get(i), 1., x.get(i));
          //
          mgSms_->get(i)->apply(b.get(i), x.get(i));
        }
      }
      // use temp as stopping criterion
      if(defectCorrection_)
        mgTrans_->getTransferOp()->apply(x.get(0), y);
    }
    if(!defectCorrection_)
      mgTrans_->getTransferOp()->apply(x.get(0), y);
  }

  /// \brief solves \f$ L y = x \f$
  /// defect correction\f$ \hat{L}u_{k+1} = f-L u_k +\hat{L}u_k \f$ and V-cylce for solving with \f$\hat{L}\f$
  void applyVFull(const DomainFieldT& rhs, RangeFieldT& y) const {

    MGFieldsT x(   mgGrids_);
    MGFieldsT temp(mgGrids_);
    MGFieldsT b(   mgGrids_);


    b.get() = rhs;
    mgTrans_->restriction(b);

    // coarse grid solution
    int i = -1;
    if(mgGrids_->participating(i)) {
      try {
        cGridSolver_->apply(b.get(i), x.get(i));
      } catch(std::logic_error& e) {
        std::cout << "error in MG on coarse grid:\n";
        cGridSolver_->print();
        b.get(i).print();
        b.get(i).write(111);
        x.get(i).write(222);
        throw(e);
      }
    }


    for(int start=-2; start>=-mgGrids_->getNGrids(); --start) {

      mgTrans_->getInterpolationOp(start)->apply(x.get(start+1), x.get(start));

      for(int j=0; j<numCycles_; ++j) {

        if(-mgGrids_->getNGrids()==start && defectCorrection_) {
          mgTrans_->getTransferOp()->apply(x.get(0), y);
          // defect correction rhs \hat{f}= b = x - L y
          mgOps_->get()->computeResidual(rhs, y, b.get());

          // transfer init y and \hat{f} to coarsest coarse
          mgTrans_->getTransferOp()->apply(y, x.get(0));
          mgTrans_->getTransferOp()->apply(b.get(), b.get(0));

          // residual temp = \hat(L) y
          mgOps_->get(0)->apply(x.get(0), temp.get(0));
          // b = x - L y + \hat{L} y
          b.get(0).add(1., b.get(0), 1, temp.get(0));
        }

        int i;
        for(i=mgGrids_->getNGrids()+start; i<mgGrids_->getNGrids()-1; ++i) {

          if(mgGrids_->participating(i)) {
            mgSms_->get(i)->apply(b.get(i), x.get(i));
            mgOps_->get(i)->computeResidual(b.get(i), x.get(i), temp.get(i));
            mgTrans_->getRestrictionOp(i)->apply(temp.get(i), b.get(i+1));
          }
        }

        // coarse grid solution
        i = -1;
        if(mgGrids_->participating(i)) {
          x.get(i).init(0.);
          try {
            cGridSolver_->apply(b.get(i), x.get(i));
          } catch(std::logic_error& e) {
            std::cout << "error in MG on coarse grid:\n";
            cGridSolver_->print();
            b.get(i).print();
            b.get(i).write(111);
            x.get(i).write(222);
            throw(e);
          }
          //x.get(i).level();
        }

        for(i=-2; i>=start; --i) {
          // interpolate/correct/smooth
          if(mgGrids_->participating(i)) {
            mgTrans_->getInterpolationOp(i)->apply(x.get(i+1), temp.get(i));
            x.get(i).add(1., temp.get(i), 1., x.get(i));
            //
            mgSms_->get(i)->apply(b.get(i), x.get(i));
          }
        }
      } // end numCycles
    }

    mgTrans_->getTransferOp()->apply(x.get(0), y);
  }

  /// \todo smoother that have to be updated should be "assigned" as well
  void assignField(const DomainFieldT& mv) {
    //std::cout << getLabel() << "assignField\n";

    MGFieldsT temp(mgGrids_);

    //mgOps_->get()->assignField(mv);
    mgTrans_->getTransferOp()->apply(mv, temp.get(0));

    for(int i=0; i<mgGrids_->getNGrids()-1; ++i)  {
      if(mgGrids_->participating(i)) {
        mgOps_->get(i)->assignField(temp.get(i));
        mgTrans_->getRestrictionOp(i)->apply(temp.get(i), temp.get(i+1));
        //std::cout << temp.get(i).norm()<< "\n";
      }
    }

    if(mgGrids_->participating(-1))
      mgOps_->get(-1)->assignField(temp.get(-1));
  };

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return mgGrids_->get();
  };

  void setParameter(const Teuchos::RCP<Teuchos::ParameterList>& para) {
    mgOps_->setParameter(para);
  }

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "MultiGrid("+mgOps_->get()->getLabel()+") ";
  };


  Teuchos::RCP<const MGTransfersT>
  getTransfers() {
    return mgTrans_;
  }

  void print(std::ostream& out=std::cout) const {
    out << "--- " << getLabel() << " ---\n";
    out << "#grids: " << numGrids_ << " numCycles: " << numCycles_ << "\n";
    out << "FOperator: " << mgOps_->get()->getLabel() << " d" << FGridT::dimNC << "\n";
    if(mgGrids_->participating(0))
      out << "COperator: " << mgOps_->get(0)->getLabel()  << " d" << CGridT::dimNC << "\n";
    if(mgGrids_->participating(0))
      out << "Smoother: " << mgSms_->get(0)->getLabel() << "\n";
    if(mgGrids_->participating(-1))
      out << "Coarse Grid Solver: " << cGridSolver_->getLabel() << "\n";
  }


}; // end of class MultiGrid



/// \relates MultiGrid
template<
  template<class> class FieldT,
  template<class, class> class TransT,
  template<class> class RestrT,
  template<class> class InterT,
  template<class> class FOperatorT,
  template<class> class COperatorT,
  template<class> class SmootherT,
  template<class> class CGridSolverT,
  class MGGridsT >
Teuchos::RCP<MultiGrid<MGGridsT, FieldT, TransT, RestrT, InterT, FOperatorT, COperatorT, SmootherT, CGridSolverT> >
createMultiGrid(
  const Teuchos::RCP<const MGGridsT>& mgGrids,
  const Teuchos::RCP<FOperatorT<typename MGGridsT::FGridT> > & fOperatr,
  const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList()) {

  return Teuchos::rcp(new
      MultiGrid<MGGridsT, FieldT, TransT, RestrT, InterT, FOperatorT, COperatorT, SmootherT, CGridSolverT>(
        mgGrids, fOperatr, pl));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIGRID_HPP
