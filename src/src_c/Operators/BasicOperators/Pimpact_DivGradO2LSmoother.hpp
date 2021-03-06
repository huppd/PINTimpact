/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_DIVGRADO2LSMOOTHER_HPP
#define PIMPACT_DIVGRADO2LSMOOTHER_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Pimpact_TeuchosTransfer.hpp"




namespace Pimpact {



/// \brief \f$\omega\f$-Line smoother for second Order DivGradOp.
///
/// \todo make stretch ready
/// \relates DivGradO2Op
/// \ingroup BaseOperator
/// \note work in progress
template<class OperatorT>
class DivGradO2LSmoother {

public:

  using GridT = typename OperatorT::GridT;

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

protected:

  using VectorT = Teuchos::SerialDenseVector<Ordinal, Scalar>;
  using OVectorT = Teuchos::SerialDenseVector<Ordinal, Ordinal>;
  using MatrixT = Teuchos::SerialDenseMatrix<Ordinal, Scalar>;

  int nIter_;
  Scalar omega_;

  Teuchos::Tuple<bool, 3> lineDirection_;

  Teuchos::Tuple<Ordinal, 3> n_;

  bool levelYes_;

  const Teuchos::RCP<const OperatorT> op_;

  Teuchos::Tuple<Teuchos::RCP<VectorT>, 3 > dl_;
  Teuchos::Tuple<Teuchos::RCP<VectorT>, 3 > d_;
  Teuchos::Tuple<Teuchos::RCP<VectorT>, 3 > du_;
  Teuchos::Tuple<Teuchos::RCP<VectorT>, 3 > du2_;

  Teuchos::Tuple<Teuchos::RCP<OVectorT> , 3 > ipiv_;

public:

  /// \brief constructor
  ///
  /// \param[in] op pointer to operator that is smoothed
  /// \param[in] pl  Parameter list of options for the multi grid solver.
  ///   These are the options accepted by the solver manager:
  ///   - "omega" - a \c Scalar damping factor. Default: for 2D 0.8 for 3D 6./7.  /
  ///   - "numIters" - a \c int number of smoothing steps. Default: 4  /
  ///   - "BC smoothing" - a \c int type of BC smoothing 0, 0: Lacbobian, else: direct. Default: 0 /
  ///   - "level" - a \c bool number of smoothing steps. Default: false  /
  DivGradO2LSmoother(
    const Teuchos::RCP<const OperatorT>& op,
    const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList()):
    nIter_(pl->get<int>("numIters", 1)),
    omega_(pl->get<Scalar>("omega", 0.75)),
    levelYes_(pl->get<bool>("level", false)),
    op_(op) {

    lineDirection_[X] = pl->get<bool>("X", false);
    lineDirection_[Y] = pl->get<bool>("Y", true );
    lineDirection_[Z] = pl->get<bool>("Z", false);

    assert(!(!lineDirection_[X] && !lineDirection_[Y] && !lineDirection_[Z]));

    for(int dir=0; dir<3; ++dir) {

      n_[dir] = grid()->ei(F::S, dir) - grid()->si(F::S, dir) + 1 ;

      if(true==lineDirection_[dir]) {

        dl_ [dir] = Teuchos::rcp(new VectorT(n_[dir]-1, true));
        d_  [dir] = Teuchos::rcp(new VectorT(n_[dir]  , true));
        du_ [dir] = Teuchos::rcp(new VectorT(n_[dir]-1, true));
        du2_[dir] = Teuchos::rcp(new VectorT(n_[dir]-2, true));

        ipiv_[dir] = Teuchos::rcp(new OVectorT(n_[dir], true));

        // --- diagonal ---
        for(int j=0; j<GridT::sdim; ++j) {
          if(j==dir)
            (*d_[dir])[grid()->si(F::S, dir)-grid()->si(F::S, dir)] += op_->getC(j, grid()->si(F::S, dir),  0);
          else
            (*d_[dir])[grid()->si(F::S, dir)-grid()->si(F::S, dir)] += 0.1*op_->getC(j, grid()->si(F::S, dir),  0);
        }
        for(Ordinal i=grid()->si(F::S, dir)+1; i<=grid()->ei(F::S, dir)-1; ++i) {
          for(int j=0; j<GridT::sdim; ++j)
            (*d_[dir])[i-grid()->si(F::S, dir)] += op_->getC(j, i,  0);
        }
        for(int j=0; j<GridT::sdim; ++j) {
          if(j==dir)
            (*d_[dir])[grid()->ei(F::S, dir)-grid()->si(F::S, dir)] += op_->getC(j, grid()->ei(F::S, dir), 0);
          else
            (*d_[dir])[grid()->ei(F::S, dir)-grid()->si(F::S, dir)] += 0.1*op_->getC(j, grid()->ei(F::S, dir), 0);
        }

        // --- off diagonal
        for(Ordinal i=grid()->si(F::S, dir); i<grid()->ei(F::S, dir); ++i) {
          (*du_[dir])[i-grid()->si(F::S, dir)] = op_->getC(dir, i+1, +1);
          (*dl_[dir])[i-grid()->si(F::S, dir)] = op_->getC(dir, i  , -1);
        }

        Ordinal lu_factorization_sucess;
        Teuchos::LAPACK<Ordinal, Scalar> lapack;
        lapack.GTTRF(
          n_[dir],
          dl_[dir]->values(),
          d_[dir]->values(),
          du_[dir]->values(),
          du2_[dir]->values(),
          ipiv_[dir]->values(),
          &lu_factorization_sucess);
        assert(!lu_factorization_sucess);
      }
    }
  }


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}(x - A y_k) \f]
  /// \note todo solve multiple problems at once
  void apply(const DomainFieldT& b, RangeFieldT& y, const Belos::ETrans
      trans=Belos::NOTRANS) const {

    DomainFieldT temp(grid());

    for(int i=0; i<nIter_; ++i) {

      for(int dir=2; dir>=0; --dir) { // why?
        //for(int dir=0; dir<3; ++dir) { // why?

        if(true==lineDirection_[dir]) {

          op_->computeResidual(b, y, temp);

          Ordinal i[3];

          int d1 = (dir + 1)%3;
          int d2 = (dir + 2)%3;
          if(d2>d1) std::swap(d2, d1);


          Teuchos::RCP<VectorT> b = Teuchos::rcp(new VectorT(n_[dir], false));

          //for(i[d1]=grid()->si(F::S, d1)+1; i[d1]<=grid()->ei(F::S, d1)-1; ++i[d1])
          //for(i[d2]=grid()->si(F::S, d2)+1; i[d2]<=grid()->ei(F::S, d2)-1; ++i[d2]) {
          for(i[d1]=grid()->si(F::S, d1)+1; i[d1]<=grid()->ei(F::S, d1)-1; ++i[d1])
            for(i[d2]=grid()->si(F::S, d2)+1; i[d2]<=grid()->ei(F::S, d2)-1; ++i[d2]) {

              // transfer
              for(i[dir]=grid()->si(F::S, dir); i[dir]<=grid()->ei(F::S, dir); ++i[dir])
                (*b)[ i[dir]-1 ] = temp(i);

              Ordinal lu_solve_sucess;
              Teuchos::LAPACK<Ordinal, Scalar> lapack;
              lapack.GTTRS(
                'N',
                n_[dir],
                1,
                dl_[dir]->values(),
                d_[dir]->values(),
                du_[dir]->values(),
                du2_[dir]->values(),
                ipiv_[dir]->values(),
                b->values(),
                b->stride(),
                &lu_solve_sucess);

              assert(!lu_solve_sucess);

              // transfer back
              for(i[dir]=grid()->si(F::S, dir); i[dir]<=grid()->ei(F::S, dir); ++i[dir])
                y(i) = y(i) + omega_*(*b)[ i[dir]-1 ];
            }

          y.changed();
        } // end of if(true==lineDirection_[dir])
      } // end of for dir = 2...0
    } // end of for i=0...nIter_

    if(levelYes_)
      y.level();
  }

  void assignField(const DomainFieldT& mv) {};

  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const Teuchos::RCP<const GridT>& grid() const {
    return op_->grid();
  };

  void setParameter(Teuchos::RCP<Teuchos::ParameterList> para) {}

  void print(std::ostream& out=std::cout) const {
    out << "--- " << getLabel() << " ---\n";
    out << "\t numIter: " << nIter_ << "\n";
    op_->print(out);
  }

  const std::string getLabel() const {
    return "DivGradO2LSmoother";
  };


}; // end of class DivGradO2LSmoother


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DIVGRADO2LSMOOTHER_HPP
