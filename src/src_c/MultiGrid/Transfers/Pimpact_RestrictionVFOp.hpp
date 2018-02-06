#pragma once
#ifndef PIMPACT_RESTRICTIONVFOP_HPP
#define PIMPACT_RESTRICTIONVFOP_HPP
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Pimpact_RestrictionBaseOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Grid.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_InterpolateS2VOp.hpp"




namespace Pimpact {




/// \brief Opetartor that restricts from a fine grid to a coarse grid
///
/// \tparam GT type of the \c Grid
template<class GT>
class RestrictionVFOp : private RestrictionBaseOp<GT> {

  static const int dimension = GT::dimension;

  using Scalar = typename GT::Scalar;
  using Ordinal = typename GT::Ordinal;

public:

  using GridT = GT;

  using FGridT = GridT;
  using CGridT = GridT;

  using DomainFieldT = ScalarField<GridT>;
  using RangeFieldT = ScalarField<GridT>;

  using StencS = Stencil<Scalar, Ordinal, 1, -1, 1 >;
  using StencV = Stencil<Scalar, Ordinal, 0,  0, 1 >;

protected:

  Teuchos::Tuple<StencS, 3> cRS_;
  Teuchos::Tuple<StencV, 3> cRV_;

  Ordinal getIF(const int dir, const Ordinal ii) const {

    Ordinal i=this->dd_[dir]*(ii-1) + 1;

    if(0<gridF()->getBCLocal()->getBCL(dir))
      i = std::max(0, i);
    if(0<gridF()->getBCLocal()->getBCU(dir))
      i = std::min(gridF()->ei(F::S, dir)-1, i);
    return i;
  }

  /// \todo mv MG_getCRVS to Base class deal with BC here
  void initVF() {

    // ------------------------- CRS, CRV
    for(int dir=0; dir<3; ++dir) {

      const Ordinal iimax = this->iimax_[dir];

      cRS_[dir] = StencS(iimax);

      MG_getCRVS(
        this->iimax_[dir],
        (this->nGather_[dir]>1)?
        gridF()->getBCLocal()->getBCL(dir):
        gridC()->getBCLocal()->getBCL(dir),
        (this->nGather_[dir]>1)?
        gridF()->getBCLocal()->getBCU(dir):
        gridC()->getBCLocal()->getBCU(dir),
        this->dd_[dir],
        gridF()->getGridSizeLocal()->get(dir),
        gridF()->bl(dir),
        gridF()->bu(dir),
        gridF()->getCoordinatesLocal()->getX(F::S, dir),
        cRS_[dir].get());


      cRV_[dir] = StencV(iimax);

      const auto& xf = gridF()->getCoordinatesLocal()->getV(dir);
      const auto& xs = gridF()->getCoordinatesLocal()->getS(dir);

      // Restriktion, linienweise, 1d
      // fine
      //     xf(i)        xs(i+1)        xf(i+1)
      //  ----->-------------o------------->-----
      //       |-----------Dx12------------|
      //  ------------------->-------------------
      //                  xc(ii)
      //
      // coarse
      for(Ordinal ii=0; ii<=iimax; ++ii) {

        Ordinal i = getIF(dir, ii);

        Scalar dx12 = xf[i+1] - xf[i];

        cRV_[dir](ii, 0) = (xf[i+1]-xs[i+1])/dx12;
        cRV_[dir](ii, 1) = (xs[i+1]-xf[i ])/dx12;
      }

      // Dirichlet boundary conditions
      if(0<gridF()->getBCLocal()->getBCL(dir)) {
        cRV_[dir](0, 0) = 1.;
        cRV_[dir](0, 1) = 0.;
      }
      if(0<gridF()->getBCLocal()->getBCU(dir)) {
        cRV_[dir](iimax, 0) = 0.;
        cRV_[dir](iimax, 1) = 1.;
      }

      // symmetric boundary conditions
      if(BC::Symmetry==gridC()->getBCLocal()->getBCL(dir)) {
        cRV_[dir](0, 0) = 0.;
        cRV_[dir](0, 1) = 0.;
      }

      if(BC::Symmetry==gridC()->getBCLocal()->getBCU(dir)) {
        cRV_[dir](iimax, 0) = 0.;
        cRV_[dir](iimax, 1) = 0.;
      }
    }
  }

public:

  RestrictionVFOp(
    const Teuchos::RCP<const GridT>& gridF,
    const Teuchos::RCP<const GridT>& gridC):
    RestrictionBaseOp<GT>(gridF, gridC) {

    initVF();
  }


  RestrictionVFOp(
    const Teuchos::RCP<const GridT>& gridF,
    const Teuchos::RCP<const GridT>& gridC,
    const Teuchos::Tuple<int, dimension>& np):
    RestrictionBaseOp<GT>(gridF, gridC, np) {

    initVF();
  }



  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    assert(x.getType()==y.getType());
    assert(x.getType()!=F::S);

    F fType  = x.getType();
    x.exchange();

    switch(fType) {
      case(F::U) : {
        if(2==GT::sdim) {
          for(Ordinal kk=gridC()->si(fType, Z, B::Y); kk<=this->iimax_[Z]; ++kk) {
            Ordinal k = kk;
            for(Ordinal jj=gridC()->si(fType, Y, B::Y); jj<=this->iimax_[Y]; ++jj) {
              Ordinal j = this->dd_[Y]*(jj - 1) + 1;
              for(Ordinal ii=gridC()->si(fType, X, B::Y); ii<=this->iimax_[X]; ++ii) {
                Ordinal i = getIF(X, ii);

                y(ii, jj, kk) = 0.;

                for(int jjj=-1; jjj<=1; ++jjj)
                  for(int iii=0; iii<=1; ++iii)
                    y(ii, jj, kk) +=
                      cRV_[X](ii, iii)*cRS_[Y](jj, jjj)*x(i+iii, j+jjj, k) ;
              }
            }
          }
        } else {
          for(Ordinal kk=gridC()->si(fType, Z, B::Y); kk<=this->iimax_[Z]; ++kk) {
            Ordinal k = this->dd_[Z]* (kk - 1) + 1;
            for(Ordinal jj=gridC()->si(fType, Y, B::Y); jj<=this->iimax_[Y]; ++jj) {
              Ordinal j = this->dd_[Y]*(jj - 1) + 1;
              for(Ordinal ii=gridC()->si(fType, X, B::Y); ii<=this->iimax_[X]; ++ii) {
                Ordinal i = getIF(X, ii);

                y(ii, jj, kk) = 0.;

                for(int kkk=-1; kkk<=1; ++kkk)
                  for(int jjj=-1; jjj<=1; ++jjj)
                    for(int iii=0; iii<=1; ++iii)
                      y(ii, jj, kk) +=
                        cRV_[X](ii, iii)*cRS_[Y](jj, jjj)*cRS_[Z](kk, kkk)*x(i+iii, j+jjj, k+kkk) ;
              }
            }
          }
        }
        break;
      }
      case(F::V) : {

        if(2==GT::sdim) {
          for(Ordinal kk=gridC()->si(fType, Z, B::Y); kk<=this->iimax_[Z]; ++kk) {
            Ordinal k = kk;
            for(Ordinal jj=gridC()->si(fType, Y, B::Y); jj<=this->iimax_[Y]; ++jj) {
              Ordinal j = getIF(Y, jj);
              for(Ordinal ii=gridC()->si(fType, X, B::Y); ii<=this->iimax_[X]; ++ii) {
                Ordinal i = this->dd_[X]*(ii - 1) + 1;

                y(ii, jj, kk) = 0.;

                for(int jjj=0; jjj<=1; ++jjj)
                  for(int iii=-1; iii<=1; ++iii)
                    y(ii, jj, kk) +=
                      cRS_[X](ii, iii)*cRV_[Y](jj, jjj)*x(i+iii, j+jjj, k) ;
              }
            }
          }
        } else {
          for(Ordinal kk=gridC()->si(fType, Z, B::Y); kk<=this->iimax_[Z]; ++kk) {
            Ordinal k = this->dd_[Z]* (kk - 1) + 1;
            for(Ordinal jj=gridC()->si(fType, Y, B::Y); jj<=this->iimax_[Y]; ++jj) {
              Ordinal j = getIF(Y, jj);
              for(Ordinal ii=gridC()->si(fType, X, B::Y); ii<=this->iimax_[X]; ++ii) {
                Ordinal i = this->dd_[X]*(ii - 1) + 1;

                y(ii, jj, kk) = 0.;

                for(int kkk=-1; kkk<=1; ++kkk)
                  for(int jjj=0; jjj<=1; ++jjj)
                    for(int iii=-1; iii<=1; ++iii)
                      y(ii, jj, kk) +=
                        cRS_[X](ii, iii)*cRV_[Y](jj, jjj)*cRS_[Z](kk, kkk)*x(i+iii, j+jjj, k+kkk) ;
              }
            }
          }
        }
        break;
      }
      case(F::W) : {

        for(Ordinal kk=gridC()->si(fType, Z, B::Y); kk<=this->iimax_[Z]; ++kk) {
          Ordinal k = getIF(Z, kk);
          for(Ordinal jj=gridC()->si(fType, Y, B::Y); jj<=this->iimax_[Y]; ++jj) {
            Ordinal j = this->dd_[Y]*(jj - 1) + 1;
            for(Ordinal ii=gridC()->si(fType, X, B::Y); ii<=this->iimax_[X]; ++ii) {
              Ordinal i = this->dd_[X]*(ii - 1) + 1;

              y(ii, jj, kk) = 0.;

              for(int kkk=0; kkk<=1; ++kkk)
                for(int jjj=-1; jjj<=1; ++jjj)
                  for(int iii=-1; iii<=1; ++iii)
                    y(ii, jj, kk) +=
                      cRS_[X](ii, iii)*cRS_[Y](jj, jjj)*cRV_[Z](kk, kkk)*x(i+iii, j+jjj, k+kkk) ;
            }
          }
        }
        break;
      }
      case(F::S) : {
        // todo: throw exception
        break;
      }
    }

    this->gather(y.getRawPtr());

    y.changed();
  }


  void print(std::ostream& out=std::cout) const {

    out << "=== Restriction OP ===\n";
    out << "nGather:\t" << this->nGather_ << "\n";
    out << "rankc2:\t" << this->rankc2_ << "\n";
    out << "comm2:\t" << this->comm2_ << "\n";

    out << " --- scalar stencil: ---";
    for(int j=0; j<3; ++j) {
      out << "\ndir: " << j << "\n";
      cRS_[j].print(out);
    }

    out << " --- velocity stencil: ---";
    for(int j=0; j<3; ++j) {
      out << "\ndir: " << j << "\n";
      cRV_[j].print(out);
    }
  }


  Teuchos::Tuple<Ordinal, dimension> getDD() const {
    return this->dd_;
  };

  Teuchos::RCP<const GridT> gridC() const {
    return this->gridC_;
  };
  Teuchos::RCP<const GridT> gridF() const {
    return this->gridF_;
  };

  const std::string getLabel() const {
    return "Restriction VF";
  };


}; // end of class RestrictionVFOp



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_RESTRICTIONVFOP_HPP
