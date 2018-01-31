#pragma once
#ifndef PIMPACT_DIVGRADNULLSPACE_HPP
#define PIMPACT_DIVGRADNULLSPACE_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialQRDenseSolver.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include "Pimpact_Stencil.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Space.hpp" // just for createOstream<>




namespace Pimpact {



/// \todo make nullspace greate again, (store vectors) not hole field
/// \todo integrate with DGProjector
template<class OperatorT>
class DivGradNullSpace {

protected:

  using SpaceT = typename OperatorT::SpaceT;

  using Scalar  = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  static const int dimNC = SpaceT::dimNC;
  static const int dim = SpaceT::dimension;

  using SW = StencilWidths<dim, dimNC>;

  using Stenc = Stencil<Scalar, Ordinal, SW::BL(0), SW::BL(0), SW::BU(0) >;

  using DomainFieldT= typename OperatorT::DomainFieldT;
  using RangeFieldT = typename OperatorT::RangeFieldT;

  using VectorT = Teuchos::SerialDenseVector<Ordinal, Scalar>;
  using MatrixT = Teuchos::SerialDenseMatrix<Ordinal, Scalar>;

  using SolverT = Teuchos::SerialQRDenseSolver<Ordinal, Scalar>;

public:

  /// \todo think about changing BC solver to proper one (inner field, setting last component one and move it to the rhs)
  void computeNullSpace(const Teuchos::RCP<const OperatorT>& div, RangeFieldT& y, const bool DJG_yes = true)  {

    Teuchos::RCP<const SpaceT> space = div->space();

    Teuchos::Tuple<Teuchos::RCP<VectorT>, SpaceT::sdim > x_ ;

    for(int dir=0; dir<SpaceT::sdim; ++dir) {
      const Ordinal N = space->nGlo(dir);
      x_[dir] = Teuchos::rcp(new VectorT(N));

      if(-1==space->getBCGlobal()->getBCL(dir)) {
        *x_[dir] = 1.;
      } else {
        // global stencil
        Ordinal nTempG = (space->nGlo(dir) + space->bu(dir) - space->bl(dir) + 1)
                         *(space->bu(dir) - space->bl(dir) + 1);

        Stenc cG1(space->nGlo(dir) + space->bu(dir));
        Stenc cG2(space->nGlo(dir) + space->bu(dir));

        for(Ordinal i = space->si(F::S, dir); i<=space->ei(F::S, dir); ++i)
          for(Ordinal ii = space->dl(dir); ii<=space->du(dir); ++ii)
            cG1(i+space->getShift(dir), ii)= div->getC(static_cast<ECoord>(dir), i, ii);

        /// \todo Allgather would be better
        MPI_Allreduce(
          cG1.get(),    		                        // const void *sendbuf,
          cG2.get(),    		                        // void *recvbuf,
          nTempG,			                              // int count,
          MPI_REAL8,	                              // MPI_Datatype datatype,
          MPI_SUM,		                              // MPI_Op op,
          space->getProcGrid()->getCommBar(dir));	// MPI_Comm comm)


        // generate global Div Stencil(copy from transposed stencil)
        // generate Matrix
        Teuchos::RCP<MatrixT> Dtrans = Teuchos::rcp(new MatrixT(N+1, N, true));
        for(Ordinal i=0; i<=N; ++i) {
          for(int ii=space->gl(dir); ii<=space->gu(dir); ++ii) {
            if(0<=i+ii-1 && i+ii-1<N)
              (*Dtrans)(i, i+ii-1) = cG2(i+ii, -ii);
          }
        }

        // generate RHS
        Teuchos::RCP<VectorT> b = Teuchos::rcp(new VectorT(N + 1, true));

        if(space->getBCGlobal()->getBCL(dir)>0) {
          if(DJG_yes) {

            for(int ii=-1; ii<=space->du(dir); ++ii)
              (*b)(ii+1) = -space->getInterpolateV2S()->getC(static_cast<Pimpact::ECoord>(dir), 1, ii);

            MPI_Bcast(
              &(*b)(0),																	// void* data,
              space->du(dir)+2,													// int count,
              MPI_DOUBLE,																// MPI_Datatype datatype,
              0,																				// int root,
              space->getProcGrid()->getCommBar(dir));	// MPI_Comm comm)
          } else
            (*b)(0) = -1.;
        }

        if(space->getBCGlobal()->getBCU(dir)>0) {
          if(DJG_yes) {

            for(int ii=space->dl(dir); ii<=0; ++ii)
              (*b)(N+ii) = space->getInterpolateV2S()->getC(static_cast<Pimpact::ECoord>(dir), space->ei(F::S, dir), ii);

            MPI_Bcast(
              &(*b)(N-3),																//void* data,
              1-space->dl(dir),													//int count,
              MPI_DOUBLE,																//MPI_Datatype datatype,
              space->getProcGrid()->getNP(dir)-1,				//int root,
              space->getProcGrid()->getCommBar(dir));	// MPI_Comm comm)
          } else
            (*b)(N) = 1.;
        }

        SolverT solv;
        solv.factorWithEquilibration(true);
        solv.setMatrix(Dtrans);
        solv.factor();
        solv.setVectors(x_[dir], b);
        solv.solve();
        // solve with QRU
        if(0==space->rankST()) {
          //std::cout <<"dir: " <<dir <<"\n";
          //std::cout <<std::scientific;
          //std::cout <<std::setprecision(std::numeric_limits<long double>::digits10 + 1);
          //std::cout <<*x_[dir] <<"\n";

          //Teuchos::RCP<std::ostream> output = Pimpact::createOstream("null.txt");
          //x_[dir]->print(*output);
        }
      }
    }
    if(3==SpaceT::sdim) {
      for(Ordinal k=space()->si(F::S, Z); k<=space()->ei(F::S, Z); ++k)
        for(Ordinal j=space()->si(F::S, Y); j<=space()->ei(F::S, Y); ++j)
          for(Ordinal i=space()->si(F::S, X); i<=space()->ei(F::S, X); ++i) {
            y(i, j, k) = (*x_[0])(i-1+space->getShift(0))* (*x_[1])(j-1+space->getShift(1))* (*x_[2])(k-1+space->getShift(2));
          }
    } else {
      for(Ordinal k=space()->si(F::S, Z); k<=space()->ei(F::S, Z); ++k)
        for(Ordinal j=space()->si(F::S, Y); j<=space()->ei(F::S, Y); ++j)
          for(Ordinal i=space()->si(F::S, X); i<=space()->ei(F::S, X); ++i)
            y(i, j, k) = (*x_[0])(i-1+space->getShift(0))* (*x_[1])(j-1+space->getShift(1));
    }
    //Scalar blup = std::sqrt(1./y.dot(y));
    //y.scale(blup);
  }


}; // end of class SimpleVectorIteration


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DIVGRADNULLSPACE_HPP
