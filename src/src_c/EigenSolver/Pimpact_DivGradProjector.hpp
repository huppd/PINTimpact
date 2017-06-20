#pragma once
#ifndef PIMPACT_DGPROJECTOR_HPP
#define PIMPACT_DGPROJECTOR_HPP


#include "Pimpact_DivOp.hpp"
#include "Pimpact_DivGradNullSpace.hpp"


namespace Pimpact {



template<class OperatorT>
class DGProjector {

  using RangeFieldT = typename OperatorT::RangeFieldT;
  using SpaceT = typename RangeFieldT::SpaceT;
  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  RangeFieldT nullspace_;

  void setCornersZero( RangeFieldT& rhs ) const {

    auto space = nullspace_.space();

    // BC XY
    if( space->getBCLocal()->getBCL(X)>0 && space->getBCLocal()->getBCL(Y)>0 ) {
      OT i = space->si(F::S,X,B::Y);
      OT j = space->si(F::S,Y,B::Y);
      for( OT k=space->si(F::S,Z,B::Y); k<=space->ei(F::S,Z,B::Y); ++k )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCL(X)>0 && space->getBCLocal()->getBCU(Y)>0 ) {
      OT i = space->si(F::S,X,B::Y);
      OT j = space->ei(F::S,Y,B::Y);
      for( OT k=space->si(F::S,Z,B::Y); k<=space->ei(F::S,Z,B::Y); ++k )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(X)>0 && space->getBCLocal()->getBCL(Y)>0 ) {
      OT i = space->ei(F::S,X,B::Y);
      OT j = space->si(F::S,Y,B::Y);
      for( OT k=space->si(F::S,Z,B::Y); k<=space->ei(F::S,Z,B::Y); ++k )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(X)>0 && space->getBCLocal()->getBCU(Y)>0 ) {
      OT i = space->ei(F::S,X,B::Y);
      OT j = space->ei(F::S,Y,B::Y);
      for( OT k=space->si(F::S,Z,B::Y); k<=space->ei(F::S,Z,B::Y); ++k )
        rhs(i,j,k) = 0.;
    }

    // BC XZ
    if( space->getBCLocal()->getBCL(X)>0 && space->getBCLocal()->getBCL(Z)>0 ) {
      OT i = space->si(F::S,X,B::Y);
      OT k = space->si(F::S,Z,B::Y);
      for( OT j=space->si(F::S,Y,B::Y); j<=space->ei(F::S,Y,B::Y); ++j )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCL(X)>0 && space->getBCLocal()->getBCU(Z)>0 ) {
      OT i = space->si(F::S,X,B::Y);
      OT k = space->ei(F::S,Z,B::Y);
      for( OT j=space->si(F::S,Y,B::Y); j<=space->ei(F::S,Y,B::Y); ++j )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(X)>0 && space->getBCLocal()->getBCL(Z)>0 ) {
      OT i = space->ei(F::S,X,B::Y);
      OT k = space->si(F::S,Z,B::Y);
      for( OT j=space->si(F::S,Y,B::Y); j<=space->ei(F::S,Y,B::Y); ++j )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(X)>0 && space->getBCLocal()->getBCU(Z)>0 ) {
      OT i = space->ei(F::S,X,B::Y);
      OT k = space->ei(F::S,Z,B::Y);
      for( OT j=space->si(F::S,Y,B::Y); j<=space->ei(F::S,Y,B::Y); ++j )
        rhs(i,j,k) = 0.;
    }

    // BC YZ
    if( space->getBCLocal()->getBCL(Y)>0 && space->getBCLocal()->getBCL(Z)>0 ) {
      OT j = space->si(F::S,Y,B::Y);
      OT k = space->si(F::S,Z,B::Y);
      for( OT i=space->si(F::S,X,B::Y); i<=space->ei(F::S,X,B::Y); ++i )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCL(Y)>0 && space->getBCLocal()->getBCU(Z)>0 ) {
      OT j = space->si(F::S,Y,B::Y);
      OT k = space->ei(F::S,Z,B::Y);
      for( OT i=space->si(F::S,X,B::Y); i<=space->ei(F::S,X,B::Y); ++i )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(Y)>0 && space->getBCLocal()->getBCL(Z)>0 ) {
      OT j = space->ei(F::S,Y,B::Y);
      OT k = space->si(F::S,Z,B::Y);
      for( OT i=space->si(F::S,X,B::Y); i<=space->ei(F::S,X,B::Y); ++i )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(Y)>0 && space->getBCLocal()->getBCU(Z)>0 ) {
      OT j = space->ei(F::S,Y,B::Y);
      OT k = space->ei(F::S,Z,B::Y);
      for( OT i=space->si(F::S,X,B::Y); i<=space->ei(F::S,X,B::Y); ++i )
        rhs(i,j,k) = 0.;
    }
    rhs.changed();
  }

public:

  DGProjector() {}
  DGProjector( const Teuchos::RCP<const OperatorT>& op):
    nullspace_( op->space() ) {

    DivGradNullSpace<DivOp<SpaceT> > compNullspace;

    compNullspace.computeNullSpace( op->getDivOp(), nullspace_ );
  }

  void operator()( RangeFieldT& rhs ) const {

    auto space = nullspace_.space();
    
    setCornersZero( rhs );
    ST bla = -nullspace_.dot( rhs );

    //if( 0==space->rankST() )
      //std::cout << "DivGrad^-1"<< ": nullspace contributtion: " << std::abs(bla)  << "\n";

    if( std::abs( bla )>0. )
      rhs.add( 1., rhs, bla, nullspace_ );

    setCornersZero( rhs );
  }

}; // end of class EmptyProjector

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DGPROJECTOR_HPP
