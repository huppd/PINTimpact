//#pragma once
#ifndef PIMPACT_PICARDPROJECTOR_HPP
#define PIMPACT_PICARDPROJECTOR_HPP


#include "Pimpact_CompoundField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_DivGradNullSpace.hpp"


namespace Pimpact {



template<class OperatorT>
class PicardProjector {

  using SpaceT = typename OperatorT::SpaceT;

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  using RangeFieldT = typename OperatorT::RangeFieldT;

  CompoundField<VectorField<SpaceT>, ScalarField<SpaceT>> nullspace_;

  void setCornersZero( ScalarField<SpaceT>& rhs ) const {

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

  PicardProjector() {}

  /// \todo compile test
  PicardProjector( const Teuchos::RCP<const OperatorT>& op):
    nullspace_( op->space() ) {

      auto space = nullspace_.space();

      DivGradNullSpace<DivOp<SpaceT> > compNullspace;

      compNullspace.computeNullSpace( op->getOpV2S()->getOperatorPtr(),
          nullspace_.getSField(), true );

      nullspace_.getVField()(Pimpact::F::U).initFromFunction(
          [&space]( ST x, ST y, ST z ) -> ST { return( ( (Pimpact::BC::Dirichlet==space->bcl(Pimpact::X)&&x<=0.)?1.:0.) + ( (Pimpact::BC::Dirichlet==space->bcu(Pimpact::X)&&1.<=x)?-1.:0.) ); } );
      nullspace_.getVField()(Pimpact::F::V).initFromFunction(
          [&space]( ST x, ST y, ST z ) -> ST { return( ( (Pimpact::BC::Dirichlet==space->bcl(Pimpact::Y)&&y<=0.)?1.:0.) + ( (Pimpact::BC::Dirichlet==space->bcu(Pimpact::Y)&&1.<=y)?-1.:0.) ); } );
      nullspace_.getVField()(Pimpact::F::W).initFromFunction(
          [&space]( ST x, ST y, ST z ) -> ST { return( ( (Pimpact::BC::Dirichlet==space->bcl(Pimpact::Z)&&z<=0.)?1.:0.) + ( (Pimpact::BC::Dirichlet==space->bcu(Pimpact::Z)&&1.<=z)?-1.:0.) ); } );

      ST blup = std::sqrt( 1./nullspace_.dot( nullspace_ ) );
      nullspace_.scale( blup );
    }


  /// \todo implment
  void operator()( RangeFieldT& rhs ) const {

    auto space = nullspace_.space();

    if( 0==space->si(F::U,3) ) {
      setCornersZero( rhs.getSField().get0Field() );

      ST bla = -nullspace_.getVField().dot( rhs.getVField().get0Field() );
      bla -= nullspace_.getSField().dot( rhs.getSField().get0Field() );

      if( 0==space->rankST() )
        std::cout << "Picard^-1"<< ": nullspace contributtion: " << std::abs(bla)  << "\n";

      if( std::abs( bla ) >= Teuchos::ScalarTraits<ST>::eps() ) {
        rhs.getVField().get0Field().add( 1., rhs.getVField().get0Field(), bla, nullspace_.getVField() );
        rhs.getSField().get0Field().add( 1., rhs.getSField().get0Field(), bla, nullspace_.getSField() );
      }

      setCornersZero( rhs.getSField().get0Field() );
    }
    for( OT i=std::max(space->si(F::U,3),1); i<=space->ei(F::U,3); ++i ) {
      {
        setCornersZero( rhs.getSField().getCField(i) );

        ST bla = -nullspace_.getVField().dot( rhs.getVField().getCField(i) );
        bla -= nullspace_.getSField().dot( rhs.getSField().getCField(i) );

        if( 0==space->rankST() )
          std::cout << "Picard^-1"<< ": nullspace contributtion: " << std::abs(bla)  << "\n";

        if( std::abs( bla ) >= Teuchos::ScalarTraits<ST>::eps() ) {
          rhs.getVField().getCField(i).add( 1., rhs.getVField().getCField(i), bla, nullspace_.getVField() );
          rhs.getSField().getCField(i).add( 1., rhs.getSField().getCField(i), bla, nullspace_.getSField() );
        }
        setCornersZero( rhs.getSField().getCField(i) );
      }
      {
        setCornersZero( rhs.getSField().getSField(i) );

        ST bla = -nullspace_.getVField().dot( rhs.getVField().getSField(i) );
        bla -= nullspace_.getSField().dot( rhs.getSField().getSField(i) );

        if( 0==space->rankST() )
          std::cout << "Picard^-1"<< ": nullspace contributtion: " << std::abs(bla)  << "\n";

        if( std::abs( bla ) >= Teuchos::ScalarTraits<ST>::eps() ) {
          rhs.getVField().getSField(i).add( 1., rhs.getVField().getSField(i), bla, nullspace_.getVField() );
          rhs.getSField().getSField(i).add( 1., rhs.getSField().getSField(i), bla, nullspace_.getSField() );
        }
        setCornersZero( rhs.getSField().getSField(i) );
      }
    }
  }

}; // end of class EmptyProjector

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_PICARDPROJECTOR_HPP
