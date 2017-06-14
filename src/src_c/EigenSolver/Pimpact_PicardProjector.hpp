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


  void project( VectorField<SpaceT>& rhs_v, ScalarField<SpaceT>& rhs_s ) const {
    setCornersZero( rhs_s );

    auto space = nullspace_.space();

    ST bla = -( nullspace_.getVField().dot( rhs_v ) +
        nullspace_.getSField().dot( rhs_s ) );

    if( 0==space->rankST() )
      std::cout << "Picard^-1"<< ": nullspace contributtion: " << std::abs(bla)  << "\n";

    if( std::abs( bla ) >= Teuchos::ScalarTraits<ST>::eps() ) {
      rhs_v.add( 1., rhs_v, bla, nullspace_.getVField() );
      rhs_s.add( 1., rhs_s, bla, nullspace_.getSField() );
    }

    setCornersZero( rhs_s );
  }

public:

  PicardProjector() {}

  /// \todo think about normal 
  PicardProjector( const Teuchos::RCP<const OperatorT>& op):
    nullspace_( op->space() ) {

      auto space = nullspace_.space();

      DivGradNullSpace<DivOp<SpaceT> > compNullspace;

      compNullspace.computeNullSpace( op->getOpV2S()->getOperatorPtr(),
          nullspace_.getSField(), true );

      nullspace_.getVField()(Pimpact::F::U).initFromFunction(
          [&space]( ST x, ST y, ST z ) -> ST { return( ( (Pimpact::BC::Dirichlet==space->bcl(Pimpact::X)&&x<=0.)?-1.:0.) + ( (Pimpact::BC::Dirichlet==space->bcu(Pimpact::X)&&1.<=x)?1.:0.) ); } );
      nullspace_.getVField()(Pimpact::F::V).initFromFunction(                                                     
          [&space]( ST x, ST y, ST z ) -> ST { return( ( (Pimpact::BC::Dirichlet==space->bcl(Pimpact::Y)&&y<=0.)?-1.:0.) + ( (Pimpact::BC::Dirichlet==space->bcu(Pimpact::Y)&&1.<=y)?1.:0.) ); } );
      nullspace_.getVField()(Pimpact::F::W).initFromFunction(                                                     
          [&space]( ST x, ST y, ST z ) -> ST { return( ( (Pimpact::BC::Dirichlet==space->bcl(Pimpact::Z)&&z<=0.)?-1.:0.) + ( (Pimpact::BC::Dirichlet==space->bcu(Pimpact::Z)&&1.<=z)?1.:0.) ); } );

      ST blup =  1./nullspace_.norm();
      nullspace_.write( 777 );
      nullspace_.scale( blup );
    }


  /// \todo implment
  void operator()( RangeFieldT& rhs ) const {

    auto space = nullspace_.space();

    if( 0==space->si(F::U,3) )
      project( rhs.getVField().get0Field(), rhs.getSField().get0Field() );

    for( OT i=std::max(space->si(F::U,3),1); i<=space->ei(F::U,3); ++i ) {
      project( rhs.getVField().getCField(i), rhs.getSField().getCField(i) );
      project( rhs.getVField().getSField(i), rhs.getSField().getSField(i) );
    }
  }

}; // end of class EmptyProjector

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_PICARDPROJECTOR_HPP
