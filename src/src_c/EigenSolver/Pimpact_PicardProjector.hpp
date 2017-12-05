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
  CompoundField<VectorField<SpaceT>, ScalarField<SpaceT>> projection_;

  ST dotNP_;

  void setCornersZero( ScalarField<SpaceT>& rhs ) const {

    auto space = nullspace_.space();
    F f = rhs.getType();

    // BC XY
    if( space->getBCLocal()->getBCL(X)>0 && space->getBCLocal()->getBCL(Y)>0 ) {
      OT i = space->si(f,X,B::Y);
      OT j = space->si(f,Y,B::Y);
      for( OT k=space->si(f,Z,B::Y); k<=space->ei(f,Z,B::Y); ++k )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCL(X)>0 && space->getBCLocal()->getBCU(Y)>0 ) {
      OT i = space->si(f,X,B::Y);
      OT j = space->ei(f,Y,B::Y);
      for( OT k=space->si(f,Z,B::Y); k<=space->ei(f,Z,B::Y); ++k )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(X)>0 && space->getBCLocal()->getBCL(Y)>0 ) {
      OT i = space->ei(f,X,B::Y);
      OT j = space->si(f,Y,B::Y);
      for( OT k=space->si(f,Z,B::Y); k<=space->ei(f,Z,B::Y); ++k )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(X)>0 && space->getBCLocal()->getBCU(Y)>0 ) {
      OT i = space->ei(f,X,B::Y);
      OT j = space->ei(f,Y,B::Y);
      for( OT k=space->si(f,Z,B::Y); k<=space->ei(f,Z,B::Y); ++k )
        rhs(i,j,k) = 0.;
    }

    // BC XZ
    if( space->getBCLocal()->getBCL(X)>0 && space->getBCLocal()->getBCL(Z)>0 ) {
      OT i = space->si(f,X,B::Y);
      OT k = space->si(f,Z,B::Y);
      for( OT j=space->si(f,Y,B::Y); j<=space->ei(f,Y,B::Y); ++j )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCL(X)>0 && space->getBCLocal()->getBCU(Z)>0 ) {
      OT i = space->si(f,X,B::Y);
      OT k = space->ei(f,Z,B::Y);
      for( OT j=space->si(f,Y,B::Y); j<=space->ei(f,Y,B::Y); ++j )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(X)>0 && space->getBCLocal()->getBCL(Z)>0 ) {
      OT i = space->ei(f,X,B::Y);
      OT k = space->si(f,Z,B::Y);
      for( OT j=space->si(f,Y,B::Y); j<=space->ei(f,Y,B::Y); ++j )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(X)>0 && space->getBCLocal()->getBCU(Z)>0 ) {
      OT i = space->ei(f,X,B::Y);
      OT k = space->ei(f,Z,B::Y);
      for( OT j=space->si(f,Y,B::Y); j<=space->ei(f,Y,B::Y); ++j )
        rhs(i,j,k) = 0.;
    }

    // BC YZ
    if( space->getBCLocal()->getBCL(Y)>0 && space->getBCLocal()->getBCL(Z)>0 ) {
      OT j = space->si(f,Y,B::Y);
      OT k = space->si(f,Z,B::Y);
      for( OT i=space->si(f,X,B::Y); i<=space->ei(f,X,B::Y); ++i )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCL(Y)>0 && space->getBCLocal()->getBCU(Z)>0 ) {
      OT j = space->si(f,Y,B::Y);
      OT k = space->ei(f,Z,B::Y);
      for( OT i=space->si(f,X,B::Y); i<=space->ei(f,X,B::Y); ++i )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(Y)>0 && space->getBCLocal()->getBCL(Z)>0 ) {
      OT j = space->ei(f,Y,B::Y);
      OT k = space->si(f,Z,B::Y);
      for( OT i=space->si(f,X,B::Y); i<=space->ei(f,X,B::Y); ++i )
        rhs(i,j,k) = 0.;
    }
    if( space->getBCLocal()->getBCU(Y)>0 && space->getBCLocal()->getBCU(Z)>0 ) {
      OT j = space->ei(f,Y,B::Y);
      OT k = space->ei(f,Z,B::Y);
      for( OT i=space->si(f,X,B::Y); i<=space->ei(f,X,B::Y); ++i )
        rhs(i,j,k) = 0.;
    }
    rhs.changed();
  }


  void project( VectorField<SpaceT>& rhs_v, ScalarField<SpaceT>& rhs_s ) const {

    //setCornersZero( rhs_s );

    auto space = nullspace_.space();

    ST bla = -( nullspace_.getVField().dot( rhs_v ) + nullspace_.getSField().dot( rhs_s )
        )/dotNP_;

    if( 0==space->rankST() )
      std::cout << "Picard^-1"<< ": nullspace contributtion: " << std::abs(bla)  << "\n";

    if( std::abs( bla ) >= Teuchos::ScalarTraits<ST>::eps() ) {
      rhs_v.add( 1., rhs_v, bla, projection_.getVField() );
      rhs_s.add( 1., rhs_s, bla, projection_.getSField() );
    }

    setCornersZero( rhs_s );
  }

public:

  PicardProjector() {}

  PicardProjector( const Teuchos::RCP<const OperatorT>& op):
    nullspace_( op->space() ), projection_( op->space() ) {

    auto space = nullspace_.space();

    DivGradNullSpace<DivOp<SpaceT> > compNullspace;

    compNullspace.computeNullSpace( op->getOpV2S()->getOperatorPtr(),
        nullspace_.getSField(), true );

    // U inflow
    if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::X) ) {
      OT si = space()->si(F::U,X,B::Y);
      for( OT k=space()->si(F::U,Z,B::Y); k<=space()->ei(F::U,Z,B::Y); ++k )
        for( OT j=space()->si(F::U,Y,B::Y); j<=space()->ei(F::U,Y,B::Y); ++j )
          nullspace_.getVField()(F::U)(si,j,k) = -1.;
    }
    // U outflow
    if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::X) ) {
      OT ei = space()->ei(F::U,X,B::Y);
      for( OT k=space()->si(F::U,Z,B::Y); k<=space()->ei(F::U,Z,B::Y); ++k )
        for( OT j=space()->si(F::U,Y,B::Y); j<=space()->ei(F::U,Y,B::Y); ++j )
          nullspace_.getVField()(F::U)(ei,j,k) =  1.;
    }

    // V inflow
    if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::Y) ) {
      OT sj = space()->si(F::V,Y,B::Y);
      for( OT k=space()->si(F::V,Z,B::Y); k<=space()->ei(F::V,Z,B::Y); ++k )
        for( OT i=space()->si(F::V,X,B::Y); i<=space()->ei(F::V,X,B::Y); ++i )
          nullspace_.getVField()(F::V)(i,sj,k) = -1.;
    }
    // V outflow
    if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::Y) ) {
      OT ej = space()->ei(F::V,Y,B::Y);
      for( OT k=space()->si(F::V,Z,B::Y); k<=space()->ei(F::V,Z,B::Y); ++k )
        for( OT i=space()->si(F::V,X,B::Y); i<=space()->ei(F::V,X,B::Y); ++i )
          nullspace_.getVField()(F::V)(i,ej,k) = 1.;
    }

    // W in/outflow
    if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::Z) ) {
      OT sk = space()->si(F::W,Z,B::Y);
      for( OT j=space()->si(F::W,Y,B::Y); j<=space()->ei(F::W,Y,B::Y); ++j )
        for( OT i=space()->si(F::W,X,B::Y); i<=space()->ei(F::W,X,B::Y); ++i )
          nullspace_.getVField()(F::W)(i,j,sk) = -1.;
    }
    // W in/outflow
    if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::Z) ) {
      OT ek = space()->ei(F::W,Z,B::Y);
      for( OT j=space()->si(F::W,Y,B::Y); j<=space()->ei(F::W,Y,B::Y); ++j )
        for( OT i=space()->si(F::W,X,B::Y); i<=space()->ei(F::W,X,B::Y); ++i )
          nullspace_.getVField()(F::W)(i,j,ek) = 1.;
    }
    nullspace_.getVField().changed();

    setCornersZero( projection_.getSField() );

    //ST blup =  1./nullspace_.norm();
    //nullspace_.scale( blup );
    //nullspace_.write( 777 );

    // outflow projection
    const ST pi = 4.*std::atan(1.);
    const ST width = 0.95;
    const ST eps = 0.;

    auto scalefunc = [=]( ST x, ST y, ST z ) ->ST {
      return (y<=width)?0.:( std::cos( pi*(y-width)/(1.-width) + pi )/2. + 0.5); };

    if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::Y) ) {
      // outflow velocity
      OT ej = space()->ei(F::V,Y,B::Y);
      for( OT k=space()->si(F::V,Z,B::Y); k<=space()->ei(F::V,Z,B::Y); ++k )
        for( OT i=space()->si(F::V,X,B::Y); i<=space()->ei(F::V,X,B::Y); ++i )
          projection_.getVField()(F::V)(i,ej,k) = 1.;

      //// pressure "outflow" 
      //OT ej = space()->ei(F::S,Y,B::Y);
      //for( OT k=space()->si(F::S,Z,B::Y); k<=space()->ei(F::S,Z,B::Y); ++k )
        //for( OT i=space()->si(F::S,X,B::Y); i<=space()->ei(F::S,X,B::Y); ++i )
          //projection_.getSField()(i,ej,k) = 1.;
    }

    projection_.getSField().initFromFunction( scalefunc );
    //setCornersZero( projection_.getVField() );
    setCornersZero( projection_.getSField() );

    dotNP_ = nullspace_.dot( projection_ );

    if( space->rankST()==0 ) std::cout << "dotNP: " << dotNP_ << "\n";
    assert( dotNP_!=0. );
  }


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
