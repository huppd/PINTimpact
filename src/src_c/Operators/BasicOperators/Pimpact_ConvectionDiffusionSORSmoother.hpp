#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP
#define PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP

#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"



namespace Pimpact {


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
template<class OperatorT>
class ConvectionDiffusionSORSmoother {

public:

  using SpaceT = typename OperatorT::SpaceT;

  using FluxFieldT = ScalarField<SpaceT>[3];

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  using SW = typename SpaceT::SW;

  ST omega_;
  int nIter_;

  int ordering_;

  Teuchos::Tuple<short int,3> dirs_;

  const Teuchos::RCP<const OperatorT> op_;

  constexpr ST getHC( const ECoord dir, const F ftype, const OT i, const OT ii ) {
    return op_->getHelmOp()->getC(dir,ftype,i,ii);
  }

public:

   /// \brief Basic constructor for ConvectionDiffusionSORSmoother.
   /// 
   /// This constructor accepts the operatore to be smoothed in addition to a parameter
   /// list of options for the solver manager. These options include the following:
   ///  - "omega" - a \c double smoothing factor. Default: 1. 
   ///  - "numIters" - a \c int specifying the maximum number of iterations the underlying solver is allowed to perform. Default: 1
   ///  - "ordering" - a \c int ordinering (0: specified by dirs, 1: all direction, 2:
   ///    SHL). Default: 1. 
   ///  - "dir X" - a \c short int if in X direction. Default: 1. 
   ///  - "dir Y" - a \c short int if in Y direction. Default: 1. 
   ///  - "dir Z" - a \c short int if in Z direction. Default: 1. 
  ConvectionDiffusionSORSmoother(
    const Teuchos::RCP<const OperatorT>& op,
    Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
    omega_( pl->get("omega", 1. ) ),
    nIter_( pl->get("numIters", 1 ) ),
    ordering_( pl->get("Ordering",1 ) ),
    op_(op) {

    if( 4==SpaceT::dimNC )
      if( 0==op->space()->rankST() )
        std::cout << "Warning!!! ConvectionDiffusionSORSmoother strange behavior for dimNC=4, problems at outflow\n";

    if( 0==ordering_ ) {
      dirs_[0] = pl->get<short int>( "dir X", 1 );
      dirs_[1] = pl->get<short int>( "dir Y", 1 );
      dirs_[2] = pl->get<short int>( "dir Z", 1 );
    }
  }



  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, const ST mulI, const ST mulC, const ST mulL, const Add add=Add::N ) const {
    std::cout << "not implmented\n";
  }

  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, const Add add=Add::N ) const {

    // testing field consistency
    assert( z.getType() == y.getType() );

    for( int i=0; i<SpaceT::sdim; ++i )
      assert( x[i].getType() == y.getType() );

    // exchange wind and "rhs"
    for( int vel_dir=0; vel_dir<SpaceT::sdim; ++vel_dir )
      x[vel_dir].exchange();

    for( int i=0; i<nIter_; ++i ) {
      switch( ordering_ ) {
        case 0 : {
          apply( x, y, z, dirs_ );
          break;
        }
        case 2 :{
          applySHL( x, y, z );
          break;
        }
        default: {
          applyNPoint( x, y, z );
          break;
        }
      }
    }
  }

protected:

  void applySHL( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z ) const {

    Teuchos::Tuple<int,SpaceT::dimension> ib = space()->getProcGrid()->getNP();

    Teuchos::Tuple<int,3>	dirS;
    Teuchos::Tuple<int,3>	inc;

    for( int i=0; i<3; ++i ) {
      if( (ib[i]-1)%2 == 0 ) {
        dirS[i] = -1;
        inc[i] = 2;
      } else {
        dirS[i] = 1;
        inc[i] = -2;
      }
    }

    dirs_[1] = 1;
    if( 3==SpaceT::sdim )
      for( dirs_[2]=dirS[2]; std::abs(dirs_[2])<=1; dirs_[2]+=inc[2] )
          for( dirs_[0]=dirS[0]; std::abs(dirs_[0])<=1; dirs_[0]+=inc[0] )
            apply( x, y, z, dirs_ );
    else {
      dirs_[2] = 1 ;
      for( dirs_[1]=-1; dirs_[1]<2; dirs_[1]+=2 )
        for( dirs_[0]=-1; dirs_[0]<2; dirs_[0]+=2 )
          apply( x, y, z, dirs_ );
    }
  }

  void applyNPoint( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z ) const {

    Teuchos::Tuple<int,SpaceT::dimension> ib = space()->getProcGrid()->getNP();

    Teuchos::Tuple<int,3>	dirS;
    Teuchos::Tuple<int,3>	inc;

    for( int i=0; i<3; ++i ) {
      if( (ib[i]-1)%2 == 0 ) {
        dirS[i] = -1;
        inc[i] = 2;
      } else {
        dirS[i] = 1;
        inc[i] = -2;
      }
    }

    dirs_[0] = -1;
    if( 3==SpaceT::sdim )
      for( dirs_[2]=dirS[2]; std::abs(dirs_[2])<=1; dirs_[2]+=inc[2] )
        for( dirs_[1]=dirS[1]; std::abs(dirs_[1])<=1; dirs_[1]+=inc[1] )
          for( dirs_[0]=dirS[0]; std::abs(dirs_[0])<=1; dirs_[0]+=inc[0] )
            apply( x, y, z, dirs_ );
    else {
      dirs_[2] = 1 ;
      for( dirs_[1]=-1; dirs_[1]<2; dirs_[1]+=2 )
        for( dirs_[0]=-1; dirs_[0]<2; dirs_[0]+=2 )
          apply( x, y, z, dirs_ );
    }
  }


  /// \brief little helper
  void apply( const FluxFieldT& wind, const DomainFieldT& b, RangeFieldT& x,
              const Teuchos::Tuple<short int,3>& dirs ) const {

    for( int i=0; i<3; ++i ) {
      assert( 1==dirs[i] || -1==dirs[i] );
    }

    const F f = x.getType();
    x.exchange();

    applyBC( b, x );

    Teuchos::Tuple<OT,3> ss;
    Teuchos::Tuple<OT,3> nn;
    for( int i=0; i<3; ++i ) {
      ss[i] = (dirs[i]>0)?(space()->si(f,i,B::N)  ):(space()->ei(f,i,B::N)  );
      nn[i] = (dirs[i]>0)?(space()->ei(f,i,B::N)+1):(space()->si(f,i,B::N)-1);
    }

    if( 3==SpaceT::sdim ) {

      for( OT k=ss[Z]; k!=nn[Z]; k+=dirs[Z] )
        for( OT j=ss[Y]; j!=nn[Y]; j+=dirs[Y] )
          for( OT i=ss[X]; i!=nn[X]; i+=dirs[X] ) {
            ST diag =
              op_->getMulI()
              + op_->getMulC() * op_->getConvSOp()->innerDiag3D(
                wind[0](i,j,k),
                wind[1](i,j,k),
                wind[2](i,j,k), f, i, j, k )
              - op_->getMulL() * op_->getHelmOp()->innerDiag3D( f, i, j, k) ;
            assert( diag!=0 );
            x(i,j,k) += omega_*( b(i,j,k)
                                 - op_->getMulI() * x(i,j,k)
                                 - op_->getMulC() * op_->getConvSOp()->innerStenc3D(
                                   wind[0](i,j,k),
                                   wind[1](i,j,k),
                                   wind[2](i,j,k), x, i, j, k )
                                 + op_->getMulL() * op_->getHelmOp()->innerStenc3D( x, f, i, j, k) ) / diag;
          }
    } else {

      for( OT k=ss[Z]; k!=nn[Z]; k+=dirs[Z] )
        for( OT j=ss[Y]; j!=nn[Y]; j+=dirs[Y] )
          for( OT i=ss[X]; i!=nn[X]; i+=dirs[X] ) {
            ST diag =
              op_->getMulI()
              + op_->getMulC() * op_->getConvSOp()->innerDiag2D(
                wind[0](i,j,k),
                wind[1](i,j,k), f, i, j, k )
              - op_->getMulL() * op_->getHelmOp()->innerDiag2D( f, i, j, k) ;
            assert( diag!=0 );
            x(i,j,k) += omega_*( b(i,j,k)
                                 - op_->getMulI() * x(i,j,k)
                                 - op_->getMulC() * op_->getConvSOp()->innerStenc2D(
                                   wind[0](i,j,k),
                                   wind[1](i,j,k), x, i, j, k )
                                 + op_->getMulL() * op_->getHelmOp()->innerStenc2D( x, f, i, j, k) ) / diag;
          }
    }
    x.changed();
  }


  /// \brief implements smoothing for Dirichlet boundary conditions as identity
  /// in tangential / velocity direction or interpolation in wand normal
  /// direction
  void applyBC( const DomainFieldT& b, RangeFieldT& y	) const {

    assert( b.getType()==y.getType() );

    const F f = y.getType();

    const ST omegaBC = omega_;
    //const ST omegaBC = 0.9;

    // U-field
    if( F::U==f ) {

      // tangential direction: Y
      if( 0<space()->bcl(Y) ) {
        OT j = space()->si(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST temp = 0.;
            for( OT jj=0; jj<=SW::BU(Y); ++jj )
              temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }
      if( 0<space()->bcu(Y) ) {
        OT j = space()->ei(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST temp = 0.;
            for( OT jj=SW::BL(Y); jj<=0; ++jj )
              temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }

      // tangential direction: Z
      if( 0<space()->bcl(Z) ) {
        OT k = space()->si(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST temp = 0.;
            for( OT kk=0; kk<=SW::BU(Z); ++kk )
              temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }
      if( 0<space()->bcu(Z) ) {
        OT k = space()->ei(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST temp = 0.;
            for( OT kk=SW::BL(Z); kk<=0; ++kk )
              temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }

      // normal direction: X
      if( 0<space()->bcl(X) ) {
        OT i = space()->si(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
            ST temp = 0.;
            for( OT ii=0; ii<=SW::BU(X); ++ii )
              temp += getHC(X,f,i,ii)*y(i+ii,j,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }
      if( 0<space()->bcu(X) ) {
        OT i = space()->ei(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
            ST temp = 0.;
            for( OT ii=SW::BL(X); ii<=0; ++ii )
              temp += getHC(X,f,i,ii)*y(i+ii,j,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }
    }

    // V-field
    if( F::V==f ) {

      // tangential direction: X
      if( 0<space()->bcl(X) ) {
        OT i = space()->si(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j ) {
            ST temp = 0.;
            for( OT ii=0; ii<=SW::BU(X); ++ii )
              temp += getHC(X,f,i,ii)*y(i+ii,j,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }
      if( 0<space()->bcu(X) ) {
        OT i = space()->ei(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j ) {
            ST temp = 0.;
            for( OT ii=SW::BL(X); ii<=0; ++ii )
              temp += getHC(X,f,i,ii)*y(i+ii,j,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }

      // tangential direction: Z
      if( 0<space()->bcl(Z) ) {
        OT k = space()->si(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT kk=0; kk<=SW::BU(Z); ++kk )
              temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }
      if( 0<space()->bcu(Z) ) {
        OT k = space()->ei(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT kk=SW::BL(Z); kk<=0; ++kk )
              temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }

      // normal direction: Y
      if( 0<space()->bcl(Y) ) {
        OT j = space()->si(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT jj=0; jj<=SW::BU(Y); ++jj )
              temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }
      if( 0<space()->bcu(Y) ) {
        OT j = space()->ei(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT jj=SW::DL(Y); jj<=SW::DU(Y); ++jj )
              temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }
    }

    // W-field
    if( F::W==f ) {

      // tangential direction: X
      if( 0<space()->bcl(X) ) {
        OT i = space()->si(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
          for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
            ST temp = 0.;
            for( OT ii=0; ii<=SW::BU(X); ++ii )
              temp += getHC(X,f,i,ii)*y(i+ii,j,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }
      if( 0<space()->bcu(X) ) {
        OT i = space()->ei(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
          for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
            ST temp = 0.;
            for( OT ii=SW::BL(X); ii<=0; ++ii )
              temp += getHC(X,f,i,ii)*y(i+ii,j,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }

      // tangential direction: Y
      if( 0<space()->bcl(Y) ) {
        OT j = space()->si(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT jj=0; jj<=SW::BU(Y); ++jj )
              temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }
      if( 0<space()->bcu(Y) ) {
        OT j = space()->ei(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT jj=SW::BL(Y); jj<=0; ++jj )
              temp += getHC(Y,f,j,jj)*y(i,j+jj,k);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }

      // normal direction: Z
      if( 0<space()->bcl(Z) ) {
        OT k = space()->si(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT kk=0; kk<=SW::BU(Z); ++kk )
              temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }
      if( 0<space()->bcu(Z) ) {
        OT k = space()->ei(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT kk=SW::BL(Z); kk<=0; ++kk )
              temp += getHC(Z,f,k,kk)*y(i,j,k+kk);
            y(i,j,k) += omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }
    }
  }

public:

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << "---\n";
    op_->print();
  }


  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op_->space();
  }

  void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}


  constexpr const std::string getLabel() const {
    return "ConvectionDiffusionSORSmoother";
  };


}; // end of class ConvectionDiffusionSORSmoother


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP
