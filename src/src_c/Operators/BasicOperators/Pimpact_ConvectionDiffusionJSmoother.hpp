#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP
#define PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP

#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"



namespace Pimpact {


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
/// \note todo merge with SORSmoother or make interface
template<class OperatorT>
class ConvectionDiffusionJSmoother {

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

  const Teuchos::RCP<const OperatorT> op_;

  constexpr ST getHC( const ECoord dir, const F ftype, OT i, OT ii ) {
    return op_->getHelmOp()->getC(dir,ftype,i,ii);
  }

public:

  /// \brief constructor
  ///
  /// These options include the following:
  /// - "omega" - damping parameter
  /// - "numIters" - an \c int specifying the maximum number of iterations the
  ConvectionDiffusionJSmoother(
    const Teuchos::RCP<const OperatorT>& op,
    Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
    omega_( pl->get<ST>("omega", 0.5 ) ),
    nIter_( pl->get("numIters", 10 ) ),
    op_(op) {}



  void apply( const FluxFieldT& wind, const DomainFieldT& x, RangeFieldT& y, ST mul, ST mulI, ST mulC, ST mulL ) const {
    std::cout << "not implmented\n";
  }


protected:

  void applyStep( const FluxFieldT& wind, const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y ) const {

    const F f = y.getType();

    x.exchange();

    //applyBC( b, x, y );

    if( 3==SpaceT::sdim ) {
      for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
        for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST diag =
              op_->getMulI()
              + op_->getMulC() * op_->getConvSOp()->innerDiag3D(
                wind[0](i,j,k),
                wind[1](i,j,k),
                wind[2](i,j,k), f, i, j, k )
              - op_->getMulL() * op_->getHelmOp()->innerDiag3D( f, i, j, k) ;
            assert( diag!=0 );
            y(i,j,k) = x(i,j,k) + omega_*( b(i,j,k)
                                           - op_->getMulI() * x(i,j,k)
                                           - op_->getMulC() * op_->getConvSOp()->innerStenc3D(
                                             wind[0](i,j,k),
                                             wind[1](i,j,k),
                                             wind[2](i,j,k), x, i, j, k )
                                           + op_->getMulL() * op_->getHelmOp()->innerStenc3D( x, f, i, j, k) ) / diag;
          }
    } else {

      for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
        for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST diag =
              op_->getMulI()
              + op_->getMulC() * op_->getConvSOp()->innerDiag2D(
                wind[0](i,j,k),
                wind[1](i,j,k), f, i, j, k )
              - op_->getMulL() * op_->getHelmOp()->innerDiag2D( f, i, j, k) ;
            assert( diag!=0 );
            y(i,j,k) = x(i,j,k) + omega_*( b(i,j,k)
                                           - op_->getMulI() * x(i,j,k)
                                           - op_->getMulC() * op_->getConvSOp()->innerStenc2D(
                                             wind[0](i,j,k),
                                             wind[1](i,j,k), x, i, j, k )
                                           + op_->getMulL() * op_->getHelmOp()->innerStenc2D( x, f, i, j, k) ) / diag;
          }
    }

    applyBC( b, x, y );
    y.changed();
  }

  /// \brief implements smoothing for Dirichlet boundary conditions as identity
  /// in tangential/velocity direction or interpolation in wand normal
  /// direction
  void applyBC( const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y ) const {

    assert( b.getType()==y.getType() );
    assert( x.getType()==y.getType() );

    const F f = y.getType();

    const ST omegaBC = omega_;

    // U-field
    if( F::U==f ) {

      // tangential direction: Y
      if( 0<space()->bcl(Y) ) {
        OT j = space()->si(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST temp = 0.;
            for( OT jj=0; jj<=SW::BU(Y); ++jj )
              temp += getHC(Y,f,j,jj)*x(i,j+jj,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }
      if( 0<space()->bcu(Y) ) {
        OT j = space()->ei(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST temp = 0.;
            for( OT jj=SW::BL(Y); jj<=0; ++jj )
              temp += getHC(Y,f,j,jj)*x(i,j+jj,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }

      // tangential direction: Z
      if( 0<space()->bcl(Z) ) {
        OT k = space()->si(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST temp = 0.;
            for( OT kk=0; kk<=SW::BU(Z); ++kk )
              temp += getHC(Z,f,k,kk)*x(i,j,k+kk);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }
      if( 0<space()->bcu(Z) ) {
        OT k = space()->ei(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
          for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
            ST temp = 0.;
            for( OT kk=SW::BL(Z); kk<=0; ++kk )
              temp += getHC(Z,f,k,kk)*x(i,j,k+kk);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }

      // normal direction: X
      if( 0<space()->bcl(X) ) {
        OT i = space()->si(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
            ST temp = 0.;
            for( OT ii=0; ii<=SW::BU(X); ++ii )
              temp += getHC(X,f,i,ii)*x(i+ii,j,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }
      if( 0<space()->bcu(X) ) {
        OT i = space()->ei(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
            ST temp = 0.;
            for( OT ii=SW::BL(X); ii<=0; ++ii )
              temp += getHC(X,f,i,ii)*x(i+ii,j,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
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
              temp += getHC(X,f,i,ii)*x(i+ii,j,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }
      if( 0<space()->bcu(X) ) {
        OT i = space()->ei(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j ) {
            ST temp = 0.;
            for( OT ii=SW::BL(X); ii<=0; ++ii )
              temp += getHC(X,f,i,ii)*x(i+ii,j,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }

      // tangential direction: Z
      if( 0<space()->bcl(Z) ) {
        OT k = space()->si(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT kk=0; kk<=SW::BU(Z); ++kk )
              temp += getHC(Z,f,k,kk)*x(i,j,k+kk);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }
      if( 0<space()->bcu(Z) ) {
        OT k = space()->ei(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT kk=SW::BL(Z); kk<=0; ++kk )
              temp += getHC(Z,f,k,kk)*x(i,j,k+kk);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }

      // normal direction: Y
      if( 0<space()->bcl(Y) ) {
        OT j = space()->si(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT jj=0; jj<=SW::BU(Y); ++jj )
              temp += getHC(Y,f,j,jj)*x(i,j+jj,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }
      if( 0<space()->bcu(Y) ) {
        OT j = space()->ei(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT jj=SW::DL(Y); jj<=SW::DU(Y); ++jj )
              temp += getHC(Y,f,j,jj)*x(i,j+jj,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
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
              temp += getHC(X,f,i,ii)*x(i+ii,j,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }
      if( 0<space()->bcu(X) ) {
        OT i = space()->ei(f,X,B::Y);
        for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
          for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
            ST temp = 0.;
            for( OT ii=SW::BL(X); ii<=0; ++ii )
              temp += getHC(X,f,i,ii)*x(i+ii,j,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(X,f,i,0);
          }
      }

      // tangential direction: Y
      if( 0<space()->bcl(Y) ) {
        OT j = space()->si(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT jj=0; jj<=SW::BU(Y); ++jj )
              temp += getHC(Y,f,j,jj)*x(i,j+jj,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }
      if( 0<space()->bcu(Y) ) {
        OT j = space()->ei(f,Y,B::Y);
        for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT jj=SW::BL(Y); jj<=0; ++jj )
              temp += getHC(Y,f,j,jj)*x(i,j+jj,k);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Y,f,j,0);
          }
      }

      // normal direction: Z
      if( 0<space()->bcl(Z) ) {
        OT k = space()->si(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT kk=0; kk<=SW::BU(Z); ++kk )
              temp += getHC(Z,f,k,kk)*x(i,j,k+kk);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }
      if( 0<space()->bcu(Z) ) {
        OT k = space()->ei(f,Z,B::Y);
        for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
          for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
            ST temp = 0.;
            for( OT kk=SW::BL(Z); kk<=0; ++kk )
              temp += getHC(Z,f,k,kk)*x(i,j,k+kk);
            y(i,j,k) = x(i,j,k) + omegaBC*( b(i,j,k) - temp )/getHC(Z,f,k,0);
          }
      }
    }
  }

public:

  void apply( const FluxFieldT& wind, const DomainFieldT& x, RangeFieldT& y, const Add add=Add::N ) const {

    const F m = y.getType();

    DomainFieldT temp( space(), true, m );

    assert( y.getType() == x.getType() );

    for( int i =0; i<SpaceT::sdim; ++i )
      assert( wind[i].getType() == x.getType() );

    for( int vel_dir=0; vel_dir<SpaceT::sdim; ++vel_dir )
      wind[vel_dir].exchange();

    for( int i=0; i<nIter_; ++i ) {
      applyStep( wind, x, y, temp );
      applyStep( wind, x, temp, y );
    }
  }

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op_->space();
  };

  void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    op_->print();
  }


  bool hasApplyTranspose() const {
    return false;
  }

  constexpr const std::string getLabel() const {
    return "ConvectionDiffusionJSmoother ";
  };

}; // end of class ConvectionDiffusionJSmoother


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP
