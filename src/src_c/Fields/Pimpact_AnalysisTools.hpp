#pragma once
#ifndef PIMPACT_ANALYSISTOOLS_HPP
#define PIMPACT_ANALYSISTOOLS_HPP


#include <cmath>
#include <vector>

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include "Pimpact_DivOp.hpp"
#include "Pimpact_InterpolateS2VOp.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {


template<class MultiHarmonicFieldT>
typename MultiHarmonicFieldT::SpaceT::Scalar
truncErrorEstimate( const MultiHarmonicFieldT& field, ENorm normType=ENorm::L2 ) {

  using ST = typename MultiHarmonicFieldT::SpaceT::Scalar;
  using OT = typename MultiHarmonicFieldT::SpaceT::Ordinal;

  auto space = field.space();

  ST u_1, u_nf;

  OT nf = space->nGlo(3);


  if( 0<nf and space->si(F::U,3)<=1 and 1<=space->ei(F::U,3) )
    u_1  = field.getField(1).norm(normType);

  if( 0<nf and space->si(F::U,3)<=nf and space->nGlo(3)<=nf )
    u_nf = field.getField(nf).norm(normType);

  int rank_1 = 0;
  if( 1==space->getProcGrid()->getNP(3) )
    rank_1 = 0;
  else if( 0==(space->nGlo(3)+1)%space->getProcGrid()->getNP(3) )
    rank_1 = 1;
  int rank_nf = space->getProcGrid()->getNP(3)-1;

  // nice nonblocking version
  MPI_Request req_1, req_nf;  

  MPI_Ibcast(
      &u_1,                                // buffer	starting address of buffer (choice)
      1,                                   // number of entries in buffer (non-negative integer)
      MPI_DOUBLE,                          // data type of buffer (handle)
      rank_1,                              // rank of broadcast root (integer)
      space->getProcGrid()->getCommBar(3), // communicator (handle)
      &req_1);                             // communication request
  MPI_Ibcast(
      &u_nf,                               // buffer	starting address of buffer (choice)
      1,                                   // number of entries in buffer (non-negative integer)
      MPI_DOUBLE,                          // data type of buffer (handle)
      rank_nf,                             // rank of broadcast root (integer)
      space->getProcGrid()->getCommBar(3), // ccommunicator (handle) ommunicator (handle)
      &req_nf);                            // communication request

  MPI_Wait(&req_1, MPI_STATUS_IGNORE); 
  MPI_Wait(&req_nf, MPI_STATUS_IGNORE); 

  ST truncError = 1.;

  if( u_nf != u_1 ) // just in case u_1=u_nf=0
    truncError = u_nf / u_1 ;

  return truncError;
}


template<class MHFieldT>
void writeMHLambda2( MHFieldT& x, int id ) {

  using SpaceT = typename MHFieldT::SpaceT;

  using OT = typename SpaceT::Ordinal;

  auto space = x.space();

  if( 0==space->si(F::U,3) )
    writeLambda2( x.get0Field(), id );

  for( OT i=std::max(space->si(F::U,3),1); i<=space->ei(F::U,3); ++i ) {
    writeLambda2( x.getCField(i), id+2*i-1 );
    writeLambda2( x.getSField(i), id+2*i   );
  }
}


template<class MHFieldT>
void writeMHLambda2evol( MHFieldT& x, int id ) {

  using SpaceT = typename MHFieldT::SpaceT;

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  auto space = x.space();

  x.exchange();

  if( space()->getProcGrid()->getIB(3)==1 ) {

    ST pi = 4.*std::atan(1.);
    OT nf = space()->nGlo(3);
    OT nt = 4*nf;

    auto temp = x.get0Field().clone( ECopy::Shallow );

    for( OT i=0; i<nt;  ++i ) {
      *temp = x.get0Field();
      for( OT j=1; j<=nf; ++j ) {
        temp->add(
            1., *temp,
            std::sin( 2.*pi*i*(static_cast<ST>(j))/nt ), x.getSField(j) );
        temp->add(
            std::cos( 2.*pi*i*(static_cast<ST>(j))/nt ), x.getCField(j),
            1., *temp );
        writeLambda2( *temp, id+2*i );
      }
    }
  }
}


template<class VectorFieldT>
void writeLambda2( VectorFieldT& x, int id ) {

  using SpaceT = typename VectorFieldT::SpaceT;

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;
  using SFT = ScalarField<SpaceT>;

  auto space = x.space();

  Teuchos::SerialDenseVector<OT, ST> TT(6);

  ST PP, QQ, RR;
  ST rho, theta;

  ST eps = 1.e-5;

  //ST temp;
  Teuchos::SerialDenseVector<OT, ST> lam(3);

  OT m, n;

  //--- pi ---
  ST pi = 2.*std::abs(std::acos(0.));

  //==============================================================================
  //=== interpolate velocities  ==================================================
  //==============================================================================
  {
    auto interpolateV2S = space->getInterpolateV2S();
    auto interpolateS2V = create<InterpolateS2V>( space );

    SFT temp( space );

    for( F i=F::U; i<SpaceT::sdim; ++i ) {
      interpolateV2S->apply( x(i), temp );
      interpolateS2V->apply( temp, x(i) );
    }
  }

  //======================================================================================
  //=== lambda ===========================================================================
  //======================================================================================

  x.exchange();
  ScalarField<SpaceT> res( space );

  auto divOp = create<DivOp>( space );

  for( OT k=space->si(F::S,Z); k<=space->ei(F::S,Z); ++k )
    for( OT j=space->si(F::S,Y); j<=space->ei(F::S,Y); ++j )
      for( OT i=space->si(F::S,X); i<=space->ei(F::S,X); ++i ) {

        Teuchos::SerialDenseMatrix<OT, ST> Dvel(3, 3, true);

        //--- d/dx ----
        for( int ii=space->dl(X); ii<=space->du(X); ++ii ) {
          Dvel(X, X) += divOp->getC(X,i,ii)*x(F::U)(i+ii,j,k);
          Dvel(Y, X) += divOp->getC(X,i,ii)*x(F::V)(i+ii,j,k);
          Dvel(Z, X) += divOp->getC(X,i,ii)*x(F::W)(i+ii,j,k);
        }

        //--- d/dy ------------------------------------------
        for( int jj=space->dl(Y); jj<=space->du(Y); ++jj ) {
          Dvel(X, Y) +=  divOp->getC(Y,j,jj)*x(F::U)(i,j+jj,k);
          Dvel(Y, Y) +=  divOp->getC(Y,j,jj)*x(F::V)(i,j+jj,k);
          Dvel(Z, Y) +=  divOp->getC(Y,j,jj)*x(F::W)(i,j+jj,k);
        }

        //--- d/dz ------------------------------------------
        for( int kk=space->dl(Z); kk<=space->du(Z); ++kk ) {
          Dvel(X, Z) +=  divOp->getC(Z,k,kk)*x(F::U)(i,j,k+kk);
          Dvel(Y, Z) +=  divOp->getC(Z,k,kk)*x(F::V)(i,j,k+kk);
          Dvel(Z, Z) +=  divOp->getC(Z,k,kk)*x(F::W)(i,j,k+kk);
        }


        //--- Tensor symmetric ----------------------------------------------------
        TT(X) = 2.*(Dvel(X,X)*Dvel(X,X) + Dvel(X,Y)*Dvel(Y,X) + Dvel(X,Z)*Dvel(Z,X));
        TT(Y) = 2.*(Dvel(Y,X)*Dvel(X,Y) + Dvel(Y,Y)*Dvel(Y,Y) + Dvel(Y,Z)*Dvel(Z,Y));
        TT(Z) = 2.*(Dvel(Z,X)*Dvel(X,Z) + Dvel(Z,Y)*Dvel(Y,Z) + Dvel(Z,Z)*Dvel(Z,Z));

        TT(3) = (Dvel(X,X) + Dvel(Y,Y))*(Dvel(X,Y) + Dvel(Y,X)) + Dvel(X,Z)*Dvel(Z,Y) + Dvel(Y,Z)*Dvel(Z,X);
        TT(4) = (Dvel(X,X) + Dvel(Z,Z))*(Dvel(X,Z) + Dvel(Z,X)) + Dvel(X,Y)*Dvel(Y,Z) + Dvel(Z,Y)*Dvel(Y,X);
        TT(5) = (Dvel(Y,Y) + Dvel(Z,Z))*(Dvel(Y,Z) + Dvel(Z,Y)) + Dvel(Y,X)*Dvel(X,Z) + Dvel(Z,X)*Dvel(X,Y);

        //--- Invariants -------------------------------------------------------------
        PP = TT(X) + TT(Y) + TT(Z);
        QQ = TT(X)*TT(Y) + TT(X)*TT(Z) + TT(Y)*TT(Z) - std::pow(TT(3), 2) -
          std::pow(TT(4), 2) - std::pow(TT(5), 2);
        RR = TT(X)*TT(Y)*TT(Z) + 2.*TT(3)*TT(4)*TT(5) - std::pow(TT(X)*TT(5), 2) -
          TT(Y)*std::pow(TT(4), 2) - TT(Z)*std::pow(TT(3), 2);


        //--- Eigenwerte --------------------------------------------------------------
        rho = std::pow((PP/3.), 2) - QQ/3.;
        if( rho <= eps ) {
          //-------------------------------------------------------------------------
          // y:=lam-PP/3.                                                            
          // y**3+p*y+q=0.                                                           
          // p:=QQ-PP**2/3.                                                          
          // q:=-2*(PP/3)**3+PP*QQ/3.-RR                                             
          // TT ist symmetrisch ==> lam(1:3) reell <==> D=(p/3.)**3+(q/2.)**2 .LE. 0.
          //                    ==> falls p=QQ-PP**2/3.=0.                           
          //                    ==> q:=-2*(PP/3)**3+PP**3/9.-RR = (PP/3)**3-RR      
          //                    ==> D=(q/2.)**2 .LE. 0. <==> q=0. <==> RR=(PP/3)**3
          //                    <==> lam(1:3)=lam(3) 
          //                    ==> y=0.                                             
          //                    ==> lam=PP/3.=RR**(1./3.)                             
          //--------------------------------------------------------------------------
          res(i,j,k) = PP/3.;
        }
        else {
          rho = std::sqrt( rho );
          QQ  = (std::pow(PP/3., 3) - QQ*(PP/6.) + RR/2.)/std::pow(rho, 3);

          if( std::abs(QQ) < 1.) 
            theta = std::acos(QQ)/3.;
          else {
            if( QQ > 0.) 
              theta = 0.;
            else
              theta = pi/3.;
          }

          lam(X) = std::cos(theta           );
          lam(Y) = std::cos(theta + 2.*pi/3.);
          lam(Z) = std::cos(theta + 4.*pi/3.);


          //--- sortieren ------------------------------------------------------------
          std::sort( lam.values(), lam.values()+3 );
          //for( int m=0; m<3; ++m )
            //for( int n=m+1; n<3; ++n )
              //if( lam(m) > lam(n) )
                //std::swap( lam(m), lam(n) );

          // Faktor 1/2, da bei TT(1:6) mit 2 durchmultipliziert wurde ...
          res(i,j,k) = rho*lam(2) + PP/6.;
        }
      }
  res.write( id );
}


/// \note has not a very good parallel efficiency( exchange) so only use for analyzing
template<class MultiHarmonicFieldT>
void writeSpectrum( const MultiHarmonicFieldT& field,
    std::ostream& out=std::cout ) {

  using OT = typename MultiHarmonicFieldT::SpaceT::Ordinal;

  auto space = field.space();

  // making sure one works on global operator, could be make more efficient to do that
  // just on rank "zero"
  Teuchos::RCP<const MultiHarmonicFieldT> y;

  if( field.global()==MultiHarmonicFieldT::Global::Y )
    y = Teuchos::rcpFromRef( field );
  else {
    Teuchos::RCP<MultiHarmonicFieldT> temp =
      Teuchos::rcp( new MultiHarmonicFieldT( space, MultiHarmonicFieldT::Global::Y ) );
    *temp = field;
    y = temp;
  }
  y->exchange();


  out << 0 << "\t" << y->get0Field().norm(ENorm::L2)*std::sqrt(2.) << "\n";
  for( OT i=1; i<=space->nGlo(3); ++i )
    out << i << "\t" << y->getField(i).norm(ENorm::L2) << "\n";
}


/// - compute energy local
/// - reduce to rank_dim=0
/// - having two vectors for y and e(z)
/// \f[ e(dir) = \int \int (u^2 + v^2 + z^2 ) \exp( -\frac{x}{\gamma}^2/2 ) \f]
template<class SpaceT>
void computeEnergyDir( const VectorField<SpaceT>& vel, std::ostream& out=std::cout,
    const ECoord dir=ECoord::Y, const typename SpaceT::Scalar gamma=0. ) {

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  auto space = vel.space();
  auto coord = space->getCoordinatesLocal();

  ScalarField<SpaceT> temp(space);

  std::vector<ST> energy( space->nLoc(dir), Teuchos::ScalarTraits<ST>::zero() );

  for( F f=F::U; f<SpaceT::sdim; ++f ) {

    space->getInterpolateV2S()->apply( vel(f), temp );

    for( OT k=space->si(F::S,Z); k<=space->ei(F::S,Z); ++k )
      for( OT j=space->si(F::S,Y); j<=space->ei(F::S,Y); ++j )
        for( OT i=space->si(F::S,X); i<=space->ei(F::S,X); ++i ) {
          ST vel = temp(i,j,k);
          switch( dir ) {
            case X: {
              ST volume = coord->dx(F::S,Y,j) * coord->dx(F::S,Z,k);
              energy[i-space->si(F::S,X)] +=
                std::pow(volume * vel * ((0.==gamma)?1.:std::exp(-0.5*std::pow( coord->getX(F::S,Z,k)/gamma, 2))), 2);
              break;
            }
            case Y: {
              ST z = coord->getX(F::S, Z, k);
              ST dz = coord->dx(F::S, Z, k);
              ST dx = coord->dx(F::S,X,i);
              energy[j-space->si(F::S,Y)] +=
                dx * std::pow( dz * vel * std::exp(-0.5*std::pow(z/gamma, 2)), 2);
              break;
            }
            case Z: {
              ST volume = coord->dx(F::S,X,i) * coord->dx(F::S,Y,j);
              energy[k-space->si(F::S,Z)] +=
                std::pow( volume * vel * ((0.==gamma)?1.:std::exp(-0.5*std::pow( coord->getX(F::S,Z,k)/gamma, 2))), 2);
              break;
            }
            case T: {
              std::cout << "Warning: not implemented!!!\n";
              break;
            }
          }
        }
  }

  MPI_Reduce(
      (0==space->getProcGrid()->getRankSlice(dir))?
      MPI_IN_PLACE:energy.data(),	                // void* send_data,
      energy.data(),                              // void* recv_data,
      space->nLoc(dir),                           // int count,
      MPI_REAL8,                                  // MPI_Datatype datatype,
      MPI_SUM,                                    // MPI_Op op,
      0,                                          // int root,
      space->getProcGrid()->getCommSlice(dir) );  // MPI_Comm communicator);

  //std::cout << space->nGlo(dir) << "\n";
  //std::cout << space->nLoc(dir) << "\n";
  //std::cout << space->getProcGrid()->getNP(dir) << "\n";

  std::vector<ST> energyGlobal( space->nGlo(dir), Teuchos::ScalarTraits<ST>::zero() );

  if( 0==space->getProcGrid()->getRankSlice(dir) ) {
    MPI_Gather(
        energy.data(),                           // void* send_data,
        space->nLoc(dir)-1,                      // int send_count,
        MPI_REAL8,                               // MPI_Datatype send_datatype,
        energyGlobal.data(),                     // void* recv_data,
        space->nLoc(dir)-1,                      // int recv_count,
        MPI_REAL8,                               // MPI_Datatype recv_datatype,
        0,                                       // int root,
        space->getProcGrid()->getCommBar(dir) ); // MPI_Comm communicator);

    if( 0==space->getProcGrid()->getRankBar(dir) )
      for( OT j=0; j<space->nGlo(dir); ++j )
        out << space->getCoordinatesGlobal()->getX(F::S,dir,j+1) << "\t" << energyGlobal[j] << "\n";
  }
}



/// - compute energy local
/// - reduce to rank_dim=0
/// - having two vectors for y and e(z)
/// \f[ e(dir) = \int \int (u^2 + v^2 + z^2 ) \exp( -\frac{x}{\gamma}^2/2 ) \f]
template<class SpaceT>
void computeHEEnergyDir(
    const VectorField<SpaceT>& vel,
    std::ostream& out=std::cout,
    const typename SpaceT::Scalar gamma=10.,
    const int n_Hermitemodes=4 ) {

  const ECoord dir=ECoord::Y;

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  auto space = vel.space();
  auto coord = space->getCoordinatesLocal();

  //----- obristd 310510: initialization of spatial modal analysis
  if( 0==space->rankST() ) {
    std::cout << "\nFourier-Hermite modal analysis\n";
    std::cout << "==============================\n";
    std::cout << " gamma     = " << gamma << "\n";
    std::cout << " n_Hermite = " << n_Hermitemodes << "\n";
  }

  OT S1p = space->si(F::S,X);
  OT N1p = space->ei(F::S,X);
  OT S2p = space->si(F::S,Y);
  OT N2p = space->ei(F::S,Y);
  OT S3p = space->si(F::S,Z);
  OT N3p = space->ei(F::S,Z);

  ST pi = 4.*std::atan(1.);

  Teuchos::SerialDenseMatrix<OT, ST> sweight(N3p-S3p+1, n_Hermitemodes);

  for( OT k=S3p; k<=N3p; ++k ) {

    ST x3p = coord->getX(F::S,Z,k);

    ST dx3p = coord->dx(F::S,Z,k);

    // He0
    sweight(k-S3p,0) = 1.;

    // He1
    if( n_Hermitemodes>1 ) 
      sweight(k-S3p,1) = x3p/gamma;

    // He2, He3, ...
    for( int n=2; n<n_Hermitemodes; ++n )
      sweight(k-S3p,n) = (x3p/gamma*sweight(k-S3p,n-1) - std::sqrt(static_cast<ST>(n-1))*sweight(k-S3p,n-2))/sqrt(static_cast<ST>(n));

    //if( 0==space->rankST() ) std::cout << "k: " << k << " x3p: " << x3p << "\t";
    for( int n=0; n<n_Hermitemodes; ++n ) {
      sweight(k-S3p,n) *= std::exp( -std::pow(x3p/gamma, 2)/2. )*dx3p/std::sqrt(gamma*std::sqrt(2*pi)); // dh: don't know where it is coming from: normalization factor for exp(...)
      //sweight(k-S3p,n) *= std::exp( -std::pow(x3p/gamma, 2)/2. )*dx3p;
      //if( 0==space->rankST() ) std::cout << sweight(k-S3p, n) << "\t";
    }
    //std::cout << "\n";
  }


  ScalarField<SpaceT> vel_dist0(space);
  ScalarField<SpaceT> vel_dist1(space);
  ScalarField<SpaceT> vel_dist2(space);

  ScalarField<SpaceT> vel_mode_0(space);
  ScalarField<SpaceT> vel_mode_1(space);
  ScalarField<SpaceT> vel_mode_2(space);

  space->getInterpolateV2S()->apply( vel(F::U), vel_dist0 );
  space->getInterpolateV2S()->apply( vel(F::V), vel_dist1 );
  space->getInterpolateV2S()->apply( vel(F::W), vel_dist2 );


  Teuchos::SerialDenseMatrix<OT, ST> He(3, n_Hermitemodes);
  Teuchos::SerialDenseMatrix<OT, ST> He_global(3, n_Hermitemodes);
  Teuchos::SerialDenseMatrix<OT, ST> energydensity( n_Hermitemodes, space->nLoc(dir) );

  for( OT j=S2p; j<=N2p; ++j ) {        // spanwise
    for( OT i=S1p; i<=N1p; ++i) {       // wall-normal

      // --- integral over z-dimension: \int_{-L_z}^{L_z} u' He(z/gamma) exp( -(z/gamma)^2 dz ---
      // ---------------- HERMITE TRANSFORM OF ALL THREE VELOCITY DISTURBANCES
      He.putScalar( 0. );                            //  initialize Hermite coefficients
      for( OT k=S3p; k<=N3p; ++k ) {  //  chordwise
        // chordwise Hermite modes by discrete He transform
        // initialization of sweight and to be found in usr_initcond.f90 at line 100
        for( OT n = 0; n<n_Hermitemodes; ++n ) {
          He(0,n) += sweight(k-S3p,n)*vel_dist0(i,j,k);
          He(1,n) += sweight(k-S3p,n)*vel_dist1(i,j,k);
          He(2,n) += sweight(k-S3p,n)*vel_dist2(i,j,k);
        }
      }
      // sum up all contributions to He along x3-direction
      // call mpi_allreduce(He, He_global, 3*n_Hermitemodes, MPI_REAL8, MPI_SUM, COMM_BAR3, merror)
      MPI_Allreduce( He.values(), He_global.values(), 3*n_Hermitemodes, MPI_REAL8, MPI_SUM, space->getProcGrid()->getCommBar(Z) );

      // --- integral over x-dimension: \int_0^\infty u_kn \dx
      //ST dx1p = coord->dx(F::S,X,i);
      //for( OT n = 0; n<n_Hermitemodes; ++n ) {
      //energydensity( n, j-space->si(F::S,Y) ) +=
      //(std::pow(He_global(0,n), 2) +
      //std::pow(He_global(1,n), 2) +
      //std::pow(He_global(2,n), 2) )*dx1p;
      //}
      //// --- integral over x-dimension: \int_0^\infty u_kn \dx with weird minus one degree
      ST dx1p = coord->dx(F::S,X,i);
      energydensity( 0, j-space->si(F::S,Y) ) +=
        std::pow(He_global(2,0), 2)*dx1p;
      for( OT n = 1; n<n_Hermitemodes; ++n ) {
        energydensity( n, j-space->si(F::S,Y) ) +=
          (std::pow(He_global(0,n-1), 2) +
           std::pow(He_global(1,n-1), 2) +
           std::pow(He_global(2,n), 2) )*dx1p;
      }
    }
  }


  Teuchos::SerialDenseMatrix<OT, ST> energyGlobal( n_Hermitemodes, space->nGlo(dir) );

  if( 0==space->getProcGrid()->getRankSlice(dir) ) {
    MPI_Gather(
        energydensity.values(),                  // void* send_data,
        (space->nLoc(dir)-1)*n_Hermitemodes,     // int send_count,
        MPI_REAL8,                               // MPI_Datatype send_datatype,
        energyGlobal.values(),                     // void* recv_data,
        (space->nLoc(dir)-1)*n_Hermitemodes,     // int recv_count,
        MPI_REAL8,                               // MPI_Datatype recv_datatype,
        0,                                       // int root,
        space->getProcGrid()->getCommBar(dir) ); // MPI_Comm communicator);

    if( 0==space->getProcGrid()->getRankBar(dir) )
      for( OT j=0; j<space->nGlo(dir); ++j ) {
        out << space->getCoordinatesGlobal()->getX(F::S,dir,j+1) << "\t" ;
        for( OT n=0; n<n_Hermitemodes; ++n ) 
          out << energyGlobal(n, j) << "\t";
        out << "\n";
      }
  }
}
} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_ANALYSISTOOLS_HPP
