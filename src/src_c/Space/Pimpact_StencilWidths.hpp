#pragma once
#ifndef PIMPACT_STENCILWIDTHS_HPP
#define PIMPACT_STENCILWIDTHS_HPP


#include <iostream>

#include "mpi.h"

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"




namespace Pimpact {


/// \brief contains the dimesnion and bounds of the different stencils.
///
/// there are three different kind of stencil the central ones for helmholtz, than the one for divergence
/// and the gradient like.
/// \ingroup SpaceObject
/// \tparam dim dimension of grid can be 3 or 4
/// \tparam dimNC dimension of stencil
/// - 4: Stabil   (xi >= 2, Re=10000, N=17)
/// - 5: Stabil  (Re=10000, N=65, leicht gestreckt, explizites Forcing)
/// - 6: Instabil (äquidistant, Re=10000, N=17)
template< int dim, int dimNC >
class StencilWidths {

  template< int d, int dnc >
  friend const Teuchos::RCP<const StencilWidths<d,dnc> > createStencilWidths(const bool&);

protected:

	const bool spectralT_;

  typedef const Teuchos::Tuple<int,dim> TO;

  typedef  Teuchos::Tuple<  Teuchos::Tuple<int,dimNC>,3> TStenc;

  /// \{ \name Konvergenzordnung der Differenzenkoeffizienten (Anzahl Koeffizienten)
  /// \{ \brief Anzahl Stencil-Koeffizienten (Rand)
  TStenc ncbC_;
  TStenc ncbD_;
  TStenc ncbG_;
  /// \}
  /// \}

  /// \{ \name Intervallgrenzen der Differenzen-Koeffizienten-Arrays

  /// \{ \brief Central
  TO bl_;
  TO bu_;
  /// \}

  /// \{ \brief Divergenz
  TO dl_;
  TO du_;
  /// \}

  /// \{ \brief Gradient
  TO gl_;
  TO gu_;
  /// \}


  /// \{ \brief upwind (nicht-linear)
  ///
  /// (aktuell wird nicht zwischen auf- und abwärtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  ///  zu können, wo KEINE upwind-Differenzierung verwendet wird)
  TO nl_;
  TO nu_;
  /// \}

  /// \}

  Teuchos::Tuple<int,3> ls_; /// < \brief Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup)

protected:

  /// \brief constructor
  StencilWidths( const bool& spectralT ):
		spectralT_(spectralT),
    ncbC_(),
    ncbD_(),
    ncbG_(),
    bl_(),
    bu_(),
    dl_(),
    du_(),
    gl_(),
    gu_(),
    nl_(),
    nu_(),
    ls_(Teuchos::tuple(-1,-1,-1) ) {

    //--- Anzahl Stencil-Koeffizienten (Rand) -------------------------------------------------------------------
    // implementation here is a little bit fuzzy, because Tuple has no nice templated constructor we copy
    for( int i=0; i<3; ++i ) {
      switch(dimNC) {
      case(2): {

        auto tempC = Teuchos::tuple( 2, 3 );
        auto tempD = Teuchos::tuple( 2, 2 );
        auto tempG = Teuchos::tuple( 0, 2 );

        for( int j=0; j<dimNC; ++j ) {
          ncbC_[i][j] = tempC[j];
          ncbD_[i][j] = tempD[j];
          ncbG_[i][j] = tempG[j];
        }
        break;
      }
      case(3): {
            auto tempC = Teuchos::tuple( 3, 4, 5 );
            auto tempD = Teuchos::tuple( 3, 4, 4 );
            auto tempG = Teuchos::tuple( 2, 3, 4 );

        for( int j=0; j<dimNC; ++j ) {
          ncbC_[i][j] = tempC[j];
          ncbD_[i][j] = tempD[j];
          ncbG_[i][j] = tempG[j];
        }
        break;
      }
      case(4): {
        // Stabil   (xi >= 2, Re=10000, N=17)
        auto tempC = Teuchos::tuple( 4, 5, 5, 7 );
        auto tempD = Teuchos::tuple( 4, 4, 6, 6 );
        auto tempG = Teuchos::tuple( 3, 4, 4, 6 );
        for( int j=0; j<dimNC; ++j ) {
          ncbC_[i][j] = tempC[j];
          ncbD_[i][j] = tempD[j];
          ncbG_[i][j] = tempG[j];
        }
        break;
      }
      case(5): {
        // Instabil (äquidistant, Re=10000, N=17)
        // auto tempC = Teuchos::tuple( 5, 6, 6, 7, 9 );
        // auto tempD = Teuchos::tuple( 5, 4, 6, 8, 8 );
        // auto tempG = Teuchos::tuple( 0, 5, 4, 6, 8 );
        // Stabil  (Re=10000, N=65, leicht gestreckt, explizites Forcing)
        auto tempC = Teuchos::tuple( 3, 7, 7, 7, 9 );
        auto tempD = Teuchos::tuple( 6, 6, 6, 8, 8 );
        auto tempG = Teuchos::tuple( 0, 6, 6, 6, 8 );

        for( int j=0; j<dimNC; ++j ) {
          ncbC_[i][j] = tempC[j];
          ncbD_[i][j] = tempD[j];
          ncbG_[i][j] = tempG[j];
        }
        break;
      }
      case(6): {
        // Instabil (äquidistant, Re=10000, N=17)
        auto tempC = Teuchos::tuple( 6, 7, 7, 7,  9, 11 );
        auto tempD = Teuchos::tuple( 6, 6, 6, 8, 10, 10 );
        auto tempG = Teuchos::tuple( 5, 6, 6, 6,  8, 10 );
        for( int j=0; j<dimNC; ++j ) {
          ncbC_[i][j] = tempC[j];
          ncbD_[i][j] = tempD[j];
          ncbG_[i][j] = tempG[j];
        }
        break;
      }
      default:
        // throw exeption
        ;
      }

    }

    // Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
    Teuchos::Tuple<int,3> ncC;
    Teuchos::Tuple<int,3> ncS;

    for( int i=0; i<3; ++i ) {
      ncC[i] = ncbC_[i][dimNC-1];
      ncS[i] = ncbG_[i][dimNC-1];
    }

    // zentral
    for( int i=0; i<3; ++i ) {
      bu_[i] = ncS[i]/2;
      bl_[i] = -bu_[i];
    }
    if( 4==dim ) bl_[3] = -1;
		if( 4==dim ) bu_[3] =  1;
//    if( 4==dim ) bu_[3] =  0;

    // divergence
    for( int i=0; i<3; ++i ) {
      dl_[i] = bl_[i];
      du_[i] = bu_[i]-1;
    }
    if( 4==dim ) dl_[3] = bl_[3];
    if( 4==dim ) du_[3] = bu_[3];

    // gradient
    for( int i=0; i<3; ++i ) {
      gl_[i] = bl_[i]+1;
      gu_[i] = bu_[i];
    }
    if( 4==dim ) gl_[3] = bl_[3];
    if( 4==dim ) gu_[3] = bu_[3];

    // upwind (nicht-linear)
    for( int i=0; i<3; ++i ) {
      nl_[i] = bl_[i];
      nu_[i] = bu_[i];
    }
    if( 4==dim ) nl_[3] = bl_[3];
    if( 4==dim ) nu_[3] = bu_[3];

};

  /// \brief constructor
  ///
  /// \param bl lower bound of storage
  /// \param bu upper bound of storage
  StencilWidths(
      //      int dimension,
      TO bl,
      TO bu ):
        bl_(bl),
        bu_(bu),
        ls_(Teuchos::tuple(-1,-1,-1))
  {};

public:

  /// prints to \c std::cout, only for debuging purpose
  void print( std::ostream& out=std::cout ) const {
    out << "\t---StencilWidths: ---\n";
    out << "\tcomput dim: " << dim << "\n";
    out << "\tncbC: " << ncbC_ << "\n";
    out << "\tncbD: " << ncbD_ << "\n";
    out << "\tncbG: " << ncbG_ << "\n";
    out << "\tbl: " << bl_ << "\n";
    out << "\tbu: " << bu_ << "\n";
    out << "\tdl: " << dl_ << "\n";
    out << "\tdu: " << du_ << "\n";
    out << "\tgl: " << gl_ << "\n";
    out << "\tgu: " << gu_ << "\n";
    out << "\tnl: " << nl_ << "\n";
    out << "\tnu: " << nu_ << "\n";
    out << "\tls: " << ls_ << "\n";
  }

	const bool& spectralT() const { return( spectralT_ ); }

  int  getDimNcbC( int i ) const { return( ncbC_[i].size() ); }
  int  getDimNcbD( int i ) const { return( ncbD_[i].size() ); }
  int  getDimNcbG( int i ) const { return( ncbG_[i].size() ); }

//  const int&  getDimNcbC( int i ) const { return( dimNC ); }
//  const int&  getDimNcbD( int i ) const { return( dimNC ); }
//  const int&  getDimNcbG( int i ) const { return( dimNC ); }

  const int* getNcbC( int i ) const { return( ncbC_[i].getRawPtr() ); }
  const int* getNcbD( int i ) const { return( ncbD_[i].getRawPtr() ); }
  const int* getNcbG( int i ) const { return( ncbG_[i].getRawPtr() ); }

  const int* getBL()        const { return( bl_.getRawPtr() ); }
  const int& getBL( int i ) const { return( bl_[i] ); }

  const int* getBU()        const { return( bu_.getRawPtr() ); }
  const int& getBU( int i ) const { return( bu_[i] ); }

  const int* getDL()        const { return( dl_.getRawPtr() ); }
  const int& getDL( int i ) const { return( dl_[i] ); }

  const int* getDU()        const { return( du_.getRawPtr() ); }
  const int& getDU( int i ) const { return( du_[i] ); }

  const int* getGL()        const { return( gl_.getRawPtr() ); }
  const int& getGL( int i ) const { return( gl_[i] ); }

  const int* getGU()        const { return( gu_.getRawPtr() ); }
  const int& getGU( int i ) const { return( gu_[i] ); }

  const int* getNL()        const { return( nl_.getRawPtr() ); }
  const int& getNL( int i ) const { return( nl_[i] ); }

  const int* getNU()        const { return( nu_.getRawPtr() ); }
  const int& getNU( int i ) const { return( nu_[i] ); }

  const int* getLS()        const { return( ls_.getRawPtr() ); }
  const int& getLS( int i ) const { return( ls_[i] ); }

}; // end of class StencilWidths



/// \brief function that creates \c Pimpact:StencilWidths
/// by getting values from \c IMPACT
/// should be changed by getting vallues from \c ProcGridSize and \c GridSize
/// \relates StencilWidths
template< int d, int dnc  >
const Teuchos::RCP<const StencilWidths<d,dnc> > createStencilWidths( const bool& spectralT ){

  return(
      Teuchos::RCP<const StencilWidths<d,dnc> > (
          new StencilWidths<d,dnc>( spectralT ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_STENCILWIDTHS_HPP
