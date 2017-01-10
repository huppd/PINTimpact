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
/// there are three different kind of stencil the central ones for Helmholtz, than the one for divergence
/// and the gradient like.
/// \ingroup SpaceObject
/// \tparam dim dimension of grid can be 3 or 4 as soon as time is own class  ->sdim
/// \tparam dimNC dimension of stencil
/// - 4: Stabil   (xi >= 2, Re=10000, N=17)
/// - 5: Stabil  (Re=10000, N=65, leicht gestreckt, explizites Forcing)
/// - 6: Instabil (äquidistant, Re=10000, N=17)
/// \todo rm dim, make it same in every direction
template< int dim, int dimNC >
class StencilWidths {

  template< int d, int dnc >
  friend const Teuchos::RCP<const StencilWidths<d,dnc> > createStencilWidths(const bool&);

protected:

	const bool spectralT_;

  using TO = const Teuchos::Tuple<int,dim>;

  using TStenc =  const Teuchos::Tuple<  Teuchos::Tuple<int,dimNC>,3>;

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

					Teuchos::Tuple<int,2> tempC = Teuchos::tuple( 2, 3 );
					Teuchos::Tuple<int,2>	tempD = Teuchos::tuple( 2, 2 );
					Teuchos::Tuple<int,2>	tempG = Teuchos::tuple( 0, 2 );

					for( int j=0; j<dimNC; ++j ) {
						ncbC_[i][j] = tempC[j];
						ncbD_[i][j] = tempD[j];
						ncbG_[i][j] = tempG[j];
					}
					break;
			 }
				case(3): {
					 Teuchos::Tuple<int,3> tempC = Teuchos::tuple( 3, 4, 5 );
					 Teuchos::Tuple<int,3> tempD = Teuchos::tuple( 3, 4, 4 );
					 Teuchos::Tuple<int,3> tempG = Teuchos::tuple( 2, 3, 4 );

					for( int j=0; j<dimNC; ++j ) {
						ncbC_[i][j] = tempC[j];
						ncbD_[i][j] = tempD[j];
						ncbG_[i][j] = tempG[j];
					}
					break;
				}
				case(4): {
					// Stabil   (xi >= 2, Re=10000, N=17)
				  Teuchos::Tuple<int,4> tempC = Teuchos::tuple( 4, 5, 5, 7 );
				  Teuchos::Tuple<int,4> tempD = Teuchos::tuple( 4, 4, 6, 6 );
				  Teuchos::Tuple<int,4> tempG = Teuchos::tuple( 3, 4, 4, 6 );
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
					Teuchos::Tuple<int,6> tempC = Teuchos::tuple( 6, 7, 7, 7,  9, 11 );
					Teuchos::Tuple<int,6> tempD = Teuchos::tuple( 6, 6, 6, 8, 10, 10 );
					Teuchos::Tuple<int,6> tempG = Teuchos::tuple( 5, 6, 6, 6,  8, 10 );
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
    //Teuchos::Tuple<int,3> ncC;
    Teuchos::Tuple<int,3> ncS;

    for( int i=0; i<3; ++i ) {
      //ncC[i] = ncbC_[i][dimNC-1];
      ncS[i] = ncbG_[i][dimNC-1];
    }

    // zentral
    for( int i=0; i<3; ++i ) {
      bu_[i] = ncS[i]/2;
      bl_[i] = -bu_[i];
    }
    if( 4==dim ) bl_[3] = -1;
		if( 4==dim ) bu_[3] =  1;

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

  static constexpr 
	int BL( const int& i ) { 
		return( (4==1)?-1:-dimNC+1 );
	}
  static constexpr 
	int BU( const int& i ) { 
		return( (4==1)?1:dimNC-1 );
	}
  static constexpr 
	int DL( const int& i ) { 
		return( (4==1)?-1:-dimNC+1 );
	}
  static constexpr 
	int DU( const int& i ) { 
		return( (4==1)?1:dimNC-2 );
	}
  static constexpr 
	int GL( const int& i ) { 
		return( (4==1)?-1:-dimNC+2 );
	}
  static constexpr 
	int GU( const int& i ) { 
		return( (4==1)?1:dimNC-1 );
	}
  static constexpr 
	int NL( const int& i ) { 
		return( (4==1)?-1:-dimNC+1 );
	}
  static constexpr 
	int NU( const int& i ) { 
		return( (4==1)?1:dimNC-1 );
	}

	constexpr const bool& spectralT() const { return( spectralT_ ); }

  constexpr int getDimNcbC( const int& i ) const { return( ncbC_[i].size() ); }
  constexpr int getDimNcbD( const int& i ) const { return( ncbD_[i].size() ); }
  constexpr int getDimNcbG( const int& i ) const { return( ncbG_[i].size() ); }

  constexpr const int* getNcbC( const int& i ) const { return( ncbC_[i].getRawPtr() ); }
  constexpr const int* getNcbD( const int& i ) const { return( ncbD_[i].getRawPtr() ); }
  constexpr const int* getNcbG( const int& i ) const { return( ncbG_[i].getRawPtr() ); }

  constexpr const int* getBL()               const { return( bl_.getRawPtr() ); }
	constexpr const int& getBL( const int& i ) const { return( bl_[i] ); }

  constexpr const int* getBU()               const { return( bu_.getRawPtr() ); }
  constexpr const int& getBU( const int& i ) const { return( bu_[i] ); }

  constexpr const int* getDL()               const { return( dl_.getRawPtr() ); }
  constexpr const int& getDL( const int& i ) const { return( dl_[i] ); }
  constexpr const TO&  getDLTuple()          const { return( dl_ ); }

  constexpr const int* getDU()               const { return( du_.getRawPtr() ); }
  constexpr const int& getDU( const int& i ) const { return( du_[i] ); }
  constexpr const TO&  getDUTuple()          const { return( du_ ); }

  constexpr const int* getGL()               const { return( gl_.getRawPtr() ); }
  constexpr const int& getGL( const int& i ) const { return( gl_[i] ); }

  constexpr const int* getGU()               const { return( gu_.getRawPtr() ); }
  constexpr const int& getGU( const int& i ) const { return( gu_[i] ); }

  constexpr const int* getNL()               const { return( nl_.getRawPtr() ); }
  constexpr const int& getNL( const int& i ) const { return( nl_[i] ); }

  constexpr const int* getNU()               const { return( nu_.getRawPtr() ); }
  constexpr const int& getNU( const int& i ) const { return( nu_[i] ); }

  constexpr const int* getLS()               const { return( ls_.getRawPtr() ); }
  constexpr const int& getLS( const int& i ) const { return( ls_[i] ); }

}; // end of class StencilWidths



/// \brief function that creates \c Pimpact:StencilWidths
/// by getting values from \c IMPACT
/// should be changed by getting vallues from \c ProcGridSize and \c GridSize
/// \relates StencilWidths
template< int d, int dnc  >
const Teuchos::RCP<const StencilWidths<d,dnc> > createStencilWidths( const bool& spectralT ){

	return( Teuchos::RCP<const StencilWidths<d,dnc> > (
				new StencilWidths<d,dnc>( spectralT ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_STENCILWIDTHS_HPP
