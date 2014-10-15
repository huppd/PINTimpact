#pragma once
#ifndef PIMPACT_FIELDSPACE_HPP
#define PIMPACT_FIELDSPACE_HPP

#include "mpi.h"

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"

#include <iostream>



namespace Pimpact {


/// \todo refector name to something like Stencil stuff
template< class Ordinal=int, int dim=3>
class FieldSpace {

  template< class OT, int dT >
  friend const Teuchos::RCP<const FieldSpace<OT,dT> > createFieldSpace();

public:

  typedef const Teuchos::Tuple<Ordinal,dim> TO;


protected:
//  static const int dim_nc = 2;
//  static const int dim_nc = 3;
  static const int dim_nc = 4;

  typedef  Teuchos::Tuple<  Teuchos::Tuple<int,dim_nc>,3> TStenc;
  /// \{ \name Konvergenzordnung der Differenzenkoeffizienten (Anzahl Koeffizienten)
  /// \{ \brief Anzahl Stencil-Koeffizienten (Rand)

  TStenc ncbC_;
  TStenc ncbD_;
  TStenc ncbG_;

  ///\}
  ///\}

  /// \{ \name Intervallgrenzen der Differenzen-Koeffizienten-Arrays

  /// \{ \brief central
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

  Teuchos::Tuple<Ordinal,3> ls_; /// < \brief Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup)

protected:

  /// \brief constructor
  ///
  /// \param bl lower bound of storage
  /// \param bu upper bound of storage
  FieldSpace():
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
    ls_(Teuchos::tuple(-1,-1,-1))
    {

    for( int i=0; i<3; ++i ) {
      //--- Anzahl Stencil-Koeffizienten (Rand) -------------------------------------------------------------------
//      // note dim_nc has to be 2
//      ncbC_[i] = Teuchos::tuple( 2, 3 );
//      ncbD_[i] = Teuchos::tuple( 2, 2 );
//      ncbG_[i] = Teuchos::tuple( 0, 2 );
//
//      // note dim_nc has to be 3
//      ncbC_[i] = Teuchos::tuple( 3, 4, 5 );
//      ncbD_[i] = Teuchos::tuple( 3, 4, 4 );
//      ncbG_[i] = Teuchos::tuple( 2, 3, 4 );

      // Stabil   (xi >= 2, Re=10000, N=17)
      // note dim_nc has to be 4
      ncbC_[i] = Teuchos::tuple( 4, 5, 5, 7 );
      ncbD_[i] = Teuchos::tuple( 4, 4, 6, 6 );
      ncbG_[i] = Teuchos::tuple( 3, 4, 4, 6 );

//      // Instabil (äquidistant, Re=10000, N=17)
//      // note dim_nc has to be 5
//      ncbC_[i] = Teuchos::tuple( 5, 6, 6, 7, 9 );
//      ncbD_[i] = Teuchos::tuple( 5, 4, 6, 8, 8 );
//      ncbG_[i] = Teuchos::tuple( 0, 5, 4, 6, 8 );
//
//      // Instabil (äquidistant, Re=10000, N=17)
//      // note dim_nc has to be 6
//      ncbC_[i] = Teuchos::tuple( 6, 7, 7, 7,  9, 11 );
//      ncbD_[i] = Teuchos::tuple( 6, 6, 6, 8, 10, 10 );
//      ncbG_[i] = Teuchos::tuple( 5, 6, 6, 6,  8, 10 );
//
//      // Stabil  (Re=10000, N=65, leicht gestreckt, explizites Forcing)
//      // note dim_nc has to be 5
//      ncbC_[i] = Teuchos::tuple( 3, 7, 7, 7, 9 );
//      ncbD_[i] = Teuchos::tuple( 6, 6, 6, 8, 8 );
//      ncbG_[i] = Teuchos::tuple( 0, 6, 6, 6, 8 );

    }

    // Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
    Teuchos::Tuple<int,3> ncC;
    Teuchos::Tuple<int,3> ncS;

    for( int i=0; i<3; ++i ) {
      ncC[i] = ncbC_[i][dim_nc-1];
      ncS[i] = ncbG_[i][dim_nc-1];
    }

    // zentral
    for( int i=0; i<3; ++i ) {
      bu_[i] = ncS[i]/2;
      bl_[i] = -bu_[i];
    }
    if( 4==dim ) bl_[3] = -1;
    if( 4==dim ) bu_[3] =  0;

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
  FieldSpace(
//      Ordinal dimension,
      TO bl,
      TO bu ):
        bl_(bl),
        bu_(bu),
        ls_(Teuchos::tuple(-1,-1,-1))
  {};

//  /// \todo necessary?
//  FieldSpace(
//      const FieldSpace& fs ):
//        bl_(fs.bl_),
//        bu_(fs.bu_),
//        ls_(Teuchos::tuple(-1,-1,-1))
//  {};

public:

  /// prints to \c std::cout, only for debuging purpose
  void print( std::ostream& out=std::cout ) const {
    out << "\t---FieldSpace: ---\n";
    out << "comput dim: " << dim << "\n";
    out << "ncbC: " << ncbC_ << "\n";
    out << "ncbD: " << ncbD_ << "\n";
    out << "ncbG: " << ncbG_ << "\n";
    out << "bl: " << bl_ << "\n";
    out << "bu: " << bu_ << "\n";
    out << "dl: " << dl_ << "\n";
    out << "du: " << du_ << "\n";
    out << "gl: " << gl_ << "\n";
    out << "gu: " << gu_ << "\n";
    out << "nl: " << nl_ << "\n";
    out << "nu: " << nu_ << "\n";
    out << "ls: " << ls_ << "\n";
  }

        Ordinal  getDimNcbC( int i ) const { return( ncbC_[i].size() ); }
        Ordinal  getDimNcbD( int i ) const { return( ncbD_[i].size() ); }
        Ordinal  getDimNcbG( int i ) const { return( ncbG_[i].size() ); }

  const Ordinal* getNcbC( int i ) const { return( ncbC_[i].getRawPtr() ); }
  const Ordinal* getNcbD( int i ) const { return( ncbD_[i].getRawPtr() ); }
  const Ordinal* getNcbG( int i ) const { return( ncbG_[i].getRawPtr() ); }

  const Ordinal* getBL()        const { return( bl_.getRawPtr() ); }
  const Ordinal& getBL( int i ) const { return( bl_[i] ); }

  const Ordinal* getBU()        const { return( bu_.getRawPtr() ); }
  const Ordinal& getBU( int i ) const { return( bu_[i] ); }

  const Ordinal* getDL()        const { return( dl_.getRawPtr() ); }
  const Ordinal& getDL( int i ) const { return( dl_[i] ); }

  const Ordinal* getDU()        const { return( du_.getRawPtr() ); }
  const Ordinal& getDU( int i ) const { return( du_[i] ); }

  const Ordinal* getGL()        const { return( gl_.getRawPtr() ); }
  const Ordinal& getGL( int i ) const { return( gl_[i] ); }

  const Ordinal* getGU()        const { return( gu_.getRawPtr() ); }
  const Ordinal& getGU( int i ) const { return( gu_[i] ); }

  const Ordinal* getNL()        const { return( nl_.getRawPtr() ); }
  const Ordinal& getNL( int i ) const { return( nl_[i] ); }

  const Ordinal* getNU()        const { return( nu_.getRawPtr() ); }
  const Ordinal& getNU( int i ) const { return( nu_[i] ); }

  const Ordinal* getLS()        const { return( ls_.getRawPtr() ); }
  const Ordinal& getLS( int i ) const { return( ls_[i] ); }

}; // end of class FieldSpace



/// \brief function that creates \c Pimpact:FieldSpace
/// by getting values from \c IMPACT
/// should be changed by getting vallues from \c ProcGridSize and \c GridSize
/// \relates FieldSpace
template< class O=int, int d=3 >
const Teuchos::RCP<const FieldSpace<O,d> > createFieldSpace(){

  return(
      Teuchos::RCP<const FieldSpace<O,d> > (
          new FieldSpace<O,d>() ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_FIELDSPACE_HPP
