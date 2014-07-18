#pragma once
#ifndef PIMPACT_INDEXSPACE_HPP
#define PIMPACT_INDEXSPACE_HPP

#include "mpi.h"

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include <iostream>

#include "Pimpact_Types.hpp"



namespace Pimpact {


extern "C" {
void SVS_get_sInd(int&,int&,int&);
void SVS_get_eInd(int&,int&,int&);

void VS_get_sIndU(int&,int&,int&);
void VS_get_eIndU(int&,int&,int&);

void VS_get_sIndV(int&,int&,int&);
void VS_get_eIndV(int&,int&,int&);

void VS_get_sIndW(int&,int&,int&);
void VS_get_eIndW(int&,int&,int&);

void VS_get_sIndUB(int&,int&,int&);
void VS_get_eIndUB(int&,int&,int&);

void VS_get_sIndVB(int&,int&,int&);
void VS_get_eIndVB(int&,int&,int&);

void VS_get_sIndWB(int&,int&,int&);
void VS_get_eIndWB(int&,int&,int&);
}


/// \brief public class, that stores neccessary information for indexing \c ScalarField and \c VectorField in Fortran
/// \todo check constant make variables protected SF friend class
template< class Ordinal = int >
class IndexSpace {

public:

  typedef Teuchos::Tuple<Ordinal,3> TO3;

  /// \brief constructor
  /// \param fieldType says which kind of field is taken
  /// \param sInd start index of gridpoints
  /// \param eInd last index of gridpoints
  IndexSpace( EFieldType fieldType,TO3 sInd, TO3 eInd ):
    fieldType_(fieldType),sInd_(sInd),eInd_(eInd) {};

  void print( std::ostream& out=std::cout ) const {
    out << "\t---IndexSpace: ---\n";
    out << "fieldType: " << fieldType_ << "\n";
    out << "sInd: " << sInd_ << "\n";
    out << "eInd: " << eInd_ << "\n";
  }

  EFieldType fieldType_;

  TO3 sInd_;
  TO3 eInd_;

}; // end of class IndexSpace



/// \brief function that creates ScaparIndexSpace
/// by getting values from \c IMPACT
/// \relates IndexSpace
template<class Ordinal=int>
Teuchos::RCP<const IndexSpace<Ordinal> >
createScalarIndexSpace(){

  typedef typename IndexSpace<Ordinal>::TO3 TO3;

  TO3 sInd;
  SVS_get_sInd( sInd[0], sInd[1], sInd[2] );

  TO3 eInd;
  SVS_get_eInd( eInd[0], eInd[1], eInd[2] );

  return( Teuchos::rcp(
      new IndexSpace<Ordinal>( EFieldType::S, sInd, eInd ) ) );

}


/// \brief function that creates ScaparIndexSpace
///
/// by getting values from impact
/// \todo make construction clean without copying
/// \relates IndexSpace
template<class Ordinal=int>
Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<Ordinal> > > createInnerFieldIndexSpaces(){

  typedef typename IndexSpace<Ordinal>::TO3 TO3;

  Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > > fIS(3);

  TO3 sInd;
  TO3 eInd;

  VS_get_sIndU( sInd[0], sInd[1], sInd[2] );
  VS_get_eIndU( eInd[0], eInd[1], eInd[2] );
  fIS[0] =	Teuchos::rcp( new IndexSpace<Ordinal>( EFieldType::U, sInd, eInd ) );

  VS_get_sIndV( sInd[0], sInd[1], sInd[2] );
  VS_get_eIndV( eInd[0], eInd[1], eInd[2] );
  fIS[1] =	Teuchos::rcp( new IndexSpace<Ordinal>( EFieldType::V, sInd, eInd ) );

  VS_get_sIndW( sInd[0], sInd[1], sInd[2] );
  VS_get_eIndW( eInd[0], eInd[1], eInd[2] );
  fIS[2] =	Teuchos::rcp( new IndexSpace<Ordinal>( EFieldType::W, sInd, eInd ) );

  return( fIS );
}


/// \brief function that creates ScaparIndexSpace
/// by getting values from impact
/// \todo make construction clean without copying
/// \relates IndexSpace
template<class Ordinal>
Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<Ordinal> > > createFullFieldIndexSpaces(){

  typedef typename IndexSpace<Ordinal>::TO3 TO3;

  Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > > fIS(3);

  TO3 sInd;
  TO3 eInd;

  VS_get_sIndUB( sInd[0], sInd[1], sInd[2] );
  VS_get_eIndUB( eInd[0], eInd[1], eInd[2] );
  fIS[0] =	Teuchos::rcp( new IndexSpace<Ordinal>( EFieldType::U, sInd, eInd ) );

  VS_get_sIndVB( sInd[0], sInd[1], sInd[2] );
  VS_get_eIndVB( eInd[0], eInd[1], eInd[2] );
  fIS[1] =	Teuchos::rcp( new IndexSpace<Ordinal>( EFieldType::V, sInd, eInd ) );

  VS_get_sIndWB( sInd[0], sInd[1], sInd[2] );
  VS_get_eIndWB( eInd[0], eInd[1], eInd[2] );
  fIS[2] =	Teuchos::rcp( new IndexSpace<Ordinal>( EFieldType::W, sInd, eInd ) );

  return( fIS );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INDEXSPACE_HPP
