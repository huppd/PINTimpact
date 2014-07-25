#pragma once
#ifndef PIMPACT_INDEXSPACE_HPP
#define PIMPACT_INDEXSPACE_HPP

#include "mpi.h"

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include <iostream>

#include "Pimpact_Types.hpp"

#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"



namespace Pimpact {


extern "C" {
void SVS_get_sInd(int&,int&,int&);
void SVS_get_eInd(int&,int&,int&);

void SVS_set_sInd(const int* const);
void SVS_set_eInd(const int* const);

void VS_get_sIndU(int&,int&,int&);
void VS_get_eIndU(int&,int&,int&);

void VS_set_sIndU(const int* const);
void VS_set_eIndU(const int* const);

void VS_set_sIndUB(const int* const);
void VS_set_eIndUB(const int* const);

void VS_get_sIndV(int&,int&,int&);
void VS_get_eIndV(int&,int&,int&);

void VS_set_sIndV(const int* const);
void VS_set_eIndV(const int* const);

void VS_set_sIndVB(const int* const);
void VS_set_eIndVB(const int* const);

void VS_get_sIndW(int&,int&,int&);
void VS_get_eIndW(int&,int&,int&);

void VS_set_sIndW(const int* const);
void VS_set_eIndW(const int* const);

void VS_set_sIndWB(const int* const);
void VS_set_eIndWB(const int* const);

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
template<class O=int, int d=3>
Teuchos::RCP<const IndexSpace<O> >
createScalarIndexSpace(
    const Teuchos::RCP<const FieldSpace<O,d> >& fS,
    const Teuchos::RCP<GridSizeLocal<O,d> >& nLoc,
    const Teuchos::RCP<BoundaryConditionsLocal>& bc ){

  typedef typename IndexSpace<O>::TO3 TO3;

  TO3 sInd;
  TO3 eInd;

  for( int i=0; i<3; ++i ) {
    sInd[i] = 2            + fS->ls_[i];
    eInd[i] = nLoc->get(i) + fS->ls_[i];
    if( bc->BCL_int_[i] > 0 )
      sInd[i] = 1;
    if( bc->BCU_int_[i] > 0 )
      eInd[i] = nLoc->get(i);
  }

  SVS_set_sInd( sInd.getRawPtr() );
  SVS_set_eInd( eInd.getRawPtr() );


  return( Teuchos::rcp(
      new IndexSpace<O>( EFieldType::S, sInd, eInd ) ) );

}


/// \brief function that creates ScaparIndexSpace
///
/// by getting values from impact
/// \todo make construction clean without copying
/// \relates IndexSpace
template<class O=int, int d=3>
Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<O> > >
createInnerFieldIndexSpaces(
    const Teuchos::RCP<const FieldSpace<O,d> >& fS,
    const Teuchos::RCP<GridSizeLocal<O,d> >& nLoc,
    const Teuchos::RCP<BoundaryConditionsLocal>& bc ){


  typedef typename IndexSpace<O>::TO3 TO3;

  Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > > fIS(3);

  TO3 sIndU;
  TO3 eIndU;

  TO3 sIndV;
  TO3 eIndV;

  TO3 sIndW;
  TO3 eIndW;

  for( int i=0; i<3; ++i ) {
    sIndU[i] = 2            + fS->ls_[i];
    eIndU[i] = nLoc->get(i) + fS->ls_[i];

    sIndV[i] = 2            + fS->ls_[i];
    eIndV[i] = nLoc->get(i) + fS->ls_[i];

    sIndW[i] = 2            + fS->ls_[i];
    eIndW[i] = nLoc->get(i) + fS->ls_[i];
  }

  // Lower index in x-direction
  int i = 0;
  if( bc->BCL_int_[i] > 0 ) {
    sIndU[i] = 1;
    sIndV[i] = 2;
    sIndW[i] = 2;
  }
  // lower index in y-direction
  i = 1;
  if( bc->BCL_int_[i] > 0 ) {
    sIndU[i] = 2;
    sIndV[i] = 1;
    sIndW[i] = 2;
  }
  // lower index in z-direction
  i = 2;
  if( bc->BCL_int_[i] > 0 ) {
    sIndU[i] = 2;
    sIndV[i] = 2;
    sIndW[i] = 1;
  }

  // lower index for symmetriBC
  for( int i=0; i<3; ++i )
    if( bc->BCL_local_[i]==SymmetryBC ) {
      sIndU[i] = 1;
      sIndV[i] = 1;
      sIndW[i] = 1;
    }

  // upper index
  for( int i=0; i<3; ++i )
    if( bc->BCU_int_[i] > 0 ) {
      eIndU[i] = nLoc->get(i)-1;
      eIndV[i] = nLoc->get(i)-1;
      eIndW[i] = nLoc->get(i)-1;
    }

  // upper index in x-direction for symmetricBC
  i=0;
  if( bc->BCU_local_[i]==SymmetryBC ) {
    eIndU[i] = nLoc->get(i)-1;
    eIndV[i] = nLoc->get(i);
    eIndW[i] = nLoc->get(i);
  }

  // upper index in y-direction for symmetricBC
  i=1;
  if( bc->BCU_local_[i]==SymmetryBC ) {
    eIndU[i] = nLoc->get(i);
    eIndV[i] = nLoc->get(i)-1;
    eIndW[i] = nLoc->get(i);
  }

  // upper index in z-direction for symmetricBC
  i=2;
  if( bc->BCU_local_[i]==SymmetryBC ) {
    eIndU[i] = nLoc->get(i);
    eIndV[i] = nLoc->get(i);
    eIndW[i] = nLoc->get(i)-1;
  }


  VS_set_sIndU( sIndU.getRawPtr() );
  VS_set_eIndU( eIndU.getRawPtr() );
  fIS[0] =  Teuchos::rcp( new IndexSpace<O>( EFieldType::U, sIndU, eIndU ) );

  VS_set_sIndV( sIndV.getRawPtr() );
  VS_set_eIndV( eIndV.getRawPtr() );
  fIS[1] =  Teuchos::rcp( new IndexSpace<O>( EFieldType::V, sIndV, eIndV ) );

  VS_set_sIndW( sIndW.getRawPtr() );
  VS_set_eIndW( eIndW.getRawPtr() );
  fIS[2] =  Teuchos::rcp( new IndexSpace<O>( EFieldType::W, sIndW, eIndW ) );

  return( fIS );
}
/// \brief function that creates ScaparIndexSpace
///
/// by getting values from impact
/// \todo make construction clean without copying
/// \relates IndexSpace
template<class O=int, int d=3>
Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<O> > >
createFullFieldIndexSpaces(
    const Teuchos::RCP<const FieldSpace<O,d> >& fS,
    const Teuchos::RCP<GridSizeLocal<O,d> >& nLoc,
    const Teuchos::RCP<BoundaryConditionsLocal>& bc ){


  typedef typename IndexSpace<O>::TO3 TO3;

  Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > > fIS(3);

  TO3 sIndU;
  TO3 eIndU;

  TO3 sIndV;
  TO3 eIndV;

  TO3 sIndW;
  TO3 eIndW;

  for( int i=0; i<3; ++i ) {
    sIndU[i] = 2            + fS->ls_[i];
    eIndU[i] = nLoc->get(i) + fS->ls_[i];

    sIndV[i] = 2            + fS->ls_[i];
    eIndV[i] = nLoc->get(i) + fS->ls_[i];

    sIndW[i] = 2            + fS->ls_[i];
    eIndW[i] = nLoc->get(i) + fS->ls_[i];
  }

  // Lower index in x-direction
  int i = 0;
  if( bc->BCL_int_[i] > 0 ) {
    sIndU[i] = 0;
    sIndV[i] = 1;
    sIndW[i] = 1;
  }
  // lower index in y-direction
  i = 1;
  if( bc->BCL_int_[i] > 0 ) {
    sIndU[i] = 1;
    sIndV[i] = 0;
    sIndW[i] = 1;
  }
  // lower index in z-direction
  i = 2;
  if( bc->BCL_int_[i] > 0 ) {
    sIndU[i] = 1;
    sIndV[i] = 1;
    sIndW[i] = 0;
  }

  // lower index for symmetriBC
  for( int i=0; i<3; ++i )
    if( bc->BCL_local_[i]==SymmetryBC ) {
      sIndU[i] = 1;
      sIndV[i] = 1;
      sIndW[i] = 1;
    }

  // upper index
  for( int i=0; i<3; ++i )
    if( bc->BCU_int_[i] > 0 ) {
      eIndU[i] = nLoc->get(i);
      eIndV[i] = nLoc->get(i);
      eIndW[i] = nLoc->get(i);
    }

  // upper index in x-direction for symmetricBC
  i=0;
  if( bc->BCU_local_[i]==SymmetryBC ) {
    eIndU[i] = nLoc->get(i)-1;
    eIndV[i] = nLoc->get(i);
    eIndW[i] = nLoc->get(i);
  }

  // upper index in y-direction for symmetricBC
  i=1;
  if( bc->BCU_local_[i]==SymmetryBC ) {
    eIndU[i] = nLoc->get(i);
    eIndV[i] = nLoc->get(i)-1;
    eIndW[i] = nLoc->get(i);
  }

  // upper index in z-direction for symmetricBC
  i=2;
  if( bc->BCU_local_[i]==SymmetryBC ) {
    eIndU[i] = nLoc->get(i);
    eIndV[i] = nLoc->get(i);
    eIndW[i] = nLoc->get(i)-1;
  }


  VS_set_sIndUB( sIndU.getRawPtr() );
  VS_set_eIndUB( eIndU.getRawPtr() );
  fIS[0] =  Teuchos::rcp( new IndexSpace<O>( EFieldType::U, sIndU, eIndU ) );

  VS_set_sIndVB( sIndV.getRawPtr() );
  VS_set_eIndVB( eIndV.getRawPtr() );
  fIS[1] =  Teuchos::rcp( new IndexSpace<O>( EFieldType::V, sIndV, eIndV ) );

  VS_set_sIndWB( sIndW.getRawPtr() );
  VS_set_eIndWB( eIndW.getRawPtr() );
  fIS[2] =  Teuchos::rcp( new IndexSpace<O>( EFieldType::W, sIndW, eIndW ) );

  return( fIS );
}


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
