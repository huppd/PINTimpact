#include "Pimpact_Utils.hpp"



std::string Pimpact::toString( Pimpact::ECoord type ) {
  switch( type ) {
  case Pimpact::ECoord::X :
    return( "X" );
  case Pimpact::ECoord::Y :
    return( "Y" );
  case Pimpact::ECoord::Z :
    return( "Z" );
  case Pimpact::ECoord::T :
    return( "T" );
  }
  //return( "" ); // prevent compiler warning
}


Pimpact::ECoord& operator++( Pimpact::ECoord& c ) {
  switch( c ) {
  case Pimpact::ECoord::X   :
    return( c = Pimpact::ECoord::Y );
  case Pimpact::ECoord::Y   :
    return( c = Pimpact::ECoord::Z );
  case Pimpact::ECoord::Z   :
    return( c = Pimpact::ECoord::end );
  case Pimpact::ECoord::end :
    return( c = Pimpact::ECoord::end );
  case Pimpact::ECoord::T   :
    return( c = Pimpact::ECoord::end );
  }
}

bool operator<( const Pimpact::ECoord& c, const int& i ) {
  return( static_cast<int>(c)<i );
}


std::string Pimpact::toString( Pimpact::F type ) {
  switch( type ) {
  case Pimpact::F::U :
    return( "U" );
  case Pimpact::F::V :
    return( "V" );
  case Pimpact::F::W :
    return( "W" );
  case Pimpact::F::S :
    return( "S" );
  }
  return( "" ); // prevent compiler warning
}


Pimpact::F& operator++( Pimpact::F& c ) {
  switch( c ) {
  case Pimpact::F::U   :
    return( c = Pimpact::F::V );
  case Pimpact::F::V   :
    return( c = Pimpact::F::W );
  case Pimpact::F::W   :
    return( c = Pimpact::F::end );
  case Pimpact::F::end :
    return( c = Pimpact::F::end );
  case Pimpact::F::S   :
    return( c = Pimpact::F::end );
  }
  return( c );
}

bool operator<( const Pimpact::F& c, const int& i ) {
  return( static_cast<int>(c)<i );
}

bool operator==( const Pimpact::F& f, const int& c ) {
  if( Pimpact::F::U==f &&  Pimpact::ECoord::X==c ) return( true );
  else if( Pimpact::F::V==f &&  Pimpact::ECoord::Y==c ) return( true );
  else if( Pimpact::F::W==f &&  Pimpact::ECoord::Z==c ) return( true );
  else return( false );
}

bool operator!=( const Pimpact::F& f, const int& c ) {
  if( Pimpact::F::U==f &&  Pimpact::ECoord::X==c ) return( false );
  else if( Pimpact::F::V==f &&  Pimpact::ECoord::Y==c ) return( false );
  else if( Pimpact::F::W==f &&  Pimpact::ECoord::Z==c ) return( false );
  else return( true );
}


bool operator==( const int& c, const Pimpact::F& f ) {
  if( Pimpact::F::U==f &&  Pimpact::ECoord::X==c ) return( true );
  else if( Pimpact::F::V==f &&  Pimpact::ECoord::Y==c ) return( true );
  else if( Pimpact::F::W==f &&  Pimpact::ECoord::Z==c ) return( true );
  else return( false );
}

bool operator!=( const int& c, const Pimpact::F& f ) {
  if( Pimpact::F::U==f &&  Pimpact::ECoord::X==c ) return( false );
  else if( Pimpact::F::V==f &&  Pimpact::ECoord::Y==c ) return( false );
  else if( Pimpact::F::W==f &&  Pimpact::ECoord::Z==c ) return( false );
  else return( true );
}

std::string Pimpact::toString( Pimpact::EScalarField type ) {
  switch( type ) {
  case Pimpact::EScalarField::ConstField :
    return( "constant" );
  case Pimpact::EScalarField::Grad2D_inX :
    return( "grad in x" );
  case Pimpact::EScalarField::Grad2D_inY :
    return( "grad in y" );
  case Pimpact::EScalarField::Grad2D_inZ :
    return( "grad in z" );
  case Pimpact::EScalarField::Poiseuille2D_inX :
    return( "poiseuille in x" );
  case Pimpact::EScalarField::Poiseuille2D_inY :
    return( "poiseuille in y" );
  case Pimpact::EScalarField::Poiseuille2D_inZ :
    return( "poiseuille in z" );
  case Pimpact::EScalarField::FPoint :
    return( "poiseuille in z" );
  }
  return( "" ); // prevent compiler warning
}


int Pimpact::getDir1( const int& dir ) {

  int dir1 = ( dir + 1 )%3;
  int dir2 = ( dir + 2 )%3;
  if( dir2<dir1 )
    return( dir2 );
  return( dir1 );
}

int Pimpact::getDir2( const int& dir ) {

  int dir1 = ( dir + 1 )%3;
  int dir2 = ( dir + 2 )%3;
  if( dir2>dir1 )
    return( dir2 );
  return( dir1 );
}


Teuchos::RCP<std::ostream>
Pimpact::createOstream( const std::string& fname, int rank ) {

  if( 0==rank )
    return( Teuchos::rcp( new std::ofstream( fname ) ) );
  else
    return( Teuchos::rcp( new Teuchos::oblackholestream ) );

}

void Pimpact::setBoundaryConditions( const
                                     Teuchos::RCP<Teuchos::ParameterList>& pl , int dtype ) {

  switch( static_cast<Pimpact::EDomainType>(dtype) ) {
  case Pimpact::AllDirichlet:
    pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::BC::Dirichlet );
    break;
  case Pimpact::AllPeriodic:
    pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::BC::Periodic );
    pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::BC::Periodic );
    pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::BC::Periodic );
    pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::BC::Periodic );
    pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::BC::Periodic );
    pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::BC::Periodic );
    break;
  case Pimpact::AllNeumann:
    pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::BC::Neumann );
    pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::BC::Neumann );
    pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::BC::Neumann );
    pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::BC::Neumann );
    pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::BC::Neumann );
    pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::BC::Neumann );
    break;
  case Pimpact::AllSymmetric:
    pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::BC::Symmetry );
    pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::BC::Symmetry );
    pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::BC::Symmetry );
    pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::BC::Symmetry );
    pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::BC::Symmetry );
    pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::BC::Symmetry );
    break;
  case Pimpact::Dirichelt2DChannel:
    pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::BC::Periodic );
    pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::BC::Periodic );
    break;
  case Pimpact::Periodic2DChannel:
    pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::BC::Periodic );
    pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::BC::Periodic );
    pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::BC::Periodic );
    pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::BC::Periodic );
    break;
  case Pimpact::Open2DChannel:
    pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::BC::Neumann   );
    pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::BC::Dirichlet );
    pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::BC::Periodic  );
    pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::BC::Periodic  );
    break;
  default:
    std::cout << "!!!Warning: unkown EDomainType:\t" <<dtype<<"\t!!!\n";
  }
}
