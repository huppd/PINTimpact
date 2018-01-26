#pragma once
#ifndef PIMPACT_UTILS_HPP
#define PIMPACT_UTILS_HPP


#include <fstream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "BelosTypes.hpp"




namespace Pimpact {


/// \brief Norm Type
enum class ENorm : int {
  One, /// < vector one norm
  Two, /// < vector two norm
  Inf, /// < vector inf norm
  L2   /// < function 2 norm
};


/// \brief Copy Type
enum class ECopy : bool {
  Deep,		///< Deep Copy, means that everything is copied including boundaries
  Shallow	///< Shallow Copy, new field is initialized with zero.
};


/// \todo change T->3, consistently
/// \todo make enum class 
/// \todo make it iterable( <, ++ )
enum ECoord {
  X = 0, ///< 0
  Y = 1, ///< 1
  Z = 2, ///< 2
  T = 4  ///< 4
};

enum class Owning : bool { Y=true, N=false };

std::string toString( ECoord type ) {
  switch( type ) {
    case ECoord::X : return "X" ;
    case ECoord::Y : return "Y" ;
    case ECoord::Z : return "Z" ;
    case ECoord::T : return "T" ;
  }
  return ""; // prevent compiler warning
}

std::ostream& operator<<( std::ostream& out, ECoord c ) {
  return out << toString(c);
}

// Special behavior for ++
ECoord& operator++( ECoord& c );

bool operator<( const ECoord c, const int& i );


enum class F : int {
  U=0, ///< 0
  V=1, ///< 1
  W=2, ///< 2
  S=4  ///< 4
};

std::string toString( F type ) {
  switch( type ) {
    case F::U : return "U";
    case F::V : return "V";
    case F::W : return "W";
    case F::S : return "S";
  }
  return ""; // prevent compiler warning
}

std::ostream& operator<<( std::ostream& out, F f ) {
  return out << toString(f);
}

F& operator++( F& c ) {
  switch( c ) {
    case F::U : return c = F::V;
    case F::V : return c = F::W;
    case F::W : return c = F::S;
    case F::S : return c = F::S;
  }
  return c;
}

bool operator<( const F& c, int i ) {
  return static_cast<int>(c)<i;
}

bool operator==( const F& f, const int& c ) {
  if( F::U==f &&  ECoord::X==c ) return true;
  else if( F::V==f &&  ECoord::Y==c ) return true;
  else if( F::W==f &&  ECoord::Z==c ) return true;
  else return false;
}

bool operator!=( const F& f, const int& c ) {
  if( F::U==f &&  ECoord::X==c ) return false;
  else if( F::V==f &&  ECoord::Y==c ) return false;
  else if( F::W==f &&  ECoord::Z==c ) return false;
  else return true;
}


bool operator==( const int& c, const F& f ) {
  if( F::U==f &&  ECoord::X==c ) return true;
  else if( F::V==f &&  ECoord::Y==c ) return true;
  else if( F::W==f &&  ECoord::Z==c ) return true;
  else return false;
}

bool operator!=( const int& c, const F& f ) {
  if( F::U==f &&  ECoord::X==c ) return false;
  else if( F::V==f &&  ECoord::Y==c ) return false;
  else if( F::W==f &&  ECoord::Z==c ) return false;
  else return true;
}

/// \brief kind of boundary conditions
/// \relates BoundaryConditionsGlobal
/// \relates BoundaryConditionsLocal
/// \todo make enumclass impl ()
enum BC {
  Symmetry  = -2, ///< -2
  Periodic  = -1, ///< -1
  Neighbor  =  0, ///<  0
  Dirichlet =  1, ///<  1
  Neumann   =  2, ///<  2
  Robin     =  3  ///<  3
};


bool operator<( const int& i, const BC& c ) {
  return i<static_cast<int>(c);
}
//bool operator<( const BC& left, const BC& right ) {
//return static_cast<int>(left)<static_cast<int>(right);
//}


/// \brief type of Domain, e.g. periodic channel, box, ...
///
/// \relates BoundaryConditionsGlobal
enum EDomainType {
  AllDirichlet       = 0,	///< 0
  AllPeriodic        = 1,	///< 1
  AllNeumann         = 2, ///< 2
  AllSymmetric       = 3, ///< 3
  Dirichelt2DChannel = 4, ///< 4 ( BC::Dirichlet,             BC::Dirichlet, BC::Periodic )
  Periodic2DChannel  = 5, ///< 5 ( BC::Periodic,              BC::Dirichelt, BC::Periodic )
  Open2DChannel      = 6  ///< 6 ( BC::Dirichlet/BC::Neumann, BC::Dirichlet, BC::Periodic )
};



/// \brief Scalar Field profile
/// \relates ScalarField::initField
enum EScalarField {
  ConstField = 0,	      ///< 0
  Grad2D_inX = 1,       ///< 1
  Grad2D_inY = 2,       ///< 2
  Grad2D_inZ = 3,       ///< 3
  Poiseuille2D_inX = 4, ///< 4
  Poiseuille2D_inY = 5, ///< 5
  Poiseuille2D_inZ = 6, ///< 6
  FPoint = 7
};



/// \brief used in mains
/// \relates initTimeVectorField
enum EFlowType {
  Zero2DFlow          =  0,	///<  0
  Const2DFlow         = 10, ///< 10
  Poiseuille_inX      =  1, ///<  1
  Poiseuille_inY      =  2, ///<  2
  Pulsatile_inX       =  3, ///<  3
  Pulsatile_inY       =  4, ///<  4
  Streaming2DFlow     =  5, ///<  5
  Streaming2DFlow2    =  6, ///<  6
  Streaming2DFlow3    =  7, ///<  7
  OscilatingDisc2D    =  8, ///<  8
  OscilatingDisc2DVel =  9, ///<  9
  ConstVel_inX        = 11  ///< 11
};


/// \brief kind of force
enum EForceType {
  Dipol                = 1, ///< 1
  Disc                 = 2, ///< 2
  RotatingDisc         = 3, ///< 3
  PseudoOscilatingDisc = 4  ///< 4
};


enum class B : bool {
  Y = true,	///< [true] including grid points on the boundary
  N = false	///< [false] only inner points (DOF)
};


enum class Add : bool {
  Y = true,	///< [true] adds applied operator to DomainField
  N = false	///< [false] sets DomainField to applied Operator
};


int getDir1( const int& dir );
int getDir2( const int& dir );

std::string toString( EScalarField type ) {
  switch( type ) {
    case EScalarField::ConstField : return "constant";
    case EScalarField::Grad2D_inX : return "grad in x";
    case EScalarField::Grad2D_inY : return "grad in y";
    case EScalarField::Grad2D_inZ : return "grad in z";
    case EScalarField::Poiseuille2D_inX : return "poiseuille in x";
    case EScalarField::Poiseuille2D_inY : return "poiseuille in y";
    case EScalarField::Poiseuille2D_inZ : return "poiseuille in z";
    case EScalarField::FPoint : return "poiseuille in z";
  }
  return ""; // prevent compiler warning
}

std::ostream& operator<<( std::ostream& out, const EScalarField& c ) {
  return out << toString(c);
}


/// \brief creates file and stream to it
///
/// \param fname name of file
/// \param rank if rank==0 then file is created otherwise ostream is a blackhole, 'default'=0
///
/// \return  pointer to a std::ofstream
Teuchos::RCP<std::ostream> createOstream( const std::string& fname, int rank=0, int restart=-1 ) {

  if( 0==rank ) {
    if( -1==restart )
      return Teuchos::rcp( new std::ofstream(fname) );
    else
      return Teuchos::rcp( new std::ofstream(fname, std::fstream::app) );
  }
  else
    return Teuchos::rcp( new Teuchos::oblackholestream );
}


void setBoundaryConditions( const Teuchos::RCP<Teuchos::ParameterList>& pl , int dtype ) {

  switch( static_cast<EDomainType>(dtype) ) {
  case AllDirichlet:
    pl->sublist("boundary conditions").set<int>( "lower X", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "upper X", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "lower Y", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "upper Y", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "lower Z", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "upper Z", static_cast<int>( BC::Dirichlet ) );
    break;
  case AllPeriodic:
    pl->sublist("boundary conditions").set<int>( "lower X", static_cast<int>( BC::Periodic ) );
    pl->sublist("boundary conditions").set<int>( "upper X", static_cast<int>( BC::Periodic ) );
    pl->sublist("boundary conditions").set<int>( "lower Y", static_cast<int>( BC::Periodic ) );
    pl->sublist("boundary conditions").set<int>( "upper Y", static_cast<int>( BC::Periodic ) );
    pl->sublist("boundary conditions").set<int>( "lower Z", static_cast<int>( BC::Periodic ) );
    pl->sublist("boundary conditions").set<int>( "upper Z", static_cast<int>( BC::Periodic ) );
    break;
  case AllNeumann:
    pl->sublist("boundary conditions").set<int>( "lower X", static_cast<int>( BC::Neumann ) );
    pl->sublist("boundary conditions").set<int>( "upper X", static_cast<int>( BC::Neumann ) );
    pl->sublist("boundary conditions").set<int>( "lower Y", static_cast<int>( BC::Neumann ) );
    pl->sublist("boundary conditions").set<int>( "upper Y", static_cast<int>( BC::Neumann ) );
    pl->sublist("boundary conditions").set<int>( "lower Z", static_cast<int>( BC::Neumann ) );
    pl->sublist("boundary conditions").set<int>( "upper Z", static_cast<int>( BC::Neumann ) );
    break;
  case AllSymmetric:
    pl->sublist("boundary conditions").set<int>( "lower X", static_cast<int>( BC::Symmetry ) );
    pl->sublist("boundary conditions").set<int>( "upper X", static_cast<int>( BC::Symmetry ) );
    pl->sublist("boundary conditions").set<int>( "lower Y", static_cast<int>( BC::Symmetry ) );
    pl->sublist("boundary conditions").set<int>( "upper Y", static_cast<int>( BC::Symmetry ) );
    pl->sublist("boundary conditions").set<int>( "lower Z", static_cast<int>( BC::Symmetry ) );
    pl->sublist("boundary conditions").set<int>( "upper Z", static_cast<int>( BC::Symmetry ) );
    break;
  case Dirichelt2DChannel:
    pl->sublist("boundary conditions").set<int>( "lower X", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "upper X", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "lower Y", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "upper Y", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "lower Z", static_cast<int>( BC::Periodic  ) );
    pl->sublist("boundary conditions").set<int>( "upper Z", static_cast<int>( BC::Periodic  ) );
    break;
  case Periodic2DChannel:
    pl->sublist("boundary conditions").set<int>( "lower X", static_cast<int>( BC::Periodic  ) );
    pl->sublist("boundary conditions").set<int>( "upper X", static_cast<int>( BC::Periodic  ) );
    pl->sublist("boundary conditions").set<int>( "lower Y", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "upper Y", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "lower Z", static_cast<int>( BC::Periodic  ) );
    pl->sublist("boundary conditions").set<int>( "upper Z", static_cast<int>( BC::Periodic  ) );
    break;
  case Open2DChannel:
    pl->sublist("boundary conditions").set<int>( "lower X", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "upper X", static_cast<int>( BC::Neumann   ) );
    pl->sublist("boundary conditions").set<int>( "lower Y", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "upper Y", static_cast<int>( BC::Dirichlet ) );
    pl->sublist("boundary conditions").set<int>( "lower Z", static_cast<int>( BC::Periodic  ) );
    pl->sublist("boundary conditions").set<int>( "upper Z", static_cast<int>( BC::Periodic  ) );
    break;
  default:
    std::cout << "!!!Warning: unkown EDomainType:\t" <<dtype<<"\t!!!\n";
  }
}

} // end of namespace Pipmact




#endif // PIMPACT_UTILS_HPP
