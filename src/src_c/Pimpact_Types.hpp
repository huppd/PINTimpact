#pragma once
#ifndef PIMPACT_TYPES_HPP
#define PIMPACT_TYPES_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "BelosTypes.hpp"




namespace Pimpact {



/// \brief Copy Type
enum class ECopy : bool {
  /// Deep Copy, means that everything is copied including boundaries
  Deep,
  /// Schallow Copy, up to now new field is initialized to zero.
  Shallow
};


enum ECoord {
	X=0, ///< 0
	Y=1, ///< 1
	Z=2, ///< 2
	T=4  ///< 4
};

std::string toString( ECoord type );

enum EField {
	U=0, ///< 0
	V=1, ///< 0
	W=2, ///< 2
	S=4  ///< 4
};

std::string toString( EField type );


/// \brief kind of boundary conditions
/// \relates BoundaryConditionsGlobal
/// \relates BoundaryConditionsLocal
enum EBCType {
  SymmetryBC  = -2, ///< -2
  PeriodicBC  = -1, ///< -1
  NeighborBC  =  0, ///<  0 
  DirichletBC =  1, ///<  1 
  NeumannBC   =  2, ///<  2 
  RobinBC     =  3  ///<  3 
};



/// \brief type of Domain, e.g. periodic channel, box, ...
/// 
/// \relates BoundaryConditionsGlobal
enum EDomainType {
  AllDirichlet       = 0,	///< 0
  AllPeriodic        = 1,	///< 1
  AllNeumann         = 2, ///< 2
  AllSymmetric       = 3, ///< 3
  Dirichelt2DChannel = 4, ///< 4
  Periodic2DChannel  = 5, ///< 5
  Open2DChannel      = 6  ///< 6
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


enum class With : bool {
	BC   = true,	///< [true] including grid points on the boundary
	noBC = false 	///< [false] only inner points (DOF)
};

int getDir1( const int& dir );
int getDir2( const int& dir );

std::string toString( Pimpact::EScalarField type );


/// \brief creates file and stream to it
///
/// \param fname name of file
/// \param rank if rank==0 then file is created otherwise ostream is a blackhole, 'default'=0
///
/// \return  pointer to a std::ofstream
Teuchos::RCP<std::ostream> createOstream( const std::string& fname, int rank=0 );


void setBoundaryConditions( const Teuchos::RCP<Teuchos::ParameterList>& pl ,
		int dtype );

}




#endif // PIMPACT_TYPES_HPP
