#pragma once
#ifndef PIMPACT_TYPES_HPP
#define PIMPACT_TYPES_HPP


#include "BelosTypes.hpp"




namespace Pimpact {



/// \brief Copy Type
enum ECopyType {
  /// Deep Copy, means that everything is copied including boundaries
  DeepCopy,
  /// Schallow Copy, up to now new field is initialized to zero.
  ShallowCopy
};


enum ECoord { X=0, Y=1, Z=2, T=4 };

std::string toString( ECoord type );

enum EField { U=0, V=1, W=2, S=4 };

std::string toString( EField type );


/// \brief kind of boundary conditions
/// \relates BoundaryConditionsGlobal
/// \relates BoundaryConditionsLocal
enum EBCType {
  SymmetryBC = -2, ///< -2
  PeriodicBC = -1, ///< -1
  NeighborBC = 0,  ///<  0 
  DirichletBC = 1, ///<  1 
  NeumannBC = 2,   ///<  2 
  RobinBC = 3      ///<  3 
};



/// \brief type of Domain, e.g. periodic channel, box, ...
/// 
/// \relates BoundaryConditionsGlobal
enum EDomainType {
  AllDirichlet = 0,
  Dirichelt2DChannel = 1,
  Periodic2DChannel = 2,
  AllNeumann2D = 3,
  AllPeriodic = 4,
  Neumann1Periodic2 = 5
};



/// \brief Scalar Field profile
/// \relates ScalarField::initField
enum EScalarField {
	ConstField = 0,
  Grad2D_inX = 1,
  Grad2D_inY = 2,
  Grad2D_inZ = 3,
	Poiseuille2D_inX = 4,
  Poiseuille2D_inY = 5,
  Poiseuille2D_inZ = 6,
  FPoint = 7
};


int getDir1( const int& dir );
int getDir2( const int& dir );

std::string toString( Pimpact::EScalarField type );

/// \brief used in mains
/// \relates initTimeVectorField
enum EFlowType {
  Zero2DFlow = 0, 
  Const2DFlow = 10, 
  Poiseuille_inX=1,
  Poiseuille_inY=2,
  Pulsatile_inX=3,
  Pulsatile_inY=4,
  Streaming2DFlow=5,
  Streaming2DFlow2=6,
  Streaming2DFlow3=7,
  OscilatingDisc2D=8,
  OscilatingDisc2DVel=9,
  ConstVel_inX=11
};


/// \brief kind of force
enum EForceType {
  Dipol        = 1,
  Disc         = 2,
  RotatingDisc = 3,
  PseudoOscilatingDisc = 4
};

/// \brief creates file and stream to it
///
/// \param fname name of file
/// \param rank if rank==0 then file is created otherwise ostream is a blackhole, 'default'=0
///
/// \return  pointer to a std::ofstream
Teuchos::RCP<std::ostream> createOstream( const std::string& fname, int rank=0 );

}


#endif // PIMPACT_TYPES_HPP
