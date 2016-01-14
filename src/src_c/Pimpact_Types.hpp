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
enum EField { U=0, V=1, W=2, S=4 };



/// \brief kind of boundary conditions
/// \relates BoundaryConditionsGlobal
/// \relates BoundaryConditionsLocal
enum EBCType {
  SymmetryBC = -2,
  PeriodicBC = -1,
  NeighborBC = 0,
  DirichletBC = 1,
  NeumannBC = 2,
  RobinBC = 3
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
	Poiseuille2D_inX = 1,
  Poiseuille2D_inY = 2,
  Poiseuille2D_inZ = 6,
  Grad2D_inX = 3,
  Grad2D_inY = 4,
  Grad2D_inZ = 5,
  FPoint = 7
};



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
  OscilatingDisc2DVel=9
};


/// \brief kind of force
enum EForceType {
  Dipol        = 1,
  Disc         = 2,
  RotatingDisc = 3,
  PseudoOscilatingDisc = 4
};

Teuchos::RCP<std::ostream> createOstream( const std::string& fname, int rank);

}


#endif // PIMPACT_TYPES_HPP
