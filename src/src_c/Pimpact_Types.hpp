#pragma once
#ifndef PIMPACT_TYPES_HPP
#define PIMPACT_TYPES_HPP



namespace Pimpact {

	enum EFieldType  { U, V, W, S };

	/// Copy Type
	enum ECopyType {
	  /// Deep Copy, means that everything is copied including boundaries
	  DeepCopy,
	  /// Schallow Copy, up to now new field is initialized to zero.
	  ShallowCopy };


	enum EFlowProfile {
	  ZeroProf=0,
		Poiseuille2D_inX=1, Poiseuille2D_inY=2,
		Pulsatile2D_inXC=3, Pulsatile2D_inXS=5,
		Pulsatile2D_inYC=4, Pulsatile2D_inYS=6,
		Streaming2D=7,
		Circle2D=8,
		RankineVortex2D=9,
		GaussianForcing1D=10,
		BoundaryFilter1D=11,
		GaussianForcing2D=12,
		BoundaryFilter2D=13
	};

	enum EFlowType { Zero2DFlow = 0,
			Poiseuille_inX=1, Poiseuille_inY=2,
			Pulsatile_inX=3, Pulsatile_inY=4,
			Streaming2DFlow=5
	};


	enum EBCType {
	  SymmetryBC = -2,
	  PeriodicBC = -1,
	  NeighborBC = 0,
	  DirichletBC = 1,
	  NeumannBC = 2,
	  RobinBC = 3
	};

	enum EDomainType {
	  AllDirichlet = 0,
	  Dirichelt2DChannel = 1,
	  Periodic2DChannel = 2,
	  AllNeumann2D = 3,
	  AllPeriodic = 4,
	  Neumann1Periodic2 = 5
	};


	struct ModeOp {};
	struct NonModeOp {};
}
#endif // PIMPACT_TYPES_HPP
