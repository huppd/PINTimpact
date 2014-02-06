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
	enum EBCType { SymmetryBC=-2, PeriodicBC=-1, NeighborBC=0, DirichletBC=1, NeumannBC=2, RobinBC=3 };

	enum EFlowProfile { ZeroProf=0,
			Poiseuille2D_inX=1, Poiseuille2D_inY=2,
			Pulsatile2D_inXC=3, Pulsatile2D_inXS=5,
			Pulsatile2D_inYC=4, Pulsatile2D_inYS=6,
			Streaming2D=7
	};

	enum EFlowType { ZeroFLow = 0,
			Poiseuille_inX=1, Poiseuille_inY=2,
			Pulsatile_inX=3, Pulsatile_inY=4,
			Streaming2DFlow=5
	};

	enum EDomainType { AllDirichlet = 0, Dirichelt2DChannel=1, Periodic2DChannel = 2 };

	struct ModeOp {};
	struct NonModeOp {};
}
#endif // PIMPACT_TYPES_HPP
