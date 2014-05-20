#pragma once
#ifndef PIMPACT_GRID_HPP
#define PIMPACT_GRID_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



namespace Pimpact{


template<class Scalar, class Ordinal>
class Grid {

public:

  Teuchos::RCP<DomainSize<Scalar> > domainSize_;

  Teuchos::RCP<BoundaryConditionsGlobal> bcGlo_;
  Teuchos::RCP<BoundaryConditionsLocal> bcLoc_;

//  Teuchos::RCP<ProcGrid> procGrid_;

};


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRID_HPP
