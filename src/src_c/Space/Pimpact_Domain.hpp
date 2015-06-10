#pragma once
#ifndef PIMPACT_DOMAIN_HPP
#define PIMPACT_DOMAIN_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_DomainSize.hpp"
#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"





namespace Pimpact{


/// \ingroup SpaceObject
template<class Scalar, int dim>
class Domain {

public:

  Domain(
      const Teuchos::RCP<const DomainSize<Scalar> >& domainSize,
      const Teuchos::RCP<const BoundaryConditionsGlobal<dim> > bcGlo,
      const Teuchos::RCP<const BoundaryConditionsLocal> bcLoc ):
        domainSize_( domainSize ),
        bcGlo_( bcGlo ),
        bcLoc_( bcLoc ) {}

  Teuchos::RCP<const DomainSize<Scalar> >      getDomainSize() const { return( domainSize_ ); }
  Teuchos::RCP<const BoundaryConditionsGlobal<dim> > getBCGlobal()   const { return( bcGlo_ ); }
  Teuchos::RCP<const BoundaryConditionsLocal>  getBCLocal()    const { return( bcLoc_ ); }

protected:

  Teuchos::RCP<const DomainSize<Scalar> > domainSize_;

  Teuchos::RCP<const BoundaryConditionsGlobal<dim> > bcGlo_;
  Teuchos::RCP<const BoundaryConditionsLocal> bcLoc_;


};



/// \relates Domain
template<class S, int d>
Teuchos::RCP<const Domain<S,d> > createDomain(
    const Teuchos::RCP<const DomainSize<S> >& domainSize,
    const Teuchos::RCP<const BoundaryConditionsGlobal<d> > bcGlo,
    const Teuchos::RCP<const BoundaryConditionsLocal> bcLoc ) {
  return(
      Teuchos::rcp(
          new Domain<S,d>( domainSize, bcGlo, bcLoc ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DOMAIN_HPP
