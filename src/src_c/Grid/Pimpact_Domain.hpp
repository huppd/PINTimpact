#pragma once
#ifndef PIMPACT_DOMAIN_HPP
#define PIMPACT_DOMAIN_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_DomainSize.hpp"
#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"



// \defgroup domain Domain
///
/// overloaded class managing indexing, grid ...


namespace Pimpact{


/// \ingroup domain
template<class Scalar>
class Domain {

public:

  Domain(
      const Teuchos::RCP<const DomainSize<Scalar> >& domainSize,
      const Teuchos::RCP<const BoundaryConditionsGlobal> bcGlo,
      const Teuchos::RCP<const BoundaryConditionsLocal> bcLoc ):
        domainSize_( domainSize ),
        bcGlo_( bcGlo ),
        bcLoc_( bcLoc ) {}

  Teuchos::RCP<const DomainSize<Scalar> >      getDomainSize() const { return( domainSize_ ); }
  Teuchos::RCP<const BoundaryConditionsGlobal> getBCGlobal()   const { return( bcGlo_ ); }
  Teuchos::RCP<const BoundaryConditionsLocal>  getBCLocal()    const { return( bcLoc_ ); }

protected:

  Teuchos::RCP<const DomainSize<Scalar> > domainSize_;

  Teuchos::RCP<const BoundaryConditionsGlobal> bcGlo_;
  Teuchos::RCP<const BoundaryConditionsLocal> bcLoc_;


};


///// \relates Domain
//template<class S=double>
//Teuchos::RCP<const Domain<S> > createDomain() {
//  return(
//      Teuchos::rcp(
//          new Domain<S>(
//              createDomainSize<S>(),
//              createBoudaryConditionsGlobal(),
//              createBoudaryConditionsLocal() ) ) );
//
//}


/// \relates Domain
template<class S>
Teuchos::RCP<const Domain<S> > createDomain(
    const Teuchos::RCP<const DomainSize<S> >& domainSize,
    const Teuchos::RCP<const BoundaryConditionsGlobal> bcGlo,
    const Teuchos::RCP<const BoundaryConditionsLocal> bcLoc ) {
  return(
      Teuchos::rcp(
          new Domain<S>( domainSize, bcGlo, bcLoc ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DOMAIN_HPP
