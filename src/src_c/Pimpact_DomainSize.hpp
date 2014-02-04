#pragma once
#ifndef PIMPACT_DOMAINSIZE_HPP
#define PIMPACT_DOMAINSIZE_HPP

#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"

extern "C" {
	void fsetDS(const double& L1, const double& L2, const double& L3 );
}

namespace Pimpact{

template<class Scalar>
class DomainSize {
public:
	typedef const Teuchos::Tuple<Scalar,3> TS3;

protected:
	TS3 domainSize_;
public:
	DomainSize( Scalar L1, Scalar L2, Scalar L3 ):
		domainSize_( Teuchos::tuple(L1, L2, L3) ) {};

	DomainSize( TS3 domainSize ):
		domainSize_( domainSize ) {};

	void set_Impact(){
		fsetDS( domainSize_[0], domainSize_[1], domainSize_[2] );
	};

	void print( std::ostream& out ) {
		out << " \tlx=" << domainSize_[0]
	      << " \tly=" << domainSize_[1]
	      << " \tlz=" << domainSize_[2] << "\n";
	};
};


template<class Scalar>
Teuchos::RCP<DomainSize<Scalar> > createDomainSize( Scalar L1=1, Scalar L2=1, Scalar L3=1 ) {
	return Teuchos::rcp( new DomainSize<Scalar>( L1, L2, L3 ) );
}

} // namespace Pimpact

#endif // PIMPACT_DOMAINSIZE_HPP
