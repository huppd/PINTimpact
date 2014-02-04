#pragma once
#ifndef PIMPACT_PROCGRIDSIZE_HPP
#define PIMPACT_PROCGRIDSIZE_HPP

#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"

extern "C" {
	void fsetPGS(const int& np1, const int& np2, const int& np3 );
}

namespace Pimpact{

template<class Ordinal>
class ProcGridSize {
public:
	typedef const Teuchos::Tuple<Ordinal,3> TO3;

protected:
	TO3 procGridSize_;
public:
	ProcGridSize( Ordinal np1, Ordinal np2, Ordinal np3 ):
		procGridSize_( Teuchos::tuple(np1, np2, np3) ) {};

	ProcGridSize( TO3 procGridSize ):
		procGridSize_( procGridSize ) {};

	void set_Impact(){
		fsetPGS( procGridSize_[0], procGridSize_[1], procGridSize_[2] );
	};

	void print( std::ostream& out ) {
		out  << " \tnpx=" << procGridSize_[0]
		     << " \tnpy=" << procGridSize_[1]
		     << " \tnpz=" << procGridSize_[2] << "\n";
	};
};


template<class Ordinal>
Teuchos::RCP<ProcGridSize<Ordinal> > createProcGridSize( Ordinal np1=1, Ordinal np2=1, Ordinal np3=1 ) {
	return Teuchos::rcp( new ProcGridSize<Ordinal>( np1, np2, np3 ) );
}

} // namespace Pimpact

#endif // PIMPACT_PROCGRIDSIZE_HPP
