#include "Pimpact_Space.hpp"

Teuchos::RCP<std::ostream>
Pimpact::createOstream( const std::string& fname, int rank=0 ) {

	if( 0==rank )
		return( Teuchos::rcp( new std::ofstream( fname ) ) );
	else
		return( Teuchos::rcp( new Teuchos::oblackholestream ) );

}

#ifdef COMPILE_ETI
template class Pimpact::Space<double,int,3,2>;
template class Pimpact::Space<double,int,3,4>;
template class Pimpact::Space<double,int,4,2>;
template class Pimpact::Space<double,int,4,4>;
#endif
