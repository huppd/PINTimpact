#include "Pimpact_Space.hpp"



std::string Pimpact::toString( Pimpact::ECoord type ) {
	switch( type ) {
		case Pimpact::ECoord::X :
			return( "X" );
			break;
		case Pimpact::ECoord::Y :
			return( "Y" );
			break;
		case Pimpact::ECoord::Z :
			return( "Z" );
			break;
		case Pimpact::ECoord::T :
			return( "T" );
			break;
	}
	return( "" ); // prevent compiler warning
}


std::string Pimpact::toString( Pimpact::EField type ) {
	switch( type ) {
		case Pimpact::EField::U :
			return( "U" );
			break;
		case Pimpact::EField::V :
			return( "V" );
			break;
		case Pimpact::EField::W :
			return( "W" );
			break;
		case Pimpact::EField::S :
			return( "S" );
			break;
	}
	return( "" ); // prevent compiler warning
}


std::string Pimpact::toString( Pimpact::EScalarField type ) {
	switch( type ) {
		case Pimpact::EScalarField::ConstField :
			return( "constant" );
			break;
		case Pimpact::EScalarField::Grad2D_inX :
			return( "grad in x" );
			break;
		case Pimpact::EScalarField::Grad2D_inY :
			return( "grad in y" );
			break;
		case Pimpact::EScalarField::Grad2D_inZ :
			return( "grad in z" );
			break;
		case Pimpact::EScalarField::Poiseuille2D_inX :
			return( "poiseuille in x" );
			break;
		case Pimpact::EScalarField::Poiseuille2D_inY :
			return( "poiseuille in y" );
			break;
		case Pimpact::EScalarField::Poiseuille2D_inZ :
			return( "poiseuille in z" );
			break;
		case Pimpact::EScalarField::FPoint :
			return( "poiseuille in z" );
			break;
	}
	return( "" ); // prevent compiler warning
}


Teuchos::RCP<std::ostream>
Pimpact::createOstream( const std::string& fname, int rank ) {

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
