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


int Pimpact::getDir1( const int& dir ) {

	int dir1 = ( dir + 1 )%3;
	int dir2 = ( dir + 2 )%3;
	if( dir2<dir1 )
		return( dir2 );
	return( dir1 );

}

int Pimpact::getDir2( const int& dir ) {

	int dir1 = ( dir + 1 )%3;
	int dir2 = ( dir + 2 )%3;
	if( dir2>dir1 )
		return( dir2 );
	return( dir1 );

}


Teuchos::RCP<std::ostream>
Pimpact::createOstream( const std::string& fname, int rank ) {

	if( 0==rank )
		return( Teuchos::rcp( new std::ofstream( fname ) ) );
	else
		return( Teuchos::rcp( new Teuchos::oblackholestream ) );

}

void Pimpact::setBoundaryConditions( const
		Teuchos::RCP<Teuchos::ParameterList>& pl , int dtype ) {

	switch( static_cast<Pimpact::EDomainType>(dtype) ) {
		case Pimpact::AllDirichlet:
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Z", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Z", Pimpact::DirichletBC );
			break;
		case Pimpact::AllPeriodic:
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower X", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper X", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Y", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Y", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Z", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Z", Pimpact::PeriodicBC );
			break;
		case Pimpact::AllNeumann:
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower X", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper X", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Y", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Y", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Z", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Z", Pimpact::NeumannBC );
			break;
		case Pimpact::AllSymmetric:
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower X", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper X", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Y", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Y", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Z", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Z", Pimpact::SymmetryBC );
			break;
		case Pimpact::Dirichelt2DChannel:
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Z", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Z", Pimpact::PeriodicBC );
			break;
		case Pimpact::Periodic2DChannel:
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower X", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper X", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Z", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Z", Pimpact::PeriodicBC );
			break;
		case Pimpact::Open2DChannel:
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper X", Pimpact::NeumannBC   );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "lower Z", Pimpact::PeriodicBC  );
			pl->sublist("boundary conditions").set<Pimpact::EBCType>( "upper Z", Pimpact::PeriodicBC  );
			break;
		default:
			std::cout << "!!!Warning: unkown EDomainType:\t" <<dtype<<"\t!!!\n";
	}
}

#ifdef COMPILE_ETI
template class Pimpact::Space<double,int,3,2>;
template class Pimpact::Space<double,int,3,4>;
template class Pimpact::Space<double,int,4,2>;
template class Pimpact::Space<double,int,4,4>;
#endif
