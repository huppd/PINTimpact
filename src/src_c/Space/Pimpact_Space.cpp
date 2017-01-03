#include "Pimpact_Space.hpp"



std::string Pimpact::toString( Pimpact::ECoord type ) {
	switch( type ) {
		case Pimpact::ECoord::X : return( "X" );
		case Pimpact::ECoord::Y : return( "Y" );
		case Pimpact::ECoord::Z : return( "Z" );
		case Pimpact::ECoord::T : return( "T" );
	}
	//return( "" ); // prevent compiler warning
}


Pimpact::ECoord& operator++( Pimpact::ECoord& c ) {
	switch( c ) {
		case Pimpact::ECoord::X   : return( c = Pimpact::ECoord::Y );
		case Pimpact::ECoord::Y   : return( c = Pimpact::ECoord::Z );
		case Pimpact::ECoord::Z   : return( c = Pimpact::ECoord::end );
		case Pimpact::ECoord::end : return( c = Pimpact::ECoord::end );
		case Pimpact::ECoord::T   : return( c = Pimpact::ECoord::end );
	}
}

bool operator<( const Pimpact::ECoord& c, const int& i ) {
	return( static_cast<int>(c)<i );
}


std::string Pimpact::toString( Pimpact::EField type ) {
	switch( type ) {
		case Pimpact::EField::U : return( "U" );
		case Pimpact::EField::V : return( "V" );
		case Pimpact::EField::W : return( "W" );
		case Pimpact::EField::S : return( "S" );
	}
	//return( "" ); // prevent compiler warning
}


Pimpact::EField& operator++( Pimpact::EField& c ) {
	switch( c ) {
		case Pimpact::EField::U   : return( c = Pimpact::EField::V );
		case Pimpact::EField::V   : return( c = Pimpact::EField::W );
		//case Pimpact::EField::W   : return( c = Pimpact::EField::end );
		//case Pimpact::EField::end : return( c = Pimpact::EField::end );
		//case Pimpact::EField::S   : return( c = Pimpact::EField::end );
	}
}

bool operator<( const Pimpact::EField& c, const int& i ) {
	return( static_cast<int>(c)<i );
}

std::string Pimpact::toString( Pimpact::EScalarField type ) {
	switch( type ) {
		case Pimpact::EScalarField::ConstField : return( "constant" );
		case Pimpact::EScalarField::Grad2D_inX : return( "grad in x" );
		case Pimpact::EScalarField::Grad2D_inY : return( "grad in y" );
		case Pimpact::EScalarField::Grad2D_inZ : return( "grad in z" );
		case Pimpact::EScalarField::Poiseuille2D_inX : return( "poiseuille in x" );
		case Pimpact::EScalarField::Poiseuille2D_inY : return( "poiseuille in y" );
		case Pimpact::EScalarField::Poiseuille2D_inZ : return( "poiseuille in z" );
		case Pimpact::EScalarField::FPoint : return( "poiseuille in z" );
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
			pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::DirichletBC );
			break;
		case Pimpact::AllPeriodic:
			pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::PeriodicBC );
			break;
		case Pimpact::AllNeumann:
			pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::NeumannBC );
			pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::NeumannBC );
			break;
		case Pimpact::AllSymmetric:
			pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::SymmetryBC );
			pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::SymmetryBC );
			break;
		case Pimpact::Dirichelt2DChannel:
			pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::PeriodicBC );
			break;
		case Pimpact::Periodic2DChannel:
			pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::PeriodicBC );
			pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::PeriodicBC );
			break;
		case Pimpact::Open2DChannel:
			pl->sublist("boundary conditions").set<int>( "lower X", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "upper X", Pimpact::NeumannBC   );
			pl->sublist("boundary conditions").set<int>( "lower Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "upper Y", Pimpact::DirichletBC );
			pl->sublist("boundary conditions").set<int>( "lower Z", Pimpact::PeriodicBC  );
			pl->sublist("boundary conditions").set<int>( "upper Z", Pimpact::PeriodicBC  );
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
