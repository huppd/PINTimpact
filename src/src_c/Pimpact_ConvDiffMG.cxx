#include <cmath>

#include "Pimpact_Operator.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"

using S = double;
using O = int;
const int d = 3;
const int dNC=4;

using SpaceT = Pimpact::Space<S,O,3,d,dNC>;

using FSpaceT = SpaceT;
using CSpaceT = Pimpact::Space<S,O,3,d,2>;


using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;

template<class T1,class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;



template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;





int main( int argi, char** argv ) {

	// intialize MPI
	MPI_Init( &argi, &argv );

	auto pl = Teuchos::parameterList();
	pl->set("npy",2);
	pl->set("npx",2);

	pl->set( "domain", 1);

	int nwinds = 64;

	S pi = 4. * std::atan(1.) ;

	pl->set<S>( "Re", 1000 );

	pl->set<O>( "nx", 257 );
	pl->set<O>( "ny", 257 );

	auto space = Pimpact::create<SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, 5 );

	Pimpact::VectorField<SpaceT> wind( space );
	Pimpact::VectorField<SpaceT> y   ( space );
	Pimpact::VectorField<SpaceT> z   ( space );
	Pimpact::VectorField<SpaceT> z2  ( space );

	auto op = Pimpact::create<ConvDiffOpT>( space );

	{
		short int dirx=1;
		short int diry=1;
		{

			auto pls = Teuchos::parameterList();
			pls->sublist("Smoother").set( "omega", 1. );
			pls->sublist("Smoother").set( "numIters", ((dirx==3)?1:4)*1 );
			pls->sublist("Smoother").set<int>( "Ordering", (dirx==3)?1:0 );
			pls->sublist("Smoother").set<short int>( "dir X", dirx );
			pls->sublist("Smoother").set<short int>( "dir Y", diry );
			pls->sublist("Smoother").set<short int>( "dir Z", 1 );

			auto smoother =
				Pimpact::createMultiGrid<
				Pimpact::VectorField,
				TransVF,
				RestrVF,
				InterVF,
				ConvDiffOpT,
				ConvDiffOpT,
				ConvDiffSORT,
				MOP
					> ( mgSpaces, pls );

			std::ofstream phifile;

			if( space()->rankST()==0 ) {
				std::string fname = "phin.txt";
				if( 3==dirx )
					fname.insert( 4, std::to_string( (long long)8 ) );
				else
					fname.insert( 4, std::to_string( (long long)dirx+diry*2+3 ) );
				phifile.open( fname, std::ofstream::out);
			}

			for( int phii=0; phii<nwinds; ++phii ) {

				S phi = 2.*pi*phii/(nwinds);

				if( space()->rankST()==0 )
					phifile << phi << "\t";

				// init solution
				y(Pimpact::U).initField( Pimpact::Grad2D_inX );
				y(Pimpact::V).initField( Pimpact::Grad2D_inY );

				auto sol = y.clone( Pimpact::ECopy::Deep );
				//sol->write(3333);

				wind(Pimpact::U).initField(  Pimpact::ConstField, std::cos( phi ) );
				wind(Pimpact::V).initField(  Pimpact::ConstField, std::sin( phi ) );

				op->assignField( wind );
				smoother->assignField( wind );

				z.initField();

				// constructing rhs
				op->apply( y, z );
				{
					y.init(0);
					auto bc = z.clone( Pimpact::ECopy::Shallow );
					op->apply( y, *bc );
					z.add( 1., z, -1., *bc );
				}
				//z->write(2222);

				y.initField();

				std::ofstream ofs;
				std::string filename = "GS.txt";
				filename.insert( 2, std::to_string( (long long)phii) );

				if( space()->rankST()==0 )
					ofs.open(filename, std::ofstream::out);

				S error;
				int iter=0;
				do {

					smoother->apply( z, y );

					z2.add( -1., *sol, 1., y );
					//					op->apply( *y, *z2 );
					//					z2->add( 1., *z2, -1., *z );

					error = z2.norm()/sol->norm();
					//          error = z2->norm();


					if( space()->rankST()==0 ) ofs << error << "\n";
					if( space()->rankST()==0 ) std::cout <<"iter: " <<iter <<" " << error << "\n";

					iter++;
					if( iter>1000) error=-1;

				}
				while( error>1.e-6 );

				if( space()->rankST()==0 )
					//					phifile << error << "\n";
					phifile << iter << "\n";


				if( space()->rankST()==0 )
					ofs.close();
			}
			if( space()->rankST()==0 )
				phifile.close();
		}

	}

	MPI_Finalize();
	return( 0 );

}
