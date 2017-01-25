#include <cmath>

#include "Pimpact_Operator.hpp"

using S = double;
using O = int;
const int sd = 2;
const int d = 3;
const int dNC = 2;

using SpaceT = Pimpact::Space<S,O,sd,d,dNC>;

template<class T> using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;


int main( int argi, char** argv ) {

	// intialize MPI
	MPI_Init( &argi, &argv );

	auto pl = Teuchos::parameterList();
	pl->set("npy",1);
	pl->set("npx",1);

	int nwinds = 64;

	S pi = (S)4. * std::atan( (S)1. ) ;

	pl->set<S>( "Re", 1000 );
	pl->set<O>( "nx", 129 );
	pl->set<O>( "ny", 129 );


	auto space = Pimpact::create< SpaceT >( pl );

	Pimpact::VectorField<SpaceT> wind( space );
	Pimpact::VectorField<SpaceT> y   ( space );
	Pimpact::VectorField<SpaceT> z   ( space );
	Pimpact::VectorField<SpaceT> z2  ( space );


	auto op = Pimpact::create<ConvDiffOpT>( space );

	for(short int dirx=-1; dirx<4; dirx+=2 ) {
		for(short int diry=-1; diry<2; diry+=2 ) {

			if( 3==dirx && diry==1 ) break;

			auto pls = Teuchos::parameterList();
			pls->set( "omega", 1. );
			pls->set( "numIters", 1 );
			pls->set<int>( "Ordering", (dirx==3)?1:0 );
			pls->set<short int>( "dir X", dirx );
			pls->set<short int>( "dir Y", diry );
			pls->set<short int>( "dir Z", 1 );

			auto smoother =
				Pimpact::create<
					Pimpact::NonlinearSmoother<
						ConvDiffOpT<SpaceT> ,
						Pimpact::ConvectionDiffusionSORSmoother > > (
						op,
						pls );

			std::ofstream phifile;

			if( space()->rankST()==0 ) {
				std::string fname = "phin.txt";
				if( 3==dirx )
					fname.insert( 4, std::to_string( static_cast<long long>(8) ) );
				else
					fname.insert( 4, std::to_string( static_cast<long long>(dirx+diry*2+3) ) );
				phifile.open( fname, std::ofstream::out);
			}

			for( int phii=0; phii<nwinds; ++phii ) {

				S phi = phii*2.*pi/(nwinds);

				if( space()->rankST()==0 )
					phifile << phi << "\t";

				// init solution
				y(Pimpact::F::U).initField( Pimpact::Grad2D_inY );
				y(Pimpact::F::V).initField( Pimpact::Grad2D_inX );

				auto sol = y.clone( Pimpact::ECopy::Deep );

				wind(Pimpact::F::U).init( std::cos( phi ) );
				wind(Pimpact::F::V).init( std::sin( phi ) );

				op->assignField( wind );

				z.init();

				// constructing rhs
				op->apply( y, z );

				y.init();

				std::ofstream ofs;
				std::string filename = "GS.txt";
				filename.insert( 2, std::to_string( static_cast<long long>(phii) ) );

				if( space()->rankST()==0 )
					ofs.open(filename, std::ofstream::out);

				S error;
				int iter=0;
				do {

					smoother->apply( z, y );

					z2.add( -1, *sol, 1, y );

					error = z2.norm()/sol->norm();

					if( iter>1000) error=-1;

					if( space()->rankST()==0 ) ofs << error << "\n";
					if( space()->rankST()==0 ) std::cout << error << "\n";

					iter++;

				}
				while( error>1.e-6 );

				if( space()->rankST()==0 )
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
