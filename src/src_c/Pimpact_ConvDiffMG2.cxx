#include <cmath>

#include "Pimpact_Operator.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_MultiGrid.hpp"

using S = double;
using O = int;
const int d = 3;
const int dNC=4;

using FSpaceT = Pimpact::Space<S,O,3,d,dNC>;
using CSpaceT = Pimpact::Space<S,O,3,d,2>;

using CS = Pimpact::CoarsenStrategy<FSpaceT,CSpaceT>;

template<class T1,class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;



template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;





int main( int argi, char** argv ) {

	// intialize MPI
	MPI_Init( &argi, &argv );

	auto pl = Teuchos::parameterList();
	pl->set("npy",1);
	pl->set("npx",1);

	pl->set( "domain", 1);

	//int nwinds = 360*2;
	//int nwinds = 360;
	//int nwinds = 180;
	//	int nwinds = 90;
	//	int nwinds = 128;
	//	int nwinds = 64;
	//	int nwinds = 32;
	//	int nwinds = 16;
	//	int nwinds = 8;
	//	int nwinds = 4;
	//	int nwinds = 2;
	int nwinds = 1;

	//  S pi = (S)4. * std::atan( (S)1. ) ;

	for( S re=1.; re<1e6; re*=10 ) {
		//	pl->set<S>( "Re", 100000 );
		//	pl->set<S>( "Re", 10000 );
		pl->set<S>( "Re", re );
		//	pl->set<S>( "Re", 100 );
		//	pl->set<S>( "Re", 10 );
		//	pl->set<S>( "Re", 1 );
		//	pl->set<S>( "Re", 0.1 );
		//	pl->set<S>( "Re", 0.01 );
		//	pl->set<S>( "Re", 0.001 );
		//	pl->set<S>( "Re", 0.0001 );


		pl->set<O>( "nx", 513 );
		pl->set<O>( "ny", 513 );
		//	pl->set<O>( "nx", 257 );
		//	pl->set<O>( "ny", 257 );
		//	pl->set<O>( "nx", 129 );
		//	pl->set<O>( "ny", 129 );
		//	pl->set<O>( "nx", 65 );
		//	pl->set<O>( "ny", 65 );
		//	pl->set<O>( "nx", 33 );
		//	pl->set<O>( "ny", 33 );
		//pl->set<O>( "nx", 17 );
		//pl->set<O>( "ny", 17 );


		auto space = Pimpact::create< Pimpact::Space<S,O,3,d,dNC> >( pl );

		auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>( space, 5 );

		auto wind = Pimpact::create<Pimpact::VectorField>( space );
		auto y = Pimpact::create<Pimpact::VectorField>( space );
		auto z = Pimpact::create<Pimpact::VectorField>( space );
		auto z2 = Pimpact::create<Pimpact::VectorField>( space );


		auto op = Pimpact::create<ConvDiffOpT>( space );

		for(short int dirx=1; dirx<4; dirx+=2 ) {
			for(short int diry=1; diry<2; diry+=2 ) {


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
					//ConvDiffJT,
					ConvDiffSORT,
					//ConvDiffSORT
					MOP > ( mgSpaces, pls );

				std::ofstream phifile;

				if( space()->rankST()==0 ) {
					std::string fname = "blaphin.txt";
					if( 3==dirx )
						fname.insert( 4, std::to_string( (long long)8 ) );
					else
						fname.insert( 4, std::to_string( (long long)dirx+diry*2+3 ) );
					phifile.open( fname, std::ofstream::out | std::ofstream::app );
				}

				for( int phii=0; phii<nwinds; ++phii ) {

					if( space()->rankST()==0 )
						phifile << re << "\t";

					// init solution
					y->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inX );
					y->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inY );

					auto sol = y->clone( Pimpact::ECopy::Deep );

					wind->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inX );
					wind->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inY );

					op->assignField( *wind );
					smoother->assignField( *wind );

					z->initField();

					// constructing rhs
					op->apply( *y, *z );
					{
						y->init(0);
						auto bc = z->clone( Pimpact::ECopy::Shallow );
						op->apply( *y, *bc );
						z->add( 1., *z, -1., *bc );
					}

					y->initField();

					std::ofstream ofs;
					std::string filename = "MG.txt";
					if( space()->rankST()==0 ) {
						if( 3==dirx )
							filename.insert( 2, std::to_string( (long long)8 ) );
						else
							filename.insert( 2, std::to_string( (long long)dirx+diry*2+3 ) );
					}

					if( space()->rankST()==0 )
						ofs.open(filename, std::ofstream::out);

					S error;
					int iter=0;
					do {

						smoother->apply( *z, *y );

						z2->add( -1, *sol, 1, *y );

						error = z2->norm()/sol->norm();

						if( space()->rankST()==0 ) ofs << error << "\n";
						if( space()->rankST()==0 ) std::cout <<"iter: " <<iter <<" " << error << "\n";

						iter++;
						if( iter>1000) error=-1;
						if( error>1e12){
							error=-1;
							iter=1000;
						}	
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
	}

	MPI_Finalize();
	return( 0 );

}
