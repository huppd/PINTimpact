#include <cmath>

#include "Pimpact_Operator.hpp"

typedef double S;
typedef int O;
const int d = 3;

template<class T> using ConvDiffOpT = Pimpact::ConvectionVOp<Pimpact::ConvectionDiffusionSOp<T> >;

int main( int argi, char** argv ) {

  // intialize MPI
  MPI_Init( &argi, &argv );

  auto pl = Teuchos::parameterList();
  pl->set("npy",1);
  pl->set("npx",1);

  pl->set( "domain", 1);

  //  int nwinds = 360*2;
//  int nwinds = 360;
    int nwinds = 360/2;
//    int nwinds = 360/4;
//    int nwinds = 360/6;
//    int nwinds = 64;
//    int nwinds = 32;
//    int nwinds = 16;
//    int nwinds = 8;
//    int nwinds = 4;
//    int nwinds = 1;

  S pi = (S)4. * std::atan( (S)1. ) ;

  pl->set<S>( "Re", 10000 );
//  pl->set<S>( "lx", 1 );
//  pl->set<S>( "ly", 1 );
//    pl->set<O>( "nx", 513 );
//    pl->set<O>( "ny", 513 );
//    pl->set<O>( "nx", 257 );
//    pl->set<O>( "ny", 257 );
    pl->set<O>( "nx", 129 );
    pl->set<O>( "ny", 129 );
//    pl->set<O>( "nx", 65 );
//    pl->set<O>( "ny", 65 );
//  pl->set<O>( "nx", 17 );
//  pl->set<O>( "ny", 17 );
  auto space = Pimpact::createSpace<S,O,d,2>( pl );


  auto wind = Pimpact::create<Pimpact::VectorField>( space );
  auto y = Pimpact::create<Pimpact::VectorField>( space );
  auto z = Pimpact::create<Pimpact::VectorField>( space );
  auto z2 = Pimpact::create<Pimpact::VectorField>( space );


  auto op = Pimpact::create<ConvDiffOpT>( space );

  for(short int dirx=-1; dirx<4; dirx+=2 ) {
    for(short int diry=-1; diry<2; diry+=2 ) {

      if( 3==dirx && diry==1 ) break;

      auto pls = Teuchos::parameterList();
      pls->set( "omega", 1. );
      pls->set( "numIter", 1 );
      pls->set<int>( "Ordering", (dirx==3)?1:0 );
      pls->set<short int>( "dir X", dirx );
      pls->set<short int>( "dir Y", diry );
      pls->set<short int>( "dir Z", 1 );

      auto smoother =
           Pimpact::create<
             Pimpact::ConvectionVSmoother<
               ConvDiffOpT<Pimpact::Space<S,O,d,2> > ,
               Pimpact::ConvectionDiffusionSORSmoother > > (
                   op,
                   pls );

      std::ofstream phifile;

      if( space()->rankST()==0 ) {
        std::string fname = "phin.txt";
        if( 3==dirx )
          fname.insert( 4, std::to_string( 8 ) );
        else
          fname.insert( 4, std::to_string( dirx+diry*2+3 ) );
        phifile.open( fname, std::ofstream::out);
      }

      for( int phii=0; phii<nwinds; ++phii ) {

        S phi = phii*2.*pi/(nwinds);
//        S phi = 2.*pi/(nwinds)*0.25;
//        S phi = 2.*pi/(nwinds)*0.;
//        S phi = phi*2.*pi/(nwinds);

        if( space()->rankST()==0 )
          phifile << phi << "\t";


        // init solution
        y->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
        y->getFieldPtr(1)->initField( Pimpact::Grad2D_inX );

        auto sol = y->clone( Pimpact::DeepCopy );

        wind->initField( Pimpact::ConstFlow, std::cos( phi ), std::sin( phi ), 0. );
//        wind->write(1111);

        z->initField( Pimpact::ConstFlow, 0., 0., 0. );

        op->assignField( *wind );
//        smoother->assignField( *wind );

        // constructing rhs
        op->apply( *y, *z );
        {
           y->init(0);
           auto bc = z->clone( Pimpact::ShallowCopy );
           op->apply( *y, *bc );
           z->add( 1., *z, -1., *bc );
        }
//        z->write(2222);

        y->initField( Pimpact::ConstFlow, 0., 0., 0. );

        std::ofstream ofs;
        std::string filename = "GS.txt";
        filename.insert( 2, std::to_string(phii) );

        if( space()->rankST()==0 )
          ofs.open(filename, std::ofstream::out);

        S error;
        int iter=0;
        do {

          smoother->apply( *z, *y, 1. );

          z2->add( -1, *sol, 1, *y );

          error = z2->norm()/sol->norm();
//          error = z2->norm();

          if( iter>1000) error=-1;

          if( space()->rankST()==0 ) ofs << error << "\n";

          iter++;

        }
        while( error>1.e-6 );

        if( space()->rankST()==0 )
//          phifile << error << "\n";
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
