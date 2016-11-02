#include <cmath>

#include "Pimpact_Operator.hpp"

using S = double;
using O = int;
const int d = 3;


template<class T> using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;


int main( int argi, char** argv ) {

  // intialize MPI
  MPI_Init( &argi, &argv );

  auto pl = Teuchos::parameterList();
  pl->set("npy",1);
  pl->set("npx",1);

  pl->set( "domain", 1);


  pl->set<S>( "Re", 10000 );
  //  pl->set<S>( "lx", 1 );
  //  pl->set<S>( "ly", 1 );
//  pl->set<O>( "nx", 513 );
//  pl->set<O>( "ny", 513 );
  //    pl->set<O>( "nx", 257 );
  //    pl->set<O>( "ny", 257 );
      pl->set<O>( "nx", 129 );
      pl->set<O>( "ny", 129 );
//      pl->set<O>( "nx", 65 );
//      pl->set<O>( "ny", 65 );
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
      //      pls->set<int>( "Ordering", 0 );
      pls->set<short int>( "dir X", dirx );
      pls->set<short int>( "dir Y", diry );
      //      pls->set<short int>( "dir X", 1 );
      //      pls->set<short int>( "dir Y", 1 );
      pls->set<short int>( "dir Z", 1 );
      //  pls->set( "numIters",10)

      auto smoother =
           Pimpact::create<
             Pimpact::NonlinearSmoother<
               ConvDiffOpT<Pimpact::Space<S,O,d,2> > ,
               Pimpact::ConvectionDiffusionSORSmoother > > (
                   op,
                   pls );

      std::ofstream phifile;

      if( space()->rankST()==0 ) {
        std::string fname = "raki.txt";
        if( 3==dirx )
          fname.insert( 4, std::to_string( (long long)8 ) );
        else
          fname.insert( 4, std::to_string( (long long)dirx+diry*2+3 ) );
        phifile.open( fname, std::ofstream::out);
      }


      // init solution
      y->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
      y->getFieldPtr(1)->initField( Pimpact::Grad2D_inX );

      auto sol = y->clone( Pimpact::ECopy::Deep );

      wind->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
      wind->getFieldPtr(1)->initField( Pimpact::Grad2D_inX );
      //        wind->write(1111);

      z->initField();

      op->assignField( *wind );

      // constructing rhs
      op->apply( *y, *z );
      {
        y->init(0);
        auto bc = z->clone( Pimpact::ECopy::Shallow );
        op->apply( *y, *bc );
        z->add( 1., *z, -1., *bc );
      }
      //        z->write(2222);

      y->initField();


      S error;
      int iter=0;
      do {

        smoother->apply( *z, *y );

        z2->add( -1, *sol, 1, *y );

        error = z2->norm()/sol->norm();

        if( iter>2000) error=-1;

        if( dirx==3 )
          iter+=4;
        else
          iter++;

        if( space()->rankST()==0 ) {
          phifile << iter << "\t" << error << "\n";
        }

      }
      while( error>1.e-6 );

    if( space()->rankST()==0 )
      phifile.close();
  }

}

MPI_Finalize();
return( 0 );

}
