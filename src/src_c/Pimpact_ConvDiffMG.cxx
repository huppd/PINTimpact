#include <cmath>

#include "Pimpact_Operator.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"

using S = double;
using O = int;
const int sd = 2;
const int d  = 3;
const int dNC= 4;

using GridT = Pimpact::Grid<S, O, sd, d, dNC>;

using FGridT = GridT;
using CGridT = Pimpact::Grid<S, O, sd, d, 2>;


using CS = Pimpact::CoarsenStrategyGlobal<FGridT, CGridT>;

template<class T1, class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1, T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;



template<class T> using MOP = Pimpact::InverseOp<T>;





int main(int argi, char** argv) {

  // intialize MPI
  MPI_Init(&argi, &argv);

  auto pl = Teuchos::parameterList();
  pl->set("npy", 2);
  pl->set("npx", 2);

  int nwinds = 64;

  S pi = 4. * std::atan(1.) ;

  pl->set<S>("Re", 1000);

  pl->set<O>("nx", 129);
  pl->set<O>("ny", 129);

  auto grid = Pimpact::create<GridT>(pl);

  auto mgGrids = Pimpact::createMGGrids<CS>(grid, 5);

  Pimpact::VectorField<GridT> wind(grid);
  Pimpact::VectorField<GridT> y   (grid);
  Pimpact::VectorField<GridT> z   (grid);
  Pimpact::VectorField<GridT> z2  (grid);

  auto op = Pimpact::create<ConvDiffOpT>(grid);

  {
    short int dirx=1;
    short int diry=1;
    {

      auto pls = Teuchos::parameterList();
      pls->sublist("Smoother").set("omega", 1.);
      pls->sublist("Smoother").set("numIters", ((dirx==3)?1:4)*1);
      pls->sublist("Smoother").set<int>("Ordering", (dirx==3)?1:0);
      pls->sublist("Smoother").set<short int>("dir X", dirx);
      pls->sublist("Smoother").set<short int>("dir Y", diry);
      pls->sublist("Smoother").set<short int>("dir Z", 1);

      auto smoother = Pimpact::createMultiGrid<
        Pimpact::VectorField,
        TransVF,
        RestrVF,
        InterVF,
        ConvDiffOpT,
        ConvDiffOpT,
        ConvDiffSORT,
        MOP > (mgGrids, op, pls);

      std::ofstream phifile;

      if(grid()->rankST()==0) {
        std::string fname = "phin.txt";
        if(3==dirx)
          fname.insert(4, std::to_string((long long)8));
        else
          fname.insert(4, std::to_string((long long)dirx+diry*2+3));
        phifile.open(fname, std::ofstream::out);
      }

      for(int phii=0; phii<nwinds; ++phii) {

        S phi = 2.*pi*phii/(nwinds);

        if(grid()->rankST()==0)
          phifile <<phi <<"\t";

        // init solution
        y(Pimpact::F::U).initField(Pimpact::Grad2D_inX);
        y(Pimpact::F::V).initField(Pimpact::Grad2D_inY);

        auto sol = y.clone(Pimpact::ECopy::Deep);
        //sol->write(3333);

        wind(Pimpact::F::U).init(std::cos(phi));
        wind(Pimpact::F::V).init(std::sin(phi));

        op->assignField(wind);
        smoother->assignField(wind);

        z.init();

        // constructing rhs
        op->apply(y, z);

        y.init();

        std::ofstream ofs;
        std::string filename = "GS.txt";
        filename.insert(2, std::to_string((long long)phii));

        if(grid()->rankST()==0)
          ofs.open(filename, std::ofstream::out);

        S error;
        int iter=0;
        do {

          smoother->apply(z, y);

          z2.add(-1., *sol, 1., y);

          error = z2.norm()/sol->norm();


          if(grid()->rankST()==0) ofs <<error <<"\n";
          if(grid()->rankST()==0) std::cout <<"iter: " <<iter <<" " <<error <<"\n";

          iter++;
          if(iter>1000) error=-1;

        } while(error>1.e-6);

        if(grid()->rankST()==0)
          //					phifile <<error <<"\n";
          phifile <<iter <<"\n";


        if(grid()->rankST()==0)
          ofs.close();
      }
      if(grid()->rankST()==0)
        phifile.close();
    }

  }

  MPI_Finalize();
  return 0;

}
