//#include <cmath>

#include <mpi.h>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Grid.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_TransferModeOp.hpp"
#include "Pimpact_ModeSmoother.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_CoarsenStrategy.hpp"



// Grid types
using ST = double;
using OT = int;

//const int sd = 2;
const int sd = 3;
const int dNC = 4;
//const int dNC = 3;
//const int dNC = 2;

using GridT = Pimpact::Grid<ST, OT, sd, 4, dNC>;

using FGridT = Pimpact::Grid<ST, OT, sd, 4, dNC>;
using CGridT = Pimpact::Grid<ST, OT, sd, 4, 2>;


using MGGridsT = Pimpact::MGGrids<FGridT, CGridT>;

using CS = Pimpact::CoarsenStrategyGlobal<FGridT, CGridT>;
//using CS = Pimpact::CoarsenStrategy<FGridT, CGridT>; // dirty fix: till gather isn't fixed


// Field types
template<class T>
using ModeField = Pimpact::ModeField<Pimpact::VectorField<T>>;

using MGFieldsT = Pimpact::MGFields<MGGridsT, ModeField>;

// Operator types
template<class T>
using ModeOp = Pimpact::ModeNonlinearOp<ConvDiffOpT<T> >;

// Smoother types
template<class OpT>
using ModeSmootherT = Pimpact::ModeSmoother<OpT>;

template<class T>
using MOP = Pimpact::InverseOp<T>;

// Transfer types
template<class T1, class T2>
using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1, T2> >;
template<class T>
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T>
using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

template<class T1, class T2>
using ModeTransVF = Pimpact::TransferModeOp<Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1, T2> > >;
template<class T>
using ModeRestrVF = Pimpact::TransferModeOp<Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> > >;
template<class T>
using ModeInterVF = Pimpact::TransferModeOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> > >;







int main(int argi, char** argv) {

  // intialize MPI
  MPI_Init(&argi, &argv);
  {
    /////////////////////////////////////////// set up parameters ///////////////////////////
    Teuchos::CommandLineProcessor my_CLP;

    std::string xmlFilename = "parameter3D.xml";
    my_CLP.setOption("filename", &xmlFilename, "file name of the input xml parameterlist");

    int realCase = 0;
    my_CLP.setOption("realCase", &realCase, "real case");

    my_CLP.recogniseAllOptions(true);
    my_CLP.throwExceptions(true);

    my_CLP.parse(argi, argv);

    auto pl = Teuchos::getParametersFromXmlFile(xmlFilename);
    ////////////////////////////////////////// end of set up parameters /////////////////////////

    int withoutput=pl->sublist("Solver").get<int>("withoutput", 1);

    Teuchos::RCP<const GridT> grid = Pimpact::create<GridT>(
        Teuchos::sublist(pl, "Grid", true));

    Pimpact::ModeField<Pimpact::VectorField<GridT> > x(grid);
    Pimpact::ModeField<Pimpact::VectorField<GridT> > y(grid);
    Pimpact::ModeField<Pimpact::VectorField<GridT> > rhs(grid);
    Pimpact::ModeField<Pimpact::VectorField<GridT> > sol(grid);
    Pimpact::ModeField<Pimpact::VectorField<GridT> > err(grid);

    ST error;

    Teuchos::RCP<ConvDiffOpT<GridT>> zeroOp = Pimpact::create<ConvDiffOpT>(grid);

    Teuchos::RCP<ModeOp<GridT>> modeOp = Teuchos::rcp(new ModeOp<GridT>(zeroOp));

    pl->sublist("M_ConvDiff").sublist("Solver").set<Teuchos::RCP<std::ostream> >(
        "Output Stream",
        Pimpact::createOstream(modeOp->getLabel()+".txt", grid->rankST()));
    auto modeInv = Pimpact::createInverseOp(modeOp, Teuchos::sublist(pl, "M_ConvDiff"));

   
    Teuchos::RCP<const MGGridsT> mgGrids =
      Pimpact::createMGGrids<CS>(grid, pl->sublist("Multi Grid").get<int>("maxGrids"));

    //pl->sublist("ConvDiff").sublist("Solver").set<Teuchos::RCP<std::ostream> >("Output Stream",
        //Pimpact::createOstream(zeroOp->getLabel()+".txt", grid->rankST()));

    pl->sublist("ConvDiff").sublist("Solver").set<Teuchos::RCP<std::ostream> >(
        "Output Stream",
        Pimpact::createOstream(zeroOp->getLabel()+".txt",
          withoutput?grid->rankST():-1));

    auto zeroInv = Pimpact::createInverseOp(
        zeroOp, Teuchos::sublist(pl, "ConvDiff"));
    //auto mgConvDiff =
      //Pimpact::create<ConvDiffSORT>(
          //zeroOp,
          //Teuchos::sublist(Teuchos::sublist(Teuchos::sublist(pl, "ConvDiff"), "Multi Grid"), "Coarse Grid Solver")); 
    auto mgConvDiff = Pimpact::createMultiGrid<
      Pimpact::VectorField,
      TransVF,
      RestrVF,
      InterVF,
      ConvDiffOpT,
      ConvDiffOpT,
      //ConvDiffJT,
      ConvDiffSORT,
      ConvDiffSORT
      //ConvDiffJT
        > (mgGrids, zeroOp, Teuchos::sublist(Teuchos::sublist(pl, "ConvDiff"), "Multi Grid")) ;

    auto mgMConvDiff = Pimpact::createMultiGrid<
      ModeField,
      ModeTransVF,
      ModeRestrVF,
      ModeInterVF,
      ModeOp,
      ModeOp,
      ModeSmootherT,
      //ModeSmootherT
      //MOP,
      MOP 
        > (mgGrids, modeOp, Teuchos::sublist(Teuchos::sublist(pl, "M_ConvDiff"), "Multi Grid")) ;

    if(0==grid->rankST()) {
      mgConvDiff->print();
      mgMConvDiff->print();
    }

    std::string convDiffPrecString =
      pl->sublist("ConvDiff").get<std::string>("preconditioner", "right");
    if("right" == convDiffPrecString)
      zeroInv->setRightPrec(Pimpact::createMultiOperatorBase(mgConvDiff));
    if("left" == convDiffPrecString)
      zeroInv->setLeftPrec(Pimpact::createMultiOperatorBase(mgConvDiff));

    //{
      //std::string modeConvDiffPrecString =
        //pl->sublist("M_ConvDiff").get<std::string>("preconditioner", "right");

      //if("none"!=modeConvDiffPrecString) {
        //auto modePrec =
          //Pimpact::createMultiOperatorBase(
              //Pimpact::create<Pimpact::ModePrec>(
                //zeroInv,
                ////mgConvDiff,
                //Teuchos::sublist(Teuchos::sublist(pl, "M_ConvDiff"), "Mode prec")));

        //if("right" == modeConvDiffPrecString)
          //modeInv->setRightPrec(modePrec);

        //if("left" == modeConvDiffPrecString)
          //modeInv->setLeftPrec(modePrec);
      //}
    //}
    std::string modeConvDiffPrecString =
      pl->sublist("M_ConvDiff").get<std::string>("preconditioner", "right");
    if("right" == modeConvDiffPrecString)
      modeInv->setRightPrec(Pimpact::createMultiOperatorBase(mgMConvDiff));

    if("left" == modeConvDiffPrecString)
      modeInv->setLeftPrec(Pimpact::createMultiOperatorBase(mgMConvDiff));

    ST iRe = 1./grid->getDomainSize()->getRe();
    ST a2 = grid->getDomainSize()->getAlpha2()*iRe;

    ST lx = grid->getDomainSize()->getSize(Pimpact::X) ;
    ST ly = grid->getDomainSize()->getSize(Pimpact::Y) ;
    ST lz = grid->getDomainSize()->getSize(Pimpact::Z) ;
    std::cout << "a2: " << a2 << "\n";
    std::cout << "iRe: " << iRe << "\n";
    std::cout << "lx: " << lx << "\n";
    std::cout << "ly: " << ly << "\n";
    std::cout << "lz: " << lz << "\n";
    //
    // computing zero mode of z
    // set paramteters
    auto para = Teuchos::parameterList();
    para->set<ST>("mulI", a2 );
    para->set<ST>("mulC", 1. );
    para->set<ST>("mulL", iRe);
    modeOp->setParameter(para);
    modeInv->setParameter(para);

    // initializtion
    {
      Pimpact::VectorField<GridT> wind(grid);
      if(realCase)
        wind.initField(pl->sublist("Base flow").sublist("0 mode"));
      else {
        wind(Pimpact::F::U).init(1.);
        wind(Pimpact::F::V).init(1.);
        wind(Pimpact::F::W).init(1.);
      }
      //if(withoutput) wind.write(1000);

      zeroOp->assignField(wind);
      mgConvDiff->assignField(wind);
      auto modeWind = x.clone(Pimpact::ECopy::Shallow);
      modeWind->getCField() = wind;
      mgMConvDiff->assignField(*modeWind);
    }

    if(realCase) {
      rhs.getCField().initField(pl->sublist("Force").sublist("cos mode"));
      rhs.getSField().initField(pl->sublist("Force").sublist("sin mode"));
      std::cout << "||rhs||_L2: " << rhs.norm(Pimpact::ENorm::L2) << "\n";
      //if(withoutput) rhs.write(100);
    }
    else{
      auto initFunC = [](ST x, ST y) ->ST { return std::pow((y-0.5), 2); };
      auto initFunS = [](ST x, ST y) ->ST { return std::pow((x-0.5), 1); };
      auto deriFunC = [=](ST y) ->ST { return 2.*(y-0.5)/ly - iRe*2./ly/ly; };
      auto deriFunS = [=](ST x) ->ST { return 1./lx; };

      x.getCField()(Pimpact::F::U).initFromFunction(
          [=](ST x, ST y, ST z) ->ST { return initFunC(x, y); });
      x.getCField()(Pimpact::F::V).initFromFunction(
          [=](ST x, ST y, ST z) ->ST { return initFunC(x, y); });
      x.getSField()(Pimpact::F::U).initFromFunction(
          [=](ST x, ST y, ST z) ->ST { return initFunS(x, y); });
      x.getSField()(Pimpact::F::V).initFromFunction(
          [=](ST x, ST y, ST z) ->ST { return initFunS(x, y); });

      sol = x;
      //if(withoutput) x.write(10);

      // solution init
      rhs.getCField()(Pimpact::F::U).initFromFunction(
          [=](ST x, ST y, ST z) ->ST {
          if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(0)>0) ||
            ( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(0)>0) ||
            ( (y  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(1)>0) ||
            ( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(1)>0) ||
            ( (z  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(2)>0) ||
            ( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(2)>0))
            return initFunC(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.));
          else
            return a2*initFunS(x, y) + deriFunC(y); });

      rhs.getCField()(Pimpact::F::V).initFromFunction(
          [=](ST x, ST y, ST z) ->ST {
          if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(0)>0) ||
            ( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(0)>0) ||
            ( (y  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(1)>0) ||
            ( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(1)>0) ||
            ( (z  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(2)>0) ||
            ( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(2)>0))
            return initFunC(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.));
          else
            return a2*initFunS(x, y) + deriFunC(y); });

      rhs.getSField()(Pimpact::F::U).initFromFunction(
          [=](ST x, ST y, ST z) ->ST {
          if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(0)>0) ||
            ( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(0)>0) ||
            ( (y  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(1)>0) ||
            ( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(1)>0) ||
            ( (z  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(2)>0) ||
            ( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(2)>0))
            return initFunS(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.));
          else
            return -a2*initFunC(x, y) +deriFunS(x); });

      rhs.getSField()(Pimpact::F::V).initFromFunction(
          [=](ST x, ST y, ST z) ->ST {
          if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(0)>0) ||
            ( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(0)>0) ||
            ( (y  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(1)>0) ||
            ( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(1)>0) ||
            ( (z  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(2)>0) ||
            ( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(2)>0))
            return initFunS(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.));
          else
            return -a2*initFunC(x, y) + deriFunS(x); });

      //if(withoutput) rhs.write(30);

      modeOp->apply(sol, y);
      //if(withoutput) y.write(20);
      //if(2==print) y.print();
      //if(3==print) rhs.print();

      err.add(1., y, -1., rhs);
      //if(withoutput) err.write(0);
      //if(1==print) err.print();

      error = err.norm(Pimpact::ENorm::Inf)/rhs.norm(Pimpact::ENorm::Inf);
      std::cout << "\nresidual: " << error << "\n";

      x.init();
    }

    //x.random(false, Pimpact::B::N);
    //x.random(false, Pimpact::B::Y);
    //x.scale(1.e-3);
    //rhs.init(0.);

    //rhs.init(-1.);
    //x.init(0.);

    //MGFieldsT xs(mgGrids);
    ////xs.get() = rhs;
    //xs.get().init(0.);

    //modeOp->computeResidual(xs.get(), rhs, x);
    //x.write(777);

    //mgMConvDiff->getTransfers()->restriction(xs);
    //xs.get(0).write(77);
    //xs.get(1).write(99);


    modeInv->apply(rhs, x);
    //modeOp->apply(rhs, x);

    if(withoutput) x.write(10);

    if(0==realCase) {
      err.add(1., sol, -1., x);

      error = err.norm(Pimpact::ENorm::Inf);
      std::cout << "\nerror: " << error << "\n";
    }
    else {
      modeOp->apply(x, sol);
      err.add(1., rhs, -1., sol);
      if(withoutput) err.write(0);

      error = err.norm(Pimpact::ENorm::Inf);
      std::cout << "\nresidual: " << error << "\n";
    }

    if(0==grid->rankST()) {
      Teuchos::writeParameterListToXmlFile(*pl, "parameterOut.xml");
    }
    Teuchos::TimeMonitor::summarize();
  }
  MPI_Finalize();
  return 0;
}
