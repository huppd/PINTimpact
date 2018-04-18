/// Pimpact 
/// \author huppd
/// \date 2018


#include <fstream>
#include <ostream>

#include <mpi.h>

#include <Teuchos_Array.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include <Teuchos_Tuple.hpp>
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "Teuchos_oblackholestream.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"

#include "NOX.H"

#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Group.hpp"
#include "NOX_Pimpact_RefinementTest.hpp"
#include "NOX_Pimpact_Vector.hpp"

#include "Pimpact_Utils.hpp"
#include "Pimpact_Grid.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinSolverParameter.hpp"
#include "BelosPimpactAdapter.hpp"




#include "Pimpact_IntResCompoundOp.hpp"
#include "Pimpact_TransferCompoundOp.hpp"
#include "Pimpact_TransferTimeOp.hpp"
#include "Pimpact_TimeStokesBSmoother.hpp"

#include "Pimpact_MultiGrid.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Operator.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"

#include "Pimpact_InterpolationTimeOp.hpp"
#include "Pimpact_RestrictionTimeOp.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"


//auto CompTime = Teuchos::TimeMonitor::getNewCounter("Pimpact:: Solving Time");

using S = double;
using O = int;

using GridT = Pimpact::Grid<S, O, 3, 4, 4>;
using CGridT = Pimpact::Grid<S, O, 3, 4, 2>;

template<class GridT> using CVF =
  Pimpact::CompoundField<
  Pimpact::TimeField<Pimpact::VectorField<GridT> >,
  Pimpact::TimeField<Pimpact::ScalarField<GridT> > >;

template<class GridT> using INT = Pimpact::IntResCompoundOp<
                                   Pimpact::InterpolationTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<GridT> > >,
                                   Pimpact::InterpolationTimeOp<                          Pimpact::InterpolationOp<GridT> > >;

template<class GridT> using RES = Pimpact::IntResCompoundOp<
                                   Pimpact::RestrictionTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<GridT> > >,
                                   Pimpact::RestrictionTimeOp<                          Pimpact::RestrictionSFOp<GridT> > >;

template<class GridT1, class GridT2> using TCO = Pimpact::TransferCompoundOp<
    Pimpact::TransferTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<GridT1, GridT2> > >,
    Pimpact::TransferTimeOp<                          Pimpact::TransferOp<GridT1, GridT2> > >;

template<class T> using MOP = Pimpact::InverseOp<T>;



using CS4L = Pimpact::CoarsenStrategy<GridT, CGridT>;
using CS4G = Pimpact::CoarsenStrategyGlobal<GridT, CGridT, 5>;

using VF = Pimpact::TimeField<Pimpact::VectorField<GridT> >;
using SF = Pimpact::TimeField<Pimpact::ScalarField<GridT> >;
using CF = Pimpact::CompoundField<VF, SF>;

using MF = Pimpact::MultiField<CF>;

using Fo = Pimpact::ForcingOp<VF>;
using BOp = Pimpact::OperatorBase<MF>;

using Inter = NOX::Pimpact::Interface<MF>;
using NV = NOX::Pimpact::Vector<typename Inter::FieldT>;




int main(int argi, char** argv) {

  // intialize MPI
  MPI_Init(&argi, &argv);

  //// get problem values form Comand line
  Teuchos::CommandLineProcessor my_CLP;

  // physical constants
  S re = 1.e0;
  my_CLP.setOption("re", &re, "Reynolds number");

  S alpha2 = 1.;
  my_CLP.setOption("alpha2", &alpha2, "introduced frequency");

  S rad = 0.1;
  my_CLP.setOption("radius", &rad, "radius of disk or sig of dipol");

  S amp = 0.1;
  my_CLP.setOption("amp", &amp, "amplitude of elongation of disc");

  S xm = 1./6.;
  my_CLP.setOption("xm", &xm, "rotation of disc");

  S ym = 0.5;
  my_CLP.setOption("ym", &ym, "rotation of disc");

  // flow type
  int flow = 3;
  my_CLP.setOption("flow", &flow,
                    "Flow type: 0=zero flow, 1=2D Poiseuille flow in x, 2=2D Poiseuille flow in y, 3=2D pulsatile flow in x, 4=2D pulsatile flow in, 5=2D streaming");

  int forcing = 0;
  my_CLP.setOption("force", &forcing,
                    "forcing, ja?");


  // domain type
  int domain = 0;
  my_CLP.setOption("domain", &domain,
                    "Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel");

  // domain size
  S l1 = 2.;
  my_CLP.setOption("lx", &l1, "length in x-direction");

  S l2 = 2.;
  my_CLP.setOption("ly", &l2, "length in y-direction");

  S l3 = 2.;
  my_CLP.setOption("lz", &l3, "length in z-direction");

  int dim = 3;
  my_CLP.setOption("dim", &dim, "dimension of problem");

  // grid size
  O n1 = 17;
  my_CLP.setOption("nx", &n1, "amount of grid points in x-direction: a*2**q+1");

  O n2 = 17;
  my_CLP.setOption("ny", &n2, "amount of grid points in y-direction: a*2**q+1");

  O n3 = 17;
  my_CLP.setOption("nz", &n3, "amount of grid points in z-direction: a*2**q+1");

  O nt = 16;
  my_CLP.setOption("nt", &nt, "amount of grid points in time-direction");

  const bool cny = true;


  // processor grid size
  O np1 = 1;
  my_CLP.setOption("npx", &np1, "amount of processors in x-direction");

  O np2 = 1;
  my_CLP.setOption("npy", &np2, "amount of processors in y-direction");

  O np3 = 1.;
  my_CLP.setOption("npz", &np3, "amount of processors in z-direction");

  O np4 = 1.;
  my_CLP.setOption("npt", &np4, "amount of processors in z-direction");

  // solver stuff
  std::string linSolName = "GMRES";
  my_CLP.setOption("linSolName", &linSolName, "name of the linear solver");

  std::string nonLinSolName = "Newton";
  my_CLP.setOption("nonLinSolName", &nonLinSolName , "name of the non linear solver");

  std::string lineSearchName = "Backtrack";
  my_CLP.setOption("linesearch", &lineSearchName, "name of the line search");

  int fixType = 1.;
  my_CLP.setOption("fixType", &fixType, "type of fixpoint iteration matrix: 1=Newton, 2=Piccard, 3=lin diag... ");

  bool isNewton = true;
  my_CLP.setOption("newton", "piccard", &isNewton, "type of fixpoint iteration matrix: 1=Newton, 2=Piccard, 3=lin diag... ");

  bool leftPrec = true;
  my_CLP.setOption("leftPrec", "rightPrec", &leftPrec, "type of fixpoint iteration matrix: 1=Newton, 2=Piccard, 3=lin diag... ");

  int precType = 0.;
  my_CLP.setOption("precType", &precType, "type of preconditioners ");

  int maxIter = 10;
  my_CLP.setOption("maxIter", &maxIter, "maximum iterations");

  S tolBelos = 1.e-4;
  my_CLP.setOption("tolBelos", &tolBelos, "tolerance for linear solver");

  S tolNOX = 1.e-2;
  my_CLP.setOption("tolNOX", &tolNOX, "tolerance for non-linear solver");

  S tolSchur = 1.e-4;
  my_CLP.setOption("tolSchur", &tolSchur, "tolerance for non-linear solver");

  S tolPrec = 1.e-4;
  my_CLP.setOption("tolPrec", &tolPrec, "tolerance for non-linear solver");


  my_CLP.recogniseAllOptions(true);
  my_CLP.throwExceptions(true);

  my_CLP.parse(argi, argv);
  // end of parsing

  // starting with ininializing
  auto pl = Teuchos::parameterList();

  pl->set("Re", re);
  pl->set("alpha2", alpha2);
  //pl->set("domain", domain);

  pl->set("lx", l1);
  pl->set("ly", l2);
  pl->set("lz", l3);

  //pl->set("dim", dim);

  pl->set("nx", n1);
  pl->set("ny", n2);
  pl->set("nz", n3);

  pl->set("nf", nt);

  // processor grid size
  pl->set("npx", np1);
  pl->set("npy", np2);
  pl->set("npz", np3);
  pl->set("npf", np4);

  auto grid = Pimpact::create<Pimpact::Grid<S, O, 3, 4, 4> >(pl);
  //  grid->print();
  int rank = grid->getProcGrid()->getRank();

  // outputs
  Teuchos::RCP<std::ostream> outPar;
  Teuchos::RCP<std::ostream> outLinSolve;
  Teuchos::RCP<std::ostream> outPrec;
  Teuchos::RCP<std::ostream> outSchur;

  if(rank==0) {
    outPar   = Teuchos::rcp(new std::ofstream("para_case.txt"));
    outLinSolve  = Teuchos::rcp(new std::ofstream("stats_linSolve.txt"));
    outPrec  = Teuchos::rcp(new std::ofstream("stats_solvPrec.txt"));
    outSchur  = Teuchos::rcp(new std::ofstream("stats_solvSchur.txt"));
  } else
    outPar = Teuchos::rcp(new Teuchos::oblackholestream()) ;


  *outPar << " \tCrankNicolson=" << cny << "\n";
  *outPar << " \tNewton=" << isNewton << "\n";
  *outPar << " \tflow=" << flow << "\n";
  *outPar << " \tforce=" << forcing << "\n";
  *outPar << " \tdomain=" << domain << "\n";

  // init grid


  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
  }
  outPar = Teuchos::null;



  // init vectors
  auto x = Pimpact::wrapMultiField(Pimpact::create<CF>(grid));


  // init Fields, init and rhs
  Pimpact::initVectorTimeField(x->getField(0).getVField(), Pimpact::EFlowType(flow));
  //x->init(0.);
  x->random();

  auto fu   = x->clone(Pimpact::ECopy::Shallow);
  fu->init(0.);

  Teuchos::RCP<VF> force=Teuchos::null;
  Teuchos::RCP<VF> forcem1=Teuchos::null;

  if(0!=forcing) {

    force = x->getField(0).getVField().clone(Pimpact::ECopy::Shallow);

    Pimpact::initVectorTimeField(
      *force,
      Pimpact::OscilatingDisc2D,
      xm*l1, ym*l2, rad, amp);
    //		force->write(500);
    //
    Pimpact::initVectorTimeField(
      fu->getField(0).getVField(),
      Pimpact::OscilatingDisc2DVel,
      xm*l1, ym*l2, rad, amp*alpha2/re);

    fu->getField(0).getVField().scale(*force);

    forcem1 = force->clone(Pimpact::ECopy::Shallow);
    forcem1->init(1.);
    forcem1->add(1., *forcem1, -1., *force);

    x->getField(0).getVField().scale(*force);

  }



  auto para = Pimpact::createLinSolverParameter(linSolName, tolBelos, 10);
  //  auto para = Pimpact::createLinSolverParameter(linSolName, tol, -1);
  //  para->set("Maximum Iterations", 30000);
  para->set("Output Stream", outLinSolve);


  Teuchos::RCP<Fo> forcingOp = Teuchos::null;
  Teuchos::RCP<Fo> forcingm1Op = Teuchos::null;
  if(0!=forcing) {
    forcingOp   = Pimpact::createForcingOp(force  );
    forcingm1Op = Pimpact::createForcingOp(forcem1);
  }

  auto op = Pimpact::create<Pimpact::TimeNSOp<GridT> >(grid);


  Teuchos::RCP<BOp> jop;
  {

    auto mgGrids = Pimpact::createMGGrids<CS4L>(grid, 10);

    auto mgPL = Teuchos::parameterList();
    mgPL->sublist("Smoother").set<std::string>("Solver name", "GMRES");
    mgPL->sublist("Smoother").set("Solver",
                                   *Pimpact::createLinSolverParameter("GMRES", 1.e-16, -1,
                                       Teuchos::rcp<std::ostream>(new Teuchos::oblackholestream()), 4));

    mgPL->sublist("Smoother").set<int>("numIters", 4);
    mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<int>("Maximum Iterations", 1000);
    mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "GMRES");
    mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver");
    mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<S>("Convergence Tolerance" , 1.e-1);

    auto mg = Pimpact::createMultiGrid<
      CVF,
      TCO,
      RES,
      INT,
      Pimpact::TimeNSOp,
      Pimpact::TimeNSOp,
      Pimpact::TimeNS4DBSmoother,
      //									Pimpact::TimeStokesBSmoother
      MOP
        > (mgGrids, op, mgPL);


    jop = Pimpact::createMultiOperatorBase(mg) ;

  }

  auto inter = NOX::Pimpact::createInterface<MF>(fu, Pimpact::createMultiOperatorBase(op), jop);

  auto nx = NOX::Pimpact::createVector(x);

  auto bla = Teuchos::parameterList();

  auto group = NOX::Pimpact::createGroup<Inter>(bla, inter, nx);

  // Set up the status tests
  Teuchos::RCP<NOX::StatusTest::Generic> statusTest =
    NOX::StatusTest::buildStatusTests(pl->sublist("NOX Status Test"), NOX::Utils());
  //auto statusTest = NOX::Pimpact::createStatusTest(maxIter, tolNOX, tolBelos*1e-4);

  // Create the solver
  Teuchos::RCP<Teuchos::ParameterList> noxSolverPara =
    Teuchos::sublist(pl, "NOX Solver");

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(group, statusTest, noxSolverPara);

  if(0==rank) std::cout << "\n\t--- Nf: 0\tdof: "<< x->getLength()<< "\t---\n";
// Solve the nonlinear system
  {
    Teuchos::TimeMonitor LocalTimer(*Teuchos::TimeMonitor::getNewCounter("Pimpact:: Solving Time"));
    //    NOX::StatusTest::StatusType status =
    solver->solve();
  }

//    // Print the parameter list
//    //    if(nt==nfs) {
//    //      if(rank==0) std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
//    //      if(rank==0) std::cout << "\n" << status << "\n";
//    //      if(rank==0) solver->getList().print(std::cout);
//    //    }

// Get the answer
  *group = solver->getSolutionGroup();

// Get a summary from the time monitor.
  Teuchos::TimeMonitor::summarize();

//    // Print the answer if(rank==0) std::cout << "\n" << "-- Final Solution From Solver --" << "\n";
//
  Teuchos::rcp_dynamic_cast<const NV>(group->getXPtr())->getConstFieldPtr()->write();
  Teuchos::rcp_dynamic_cast<const NV>(group->getFPtr())->getConstFieldPtr()->write(1000);

  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>(outPrec)->close();
    Teuchos::rcp_static_cast<std::ofstream>(outSchur)->close();
  }

  MPI_Finalize();
  return 0;

}
