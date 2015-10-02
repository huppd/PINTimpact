#include <ostream>
#include <fstream>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include "Teuchos_Range1D.hpp"
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"
#include "Teuchos_oblackholestream.hpp"

#include "pimpact.hpp"
#include "Pimpact_Types.hpp"

#include "Pimpact_Space.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "BelosPimpactAdapter.hpp"




int main(int argi, char** argv ) {

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; // Save some typing

  typedef double S;
  typedef int O;

  typedef Pimpact::Space<S,O,3,4> SpaceT;

  typedef Pimpact::VectorField<SpaceT> VF;
  typedef Pimpact::ScalarField<SpaceT> SF;
  typedef Pimpact::MultiField<VF>  MVF;
  typedef Pimpact::MultiField<SF>  MSF;

  // intialize MPI
  MPI_Init( &argi, &argv );

  //// get problem values form Comand line
  Teuchos::CommandLineProcessor my_CLP;

  // physical constants
  S re = 1.;
  my_CLP.setOption( "re", &re, "Reynolds number" );

  S alpha2 = 1.;
  my_CLP.setOption( "omega", &alpha2, "introduced frequency" );

  // flow type \todo + Boundary conditions
  int flow = 1;
  my_CLP.setOption( "flow", &flow, "Flow type: 0=zero flow, 1=2D Poiseuille flow in x, 2=2D Poiseuille flow in y" );

  int domain = 1;
  my_CLP.setOption( "domain", &domain, "length in x-direction" );

  // domain size
  int dim = 2;
  my_CLP.setOption( "dim", &dim, "length in x-direction" );

  S l1 = 2.;
  my_CLP.setOption( "lx", &l1, "length in x-direction" );

  S l2 = 2.;
  my_CLP.setOption( "ly", &l2, "length in y-direction" );

  S l3 = 2.;
  my_CLP.setOption( "lz", &l3, "length in z-direction" );

  // grid size
  O n1 = 33;
  my_CLP.setOption( "nx", &n1, "amount of grid points in x-direction: a*2**q+1" );

  O n2 = 33;
  my_CLP.setOption( "ny", &n2, "amount of grid points in y-direction: a*2**q+1" );

  O n3 = 2.;
  my_CLP.setOption( "nz", &n3, "amount of grid points in z-direction: a*2**q+1" );

  // processor grid size
  O np1 = 1;
  my_CLP.setOption( "npx", &np1, "amount of processors in x-direction" );

  O np2 = 1;
  my_CLP.setOption( "npy", &np2, "amount of processors in y-direction" );

  O np3 = 1.;
  my_CLP.setOption( "npz", &np3, "amount of processors in z-direction" );


  my_CLP.recogniseAllOptions(true);
  my_CLP.throwExceptions(true);

  my_CLP.parse(argi,argv);
  // end of parsing

  // starting with ininializing
  // starting with ininializing
   auto pl = Teuchos::parameterList();

   pl->set( "Re", re );
   pl->set( "alpha2", alpha2 );
   pl->set( "domain", domain );

   pl->set( "lx", l1 );
   pl->set( "ly", l2 );
   pl->set( "lz", l3 );

   pl->set( "dim", dim );

   pl->set("nx", n1 );
   pl->set("ny", n2 );
   pl->set("nz", n3 );


   // processor grid size
   pl->set("npx", np1 );
   pl->set("npy", np2 );
   pl->set("npz", np3 );

   auto space = Pimpact::createSpace<S,O,3>( pl );

   space->print();

   int rank = space->getProcGrid()->getRank();

  // outputs
  Teuchos::RCP<std::ostream> outPar;
  Teuchos::RCP<std::ostream> outLap1;
  Teuchos::RCP<std::ostream> outLap2;
  Teuchos::RCP<std::ostream> outSchur;

  if(rank==0) {
    outPar   = Teuchos::rcp( new std::ofstream("para_case.txt") );
    outLap1  = Teuchos::rcp( new std::ofstream("stats_solvLap1.txt") );
    outLap2  = Teuchos::rcp( new std::ofstream("stats_solvLap2.txt") );
    outSchur = Teuchos::rcp( new std::ofstream("stats_solvSchur.txt") );
  }	else
    outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;

  *outPar << " \tflow=" << flow << "\n";
  space->print( *outPar );


  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
  }
  outPar = Teuchos::null;


  // init vectors
  auto sca = Pimpact::createScalarField( space );
  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  auto p     = Pimpact::createMultiField<SF>(*sca,1);
  auto temps = Pimpact::createMultiField<SF>(*sca,1);
  auto u     = Pimpact::createMultiField<VF>(*vel,1);
  auto fu    = Pimpact::createMultiField<VF>(*vel,1);
  auto fp    = Pimpact::createMultiField<SF>(*sca,1);
  auto tempv = Pimpact::createMultiField<VF>(*vel,1);

  sca = Teuchos::null;
  vel = Teuchos::null;

  p->init(0.);
  u->getField(0).initField( Pimpact::EVectorField(flow) );
  u->init(0);
  u->write(1);

  tempv->getField(0).initField( Pimpact::EVectorField(flow) );
  fu->getField(0).initField( Pimpact::ZeroFlow );


  // init operators
  auto lap  =
      Pimpact::createMultiOperatorBase(
          Pimpact::create<Pimpact::HelmholtzOp>( space ) );

  auto div  = Pimpact::createMultiOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );
  auto grad = Pimpact::createMultiOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );

  // init rhs
  lap->apply( *u, *fu );
  fu->scale(-1.);
  div->apply( *u, *fp );
  fp->scale(-1.);

  fu->write( 9000 );
  fp->write( 8000 );

  u->getField(0).initField( Pimpact::ZeroFlow );
  tempv->getField(0).initField( Pimpact::ZeroFlow );

  // solve parameter for GMRES
  RCP<ParameterList> solveParaGMRES = Pimpact::createLinSolverParameter( "GMRES", 1.0e-6 );
  Teuchos::writeParameterListToXmlFile( *solveParaGMRES, "para_solverGMRES.xml" );

  // solve parameter for CG
  RCP<ParameterList> solveParaCG = Pimpact::createLinSolverParameter( "GMRES", 1.e-6 );
  Teuchos::writeParameterListToXmlFile( *solveParaCG, "para_solverCG.xml" );

  //	solveParaGMRES = Teuchos::getParametersFromXmlFile("solver1.xml");

  // create linear problems/solvers
  solveParaCG->set( "Output Stream", outLap1 );
  auto lap_problem = Pimpact::createLinearProblem<MVF>( lap, u, fu, solveParaCG, "GMRES" );

  auto lap_inv = Pimpact::createInverseOperator( lap_problem );

  auto schur = Pimpact::createOperatorBase(
      Pimpact::createTripleCompositionOp(
          div,
					lap_inv,
					grad ) );

  solveParaGMRES->set( "Output Stream", outSchur );
  auto schur_prob = Pimpact::createLinearProblem<MSF>( schur, p, temps, solveParaGMRES, "GMRES" );


  // solve stationary stokes
  lap_problem->solve( tempv, fu );
  tempv->write( 7000 );

  div->apply( *tempv, *temps );
  temps->write( 6000 );


  //	solveParaCG->set( "Output Stream", Teuchos::rcpFromRef( &std::cout ) );
  solveParaCG->set ("Verbosity",  int(Belos::Errors) );
  lap_problem->setParameters( solveParaCG );

  temps->add( -1., *fp, 1., *temps );

  schur_prob->solve( p, temps );
  p->write();

  grad->apply( *p, *tempv );
  tempv->write(2003);

  tempv->add( -1., *tempv, 1., *fu );

  solveParaCG->set ("Verbosity",  int(Belos::Errors + Belos::Warnings + Belos::IterationDetails +
      Belos::OrthoDetails + Belos::FinalSummary +	Belos::TimingDetails +
      Belos::StatusTestDetails ) );
  solveParaCG->set( "Output Stream", outLap2 );
  lap_problem->setParameters( solveParaCG );

  lap_problem->solve(u,tempv);

  u->write();

  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>( outLap1)->close();
    Teuchos::rcp_static_cast<std::ofstream>( outLap2)->close();
    Teuchos::rcp_static_cast<std::ofstream>(outSchur)->close();
  }

  MPI_Finalize();
  return( 0 );
}
