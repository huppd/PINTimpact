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
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_ProcGridSize.hpp"

#include "Pimpact_Space.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorMV.hpp"
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
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::MultiField<VF> MVF;
  typedef Pimpact::MultiField<SF> MSF;
  typedef Pimpact::MultiOpWrap< Pimpact::HelmholtzOp<S,O> >  Lap;
  typedef Pimpact::MultiOpWrap< Pimpact::DivHinvGrad<S,O> >  Schur;
  typedef Pimpact::MultiOpWrap< Pimpact::Grad<S,O> >  G;

  // intialize MPI
  MPI_Init( &argi, &argv );

  //// get problem values form Comand line
  Teuchos::CommandLineProcessor my_CLP;

  // physical constants
  S re = 1.;
  my_CLP.setOption( "re", &re, "Reynolds number" );

  S omega = 1.;
  my_CLP.setOption( "omega", &omega, "introduced frequency" );

  // flow type \todo + Boundary conditions
  int flow = 1;
  my_CLP.setOption( "flow", &flow, "Flow type: 0=zero flow, 1=2D Poiseuille flow in x, 2=2D Poiseuille flow in y" );

  // domain size
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
  int rank = Pimpact::init_impact_pre();

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
    //		outPar = Teuchos::rcp( &blackhole, false) ;
    outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;

  *outPar << " \tflow=" << flow << "\n";
  *outPar << " \tre=" << re << "\n";
  *outPar << " \tomega=" << omega << "\n";

  auto ds = Pimpact::createDomainSize<S>(l1,l2,l3);
  ds->set_Impact();
  ds->print( *outPar );

  auto bc = Pimpact::createBoudaryConditionsGlobal( Pimpact::Dirichelt2DChannel );
  bc->set_Impact();

  auto gs = Pimpact::createGridSizeGlobal(n1,n2,n3);
  gs->set_Impact();
  gs->print( *outPar );

  auto pgs = Pimpact::createProcGridSize<O>(np1,np2,np3);
  pgs->set_Impact();
  pgs->print( *outPar );

  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
  }
  outPar = Teuchos::null;

  Pimpact::init_impact_post();

  // init Spaces
//  auto fS = Pimpact::createFieldSpace<O>();
//
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();
  auto space = Pimpact::createSpace();

  // init vectors
  auto sca = Pimpact::createScalarField<S,O>( space );
  auto vel = Pimpact::createVectorField<S,O>( space );

  auto p     = Pimpact::createMultiField<SF>(*sca,1);
  auto temps = Pimpact::createMultiField<SF>(*sca,1);
  auto u     = Pimpact::createMultiField<VF>(*vel,1);
  auto fu     = Pimpact::createMultiField<VF>(*vel,1);
  auto fp     = Pimpact::createMultiField<SF>(*sca,1);
  auto tempv = Pimpact::createMultiField<VF>(*vel,1);

  sca = Teuchos::null;
  vel = Teuchos::null;

  p->init(0.);
  u->getField(0).initField( Pimpact::EFlowProfile(flow) );
  u->init(0);

  tempv->getField(0).initField( Pimpact::EFlowProfile(flow) );
  fu->getField(0).initField( Pimpact::ZeroProf );


  // init operators
  auto lap  = Pimpact::createMultiOperatorBase<MVF,Pimpact::HelmholtzOp<S,O> >(Pimpact::createHelmholtzOp<S,O>( 0., 1./re ) );
  auto div  = Pimpact::createMultiOpWrap<Pimpact::Div<S,O> >();
  auto grad = Pimpact::createMultiOpWrap<Pimpact::Grad<S,O> >();

  // init rhs
  lap->apply( *u, *fu );
  fu->scale(-1.);
  div->apply( *u, *fp );
  fp->scale(-1.);

  fu->write( 9000 );
  fp->write( 8000 );

  u->getField(0).initField( Pimpact::ZeroProf );
  tempv->getField(0).initField( Pimpact::ZeroProf );

  // solve parameter for GMRES
  RCP<ParameterList> solveParaGMRES = Pimpact::createLinSolverParameter( "GMRES", 1.0e-6 );

  Teuchos::writeParameterListToXmlFile( *solveParaGMRES, "para_solverGMRES.xml" );

  // solve parameter for CG
  RCP<ParameterList> solveParaCG = Pimpact::createLinSolverParameter( "GMRES", 1.e-6 );

  Teuchos::writeParameterListToXmlFile( *solveParaCG, "para_solverCG.xml" );

  //	solveParaGMRES = Teuchos::getParametersFromXmlFile("solver1.xml");

  // create problems/solvers
  solveParaCG->set( "Output Stream", outLap1 );
  auto lap_problem = Pimpact::createLinearProblem<MVF>( lap, u, fu, solveParaCG, "GMRES" );

  auto schur = Pimpact::createMultiOperatorBase<MSF, Pimpact::DivHinvGrad<S,O> >(Pimpact::createDivHinvGrad<S,O>( u, lap_problem ) );

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
