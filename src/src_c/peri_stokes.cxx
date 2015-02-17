#include <ostream>
#include <fstream>
#include <mpi.h>

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
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"

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

  typedef Pimpact::ModeField<Pimpact::VectorField<SpaceT> > VF;
  typedef Pimpact::ModeField<Pimpact::ScalarField<SpaceT> >  SF;
  typedef Pimpact::MultiField<VF> MVF;
  typedef Pimpact::MultiField<SF> MSF;


  // intialize MPI
  MPI_Init( &argi, &argv );

  //// get problem values form Comand line
  Teuchos::CommandLineProcessor my_CLP;

  // physical constants
  S re = 1.;
  my_CLP.setOption( "re", &re, "Reynolds number" );

  S alpha2 = 1.;
  my_CLP.setOption( "alpha2", &alpha2, "introduced frequency" );

  S px = 1.;
  my_CLP.setOption( "px", &px, "pressure gradient(only necessary for pulsatile flows)" );

  // flow type
  int flow = 1;
  my_CLP.setOption( "flow", &flow,
      "Flow type: 0=zero flow, 1=2D Poiseuille flow in x, 2=2D Poiseuille flow in y, 3=2D pulsatile flow in x, 4=2D pulsatile flow in, 5=2D streaming" );

  // domain type
  int domain = 1;
  my_CLP.setOption( "domain", &domain,
      "Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel" );

  // domain size
  int dim = 2.;
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
  O np1 = 2;
  my_CLP.setOption( "npx", &np1, "amount of processors in x-direction" );

  O np2 = 2;
  my_CLP.setOption( "npy", &np2, "amount of processors in y-direction" );

  O np3 = 1.;
  my_CLP.setOption( "npz", &np3, "amount of processors in z-direction" );

  // solver name
  std::string solver_name_1 = "GMRES";
  my_CLP.setOption( "solver1", &solver_name_1, "name of the solver for H" );

  std::string solver_name_2 = "GMRES";
  my_CLP.setOption( "solver2", &solver_name_2, "name of the solver for Schur complement" );

  S tol= 1.e-6;
  my_CLP.setOption( "tol", &tol, "name of the solver for Schur complement" );

  // preconditioner type
  int precType = 0;
  my_CLP.setOption( "prec", &precType, "Type of the preconditioner." );

  int precTypeSchur = 0;
  my_CLP.setOption( "precTypeSchur", &precTypeSchur, "Type of the preconditioner for the SchurComplement." );

  bool leftPrec = true;
  my_CLP.setOption( "leftPrec","rightPrec", &leftPrec, "type of fixpoint iteration matrix: 1=Newton, 2=Piccard, 3=lin diag... " );

  my_CLP.recogniseAllOptions(true);
  my_CLP.throwExceptions(true);

  my_CLP.parse(argi,argv);
  // end of parsing

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
    //		outPar = Teuchos::rcp( &blackhole, false) ;
    outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;

  *outPar << " \tflow=" << flow << "\n";
  *outPar << " \tdomain=" << domain << "\n";
  *outPar << " \tpx=" << px << "\n";

  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
  }
  outPar = Teuchos::null;


  // init vectors
  auto p     = Pimpact::createInitMSF( space );
  auto temps = p->clone();
  auto fp    = p->clone();

  auto u     = Pimpact::createInitMVF( Pimpact::EFlowType(flow), space, alpha2, px );
  auto tempv = Pimpact::createInitMVF( Pimpact::Zero2DFlow, space );
  auto fu    = Pimpact::createInitMVF( Pimpact::Zero2DFlow, space );

  u->write(3);

  p->init(0.);
  u->init(0);
  tempv->init(0);



  // init Belos operators
  auto H  =
      Pimpact::createMultiOperatorBase(
          Pimpact::createDtLapOp( space, alpha2, 1./re ) );


  // create choosen preconditioner
  Teuchos::RCP<Pimpact::OperatorBase<MVF> > lprec = Teuchos::null;

  switch(precType) {
		
  //case 1: { // didn't realy work
    //if(rank==0) std::cout << "\n\tprecType: 1, -Dt/omega\n";

    //lprec =
        //Pimpact::createMultiOperatorBase(
            //Pimpact::createDtModeOp<SpaceT>( -1./alpha2 ) );
    //break;
  //}
  case 2: {
    if(rank==0) std::cout << "\n\tprecType: 2, Lap^-1\n";

    auto solverParams = Pimpact::createLinSolverParameter( "CG", tol*l1*l2/n1/n2/1000 );
    solverParams->set ("Verbosity", int( Belos::Errors) );
    auto op =
        Pimpact::createMultiModeOperatorBase<MVF>(
            Pimpact::create<Pimpact::HelmholtzOp>( space )
        );
    // Create the Pimpact::LinearSolver solver.
    auto prob =
        Pimpact::createLinearProblem<MVF>(
            op,
            fu,
            fu, solverParams, "CG" );

    lprec = Pimpact::createOperatorBase( Pimpact::createInverseOperator<MVF>( prob ) );
    break;
  }
  case 3: {
    if(rank==0) std::cout << "\n\tprecType: 3, EddyPrec(CG)\n";

    auto solverParams = Pimpact::createLinSolverParameter( "CG", tol*l1*l2/n1/n2/1000 );
    solverParams->set ("Verbosity", int( Belos::Errors) );
    auto A =
        Pimpact::createMultiModeOperatorBase<MVF>(
            Pimpact::create<Pimpact::HelmholtzOp>( space )
        );
    auto prob2 = Pimpact::createLinearProblem<MVF>( A, fu->clone(), fu->clone(), solverParams, "CG" );
    auto op2 = Pimpact::createEddyPrec<SpaceT>( fu->clone(), Pimpact::createInverseOperatorBase<MVF>(prob2) ) ;
    lprec = Pimpact::createMultiOperatorBase( op2 );
    break;
  }
  //  case 4: {
  //    if(rank==0) std::cout << "\n\tprecType: 4, EddyPrec(ML)\n";
  //
  //    auto bla = Pimpact::createMultiModeOperatorBase<MVF>(
  //        Pimpact::createMLHelmholtzOp( space, 20, alpha2, 1./re, tol*l1*l2/n1/n2/1000 ) );
  //    auto op2 = Pimpact::createEddyPrec<SpaceT>( fu->clone(), bla ) ;
  //    lprec = Pimpact::createMultiOperatorBase<MVF>( op2 );
  //    break;
  //  }
  case 5: {
    if(rank==0) std::cout << "\n\tprecType: 5, EddyPrec(ImpactSolver)\n";

//    auto bla = Pimpact::createMultiModeOperatorBase<MVF>(
//        Pimpact::createInverseHelmholtzOp( alpha2, 1./re, tol*l1*l2/n1/n2/1000, 1000, true, true, false ) );
//    auto op2 = Pimpact::createEddyPrec<SpaceT>( fu->clone(), bla ) ;
//    lprec = Pimpact::createMultiOperatorBase<MVF>( op2 );
//    break;
  }
  case 6: {
    if(rank==0) std::cout << "\n\tprecType: 6, EddyPrec(CG+ImpactSolver)\n";

//    auto bla = Pimpact::createMultiModeOperatorBase<MVF>(
//        Pimpact::createMGVHelmholtzOp( alpha2, 1./re, false ) );
//
//    auto solverParams = Pimpact::createLinSolverParameter( "CG", tol*l1*l2/n1/n2/100 );
//    solverParams->set ("Verbosity", int( Belos::Errors) );
//
//    auto A = Pimpact::createMultiModeOperatorBase<MVF,Pimpact::HelmholtzOp >( Pimpact::createHelmholtzOp( space, alpha2, 1./re ) );
//
//    auto prob2 = Pimpact::createLinearProblem<MVF>( A, fu->clone(), fu->clone(), solverParams, "CG" );
//    if( leftPrec )
//      prob2->setLeftPrec( bla );
//    else
//      prob2->setRightPrec( bla );
//
//    auto op2 = Pimpact::createEddyPrec<SpaceT>( fu->clone(), Pimpact::createInverseOperatorBase<MVF>(prob2) ) ;
//    lprec = Pimpact::createMultiOperatorBase<MVF >( op2 );
    break;
  }
  default:
    break;
  }



  // init MV operators
  auto div  =
      Pimpact::createMultiModeOpWrap( Pimpact::create<Pimpact::DivOp>(space) );
  auto grad =
      Pimpact::createMultiModeOpWrap( Pimpact::create<Pimpact::GradOp>(space) );


  // init boundary conditions
  H->apply( *u, *fu );
  div->apply( *u, *fp );
  fu->scale(-1.);
  fp->scale(-1.);

  fu->write(6);
  fp->write(6);

  u = Pimpact::createInitMVF( Pimpact::Zero2DFlow, space );

  // create parameter for linsovlers
  auto solverParams = Pimpact::createLinSolverParameter( solver_name_1, tol*l1*l2/n1/n2 );

  solverParams->set ("Output Stream", outLap1 );
  //  if(precType==0) {
  //    solverParams->set( "Num Blocks", 100 );
  //    solverParams->set( "Maximum Iterations", 1000  );
  //  }
  //  else {
  //    solverParams->set( "Num Blocks", 10 );
  //    solverParams->set( "Maximum Iterations", 10  );
  //  }

  Teuchos::writeParameterListToXmlFile( *solverParams, "para_solver.xml" );

  // create problems/solvers
  auto H_prob = Pimpact::createLinearProblem<MVF>( H, u, fu, solverParams, solver_name_1 );
  if( leftPrec )
    H_prob->setLeftPrec( lprec );
  else
    H_prob->setRightPrec( lprec );

  auto H_inv = Pimpact::createInverseOperator( H_prob );

  auto schur =
      Pimpact::createOperatorBase(
          Pimpact::createTripleCompositionOp(
              div,
              H_inv,
              grad
          )
      );

  solverParams = Pimpact::createLinSolverParameter( solver_name_2, tol*l1*l2/n1/n2  );
  solverParams->set( "Output Stream", outSchur );
  //  solverParams->set( "Num Blocks", 10 );
  //  solverParams->set( "Maximum Iterations", 10  );

  auto schurProb = Pimpact::createLinearProblem<MSF>( schur, p,temps, solverParams, solver_name_2);

  // create choosen preconditioner
  Teuchos::RCP<Pimpact::OperatorBase<MSF> > precSchur = Teuchos::null;
  switch(precTypeSchur) {
  case 0:
    break;
  case 1: {

//    //--- inverse DivGrad
//    auto divGradPrec =
//        Pimpact::createMultiOperatorBase< MSF >(
//            Pimpact::createModeOpWrap( Pimpact::createMGVDivGradOp<S,O,3>(true) ) );
//
//    auto divGradProb =
//        Pimpact::createLinearProblem< MSF >(
//            Pimpact::createMultiOperatorBase< MSF >(
//                Pimpact::createModeOpWrap(
//                    Pimpact::createDivGradOp<S,O,3>(
//                        u->getConstFieldPtr(0)->getConstCFieldPtr()->clone(),
//                        Pimpact::createDivOp( space ),
//                        Pimpact::createGradOp( space )
//                    )
//                )
//            ),
//            Teuchos::null,
//            Teuchos::null,
//            Teuchos::parameterList(),
//            "GMRES" );
//    divGradProb->setRightPrec( divGradPrec );
//    auto divGradInv = Pimpact::createInverseOperatorBase( divGradProb );
//
//    precSchur =
//        Pimpact::createOperatorBase(
//            Pimpact::createTripleCompositionOp(
//                divGradInv,
//                Pimpact::createInverseOperator( schurProb ),
//                divGradInv ) );

    break;
  }
  default:
    break;
  }
  if( leftPrec )
    schurProb->setLeftPrec( precSchur );
  else
    schurProb->setRightPrec( precSchur );

  // solve stationary stokes
  H_prob->solve( tempv, fu );

  div->apply( *tempv, *temps );
  temps->add( -1., *fp, 1., *temps );

  solverParams = Pimpact::createLinSolverParameter( solver_name_1, tol*l1*l2/n1/n2 );
  solverParams->set( "Output Stream", outLap2 );
  solverParams->set("Verbosity", int( Belos::Errors) );
  H_prob->setParameters( solverParams );
  schurProb->solve( p, temps );
  p->write();

  grad->apply( *p, *tempv );

  tempv->add( -1., *tempv, 1., *fu );

  solverParams->set ("Verbosity",  Belos::Errors + Belos::Warnings + Belos::IterationDetails +
      Belos::OrthoDetails + Belos::FinalSummary +	Belos::TimingDetails +
      Belos::StatusTestDetails + Belos::Debug );
  solverParams->set( "Output Stream", outLap2 );
  H_prob->setParameters( solverParams );

  H_prob->solve(u,tempv);

  u->write();

  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>( outLap1)->close();
    Teuchos::rcp_static_cast<std::ofstream>( outLap2)->close();
    Teuchos::rcp_static_cast<std::ofstream>(outSchur)->close();
  }


  MPI_Finalize();
  return( 0 );

}
