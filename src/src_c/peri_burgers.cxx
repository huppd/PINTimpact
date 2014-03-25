#include <mpi.h>

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
#include "Pimpact_GridSize.hpp"
#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinSolverParameter.hpp"
#include "BelosPimpactAdapter.hpp"

#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_LinearStokes.hpp"
#include "NOX_Pimpact_SimpleLinear.hpp"
#include "NOX_Pimpact_SimpleNonlinear.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Group.hpp"

#include "NOX.H"



int main(int argi, char** argv ) {

//  using Teuchos::ParameterList;
//  using Teuchos::parameterList;
//  using Teuchos::RCP;
//  using Teuchos::rcp; // Save some typing

  typedef double S;
  typedef int O;
  typedef Pimpact::MultiHarmonicField< Pimpact::VectorField<S,O> > VF;
  typedef Pimpact::MultiField<VF> MVF;


  typedef Pimpact::OperatorBase<MVF> BVOp;

  typedef Pimpact::MultiHarmonicNonlinear<S,O>  Op1;
  typedef Pimpact::MultiDtHelmholtz<S,O>  Op2;
  typedef Pimpact::MultiHarmonicOpWrap< Pimpact::ForcingOp<S,O> > Op3;
//  typedef Pimpact::MultiOpWrap<Pimpact::AddOp<Op1,Op2> >  Op;
  typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<Op1,Op2>, Op3 > > Op;
  typedef Pimpact::MultiHarmonicNonlinearJacobian<S,O>  JOp1;
  typedef Op2 JOp2;
//  typedef Pimpact::MultiOpWrap<Pimpact::AddOp<JOp1,JOp2> > JOp;
//  typedef Pimpact::MultiOpWrap<Pimpact::AddOp<JOp1,JOp2> > JOp;
  typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<JOp1,Op2>, Op3 > > JOp;
//  typedef Op JOp;
  typedef Pimpact::OperatorBase<MVF>  BOp;

  typedef NOX::Pimpact::Interface<MVF> Inter;
  typedef NOX::Pimpact::Vector<typename Inter::Field> NV;

  // intialize MPI
  MPI_Init( &argi, &argv );

  //// get problem values form Comand line
  Teuchos::CommandLineProcessor my_CLP;

  // physical constants
  S re = 1.e-2;
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
  S l1 = 1.;
  my_CLP.setOption( "lx", &l1, "length in x-direction" );

  S l2 = 1.;
  my_CLP.setOption( "ly", &l2, "length in y-direction" );

  S l3 = 1.;
  my_CLP.setOption( "lz", &l3, "length in z-direction" );

  // grid size
  O n1 = 33;
  my_CLP.setOption( "nx", &n1, "amount of grid points in x-direction: a*2**q+1" );

  O n2 = 7;
  my_CLP.setOption( "ny", &n2, "amount of grid points in y-direction: a*2**q+1" );

  O n3 = 2.;
  my_CLP.setOption( "nz", &n3, "amount of grid points in z-direction: a*2**q+1" );

  O nf = 4.;
  my_CLP.setOption( "nf", &nf, "amount of grid points in f-direction" );

  // processor grid size
  O np1 = 4;
  my_CLP.setOption( "npx", &np1, "amount of processors in x-direction" );

  O np2 = 1;
  my_CLP.setOption( "npy", &np2, "amount of processors in y-direction" );

  O np3 = 1.;
  my_CLP.setOption( "npz", &np3, "amount of processors in z-direction" );

  // solver name
  std::string solver_name_1 = "GMRES";
  my_CLP.setOption( "solver1", &solver_name_1, "name of the solver for H" );

  std::string solver_name_2 = "GMRES";
  my_CLP.setOption( "solver2", &solver_name_2, "name of the solver for Schur complement" );

  // preconditioner type
  int precType = 0;
  my_CLP.setOption( "prec", &precType, "Type of the preconditioner." );

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
  } else
//    outPar = Teuchos::rcp( &blackhole, false) ;
    outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;

  *outPar << " \tflow=" << flow << "\n";
  *outPar << " \tdomain=" << domain << "\n";
  *outPar << " \tre=" << re << "\n";
  *outPar << " \tpx=" << px << "\n";
  *outPar << " \talpha2=" << alpha2 << "\n";

  auto ds = Pimpact::createDomainSize<S>( l1, l2, l3 );
  ds->set_Impact();
  ds->print( *outPar );

//  auto bc = Pimpact::createBC( Pimpact::AllPeriodic );
  auto bc = Pimpact::createBC( Pimpact::EDomainType(domain) );
  bc->set_Impact();

  auto gs = Pimpact::createGridSize<O>(n1,n2,n3);
  gs->set_Impact();
  gs->print( *outPar );

  auto pgs = Pimpact::createProcGridSize<O>(np1,np2,np3);
  pgs->set_Impact();
  pgs->print( *outPar );

  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
  }
  outPar = Teuchos::null;


  // init IMPACT
  Pimpact::init_impact_post();


  // init Spaces
  auto fS = Pimpact::createFieldSpace<O>();

  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();


  // init vectors
//  auto p     = Pimpact::createInitMSF<S,O>( fS );
//  auto temps = Pimpact::createInitMSF<S,O>( fS );
//  auto fp    = Pimpact::createInitMSF<S,O>( fS );

  auto x    = Pimpact::createMultiField( Pimpact::createMultiHarmonicVectorField<S,O>( fS, iIS, fIS, nf ) );
  auto temp = Pimpact::createMultiField( Pimpact::createMultiHarmonicVectorField<S,O>( fS, iIS, fIS, nf ) );
  auto fu   = Pimpact::createMultiField( Pimpact::createMultiHarmonicVectorField<S,O>( fS, iIS, fIS, nf ) );

  x->init(0);
//  x->random();
//  x->scale();
  temp->init(0);

  auto force = x->getConstFieldPtr(0)->getConst0FieldPtr()->clone();
  force->initField( Pimpact::BoundaryFilter1D );

//  auto op = Pimpact::createOperatorBase<MVF,Op>(
//       Pimpact::createMultiOpWrap(
//           Pimpact::createAddOp<Op1,Op2>(
//               Pimpact::createMultiHarmonicNonlinear<S,O>( x->getConstFieldPtr(0)->getConst0FieldPtr()->clone() ),
//               Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
//               temp->getFieldPtr(0)->clone() ) )
//  );
  auto op = Pimpact::createOperatorBase<MVF,Op>(
       Pimpact::createMultiOpWrap(
           Pimpact::createAddOp<Pimpact::AddOp<Op1,Op2>,Op3>(
             Pimpact::createAddOp<Op1,Op2>(
               Pimpact::createMultiHarmonicNonlinear<S,O>( x->getConstFieldPtr(0)->getConst0FieldPtr()->clone() ),
               Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
               temp->getFieldPtr(0)->clone() ) ,
             Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
             temp->getFieldPtr(0)->clone() ) )
  );

  auto jop = Pimpact::createOperatorBase<MVF,JOp>(
       Pimpact::createMultiOpWrap(
           Pimpact::createAddOp<Pimpact::AddOp<JOp1,Op2>,Op3>(
             Pimpact::createAddOp<JOp1,Op2>(
               Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
                   x->getConstFieldPtr(0)->getConst0FieldPtr()->clone(), x->getConstFieldPtr(0)->clone() ),
               Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
               temp->getFieldPtr(0)->clone() ) ,
             Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
             temp->getFieldPtr(0)->clone() ) )
  );


//  auto jop = Pimpact::createOperatorBase<MVF,JOp>(
//       Pimpact::createMultiOpWrap(
//           Pimpact::createAddOp<JOp1,JOp2>(
//               Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
//                   x->getConstFieldPtr(0)->getConst0FieldPtr()->clone(), x->getConstFieldPtr(0)->clone() ),
//               Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
//               temp->getFieldPtr(0)->clone() ) ) );

  // init Fields, init and rhs
//  x ->getFieldPtr(0)->get0FieldPtr()->initField( Pimpact::RankineVortex2D );
//  x ->getFieldPtr(0)->get0FieldPtr()->initField( Pimpact::RankineVortex2D );
//  fu->getFieldPtr(0)->get0FieldPtr()->initField( Pimpact::ZeroProf );
//  fu->getFieldPtr(0)->getCFieldPtr(1)->initField( Pimpact::BoundaryFilter1D );
  fu->getFieldPtr(0)->getCFieldPtr(0)->initField( Pimpact::BoundaryFilter1D );
  fu->getFieldPtr(0)->getCFieldPtr(0)->scale( 0.5 );
  fu->getFieldPtr(0)->get0FieldPtr()->initField( Pimpact::BoundaryFilter1D );
//  fu->getFieldPtr(0)->get0FieldPtr()->initField( Pimpact::GaussianForcing2D );
//  fu->scale(re);

  x->init( 0. );
//  x->write(0);
//  op->apply(*x,*fu);
  fu->write(100);


  auto para = Pimpact::createLinSolverParameter("GMRES",1.e-5 )->get();
 //  auto para = Teuchos::parameterlist();
   para->set( "Num Blocks",          800/2  );
   para->set( "Maximum Iterations", 1600/2 );
 //  para->set( "Num Recycled Blocks",  20  );
   para->set( "Implicit Residual Scaling", "Norm of RHS");
   para->set( "Explicit Residual Scaling", "Norm of RHS" );


   auto lp = Pimpact::createLinearProblem<MVF>(
       jop, x->clone(), fu->clone(), para , "GMRES" );

   auto inter = NOX::Pimpact::createInterface<MVF>( fu, op, lp );

   auto nx = NOX::Pimpact::createVector(x);

   auto bla = Teuchos::parameterList();

   auto group = NOX::Pimpact::createGroup<Inter>( bla, inter, nx );


   // Set up the status tests
   Teuchos::RCP<NOX::StatusTest::NormF> statusTestNormF =
     Teuchos::rcp(new NOX::StatusTest::NormF( 1.0e-10 ) );
   Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestMaxIters =
     Teuchos::rcp(new NOX::StatusTest::MaxIters( 10 ) );
   Teuchos::RCP<NOX::StatusTest::NormUpdate> statusTestNormUpdate =
     Teuchos::rcp(new NOX::StatusTest::NormUpdate( 1.e-5,  NOX::Abstract::Vector::TwoNorm ) );
   Teuchos::RCP<NOX::StatusTest::Combo> statusTestsCombo =
     Teuchos::rcp(new NOX::StatusTest::Combo( NOX::StatusTest::Combo::OR,
     Teuchos::rcp(new NOX::StatusTest::Combo( NOX::StatusTest::Combo::OR,
                                             statusTestNormF,
                                             statusTestMaxIters) ),
         statusTestNormUpdate ) ) ;

   // Create the list of solver parameters
   Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr =
     Teuchos::rcp(new Teuchos::ParameterList);

 //   Select the solver (this is the default)
 //  solverParametersPtr->set("Nonlinear Solver", "Line Search Based");

   // Create the directions parameters sublist
 //  Teuchos::ParameterList&  sl = solverParametersPtr->sublist("Direction");
 //  sl.set("Method","NonlinearCG");
 //  Teuchos::ParameterList&  sll = sl.sublist("Nonlinear CG");
 //  sll.set( "Precondition", "On" );
 ////  sll.set( "Restart Frequency", 10  );

   // Create the line search parameters sublist
 //  solverParametersPtr->sublist("Line Search").set("Method","Polynomial");
 //  solverParametersPtr->sublist("Line Search").set("Method","Backtrack");
   Teuchos::ParameterList& lineSearchParameters = solverParametersPtr->sublist("Line Search");
   lineSearchParameters.set("Method","Backtrack");
   lineSearchParameters.sublist("Backtrack").set("Recovery Step",1.e-2);

   // Set the line search method
 //  lineSearchParameters.set("Method","More'-Thuente");


   // Create the solver
   Teuchos::RCP<NOX::Solver::Generic> solver =
     NOX::Solver::buildSolver( group, statusTestsCombo, solverParametersPtr);

   // Solve the nonlinear system
   NOX::StatusTest::StatusType status = solver->solve();

   // Print the parameter list
   if(rank==0) std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
   if(rank==0) std::cout << "\n" << status << "\n";
   if(rank==0) solver->getList().print(std::cout);

   // Get the answer
   *group = solver->getSolutionGroup();

   // Print the answer
   if(rank==0) std::cout << "\n" << "-- Final Solution From Solver --" << "\n";
   Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->write(800);
//   auto blabla = Teuchos::rcp_const_cast<Pimpact::VectorField<double,int> >(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->getConstFieldPtr(0));
//   vel->initField( Pimpact::RankineVortex2D );
//   auto er = vel->clone();
//   er->initField( Pimpact::ZeroProf );
//   er->add( 1., *blabla, -1., *vel );
//   er->write(100);

   Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(900);

//  auto schur =
//      Pimpact::createMultiOperatorBase<MSF,Pimpact::DivOpGrad<S,O> >(
//          Pimpact::createDivOpGrad<S,O>( u, lap_problem ) );
//
//  solverParams = Pimpact::createLinSolverParameter( solver_name_2, 1.e-7*l1*l2/n1/n2  );
//  solverParams->get()->set( "Output Stream", outSchur );
//  solverParams->get()->set( "Num Blocks", 100 );
//
//  auto schur_prob = Pimpact::createLinearProblem<MSF>( schur, p,temps, solverParams->get(), solver_name_2);
//
//  // solve stationary stokes
//  lap_problem->solve( tempv, fu );
//  //  tempv->write( 2000 );
//
//  div->apply( *tempv, *temps );
//  //temps->write( 2002 );
//  temps->add( -1., *fp, 1., *temps );
//
//  solverParams = Pimpact::createLinSolverParameter( solver_name_1, 1.e-7*l1*l2/n1/n2 );
//  solverParams->get()->set( "Output Stream", outLap2 );
//  solverParams->get()->set ("Verbosity", int( Belos::Errors) );
//  lap_problem->setParameters( solverParams->get() );
//  schur_prob->solve( p, temps );
//  p->write();
//
//  grad->apply( *p, *tempv );
//  //tempv->write(2006);
//
//  tempv->add( -1., *tempv, 1., *fu );
//
//  solverParams->get()->set ("Verbosity",  Belos::Errors + Belos::Warnings + Belos::IterationDetails +
//      Belos::OrthoDetails + Belos::FinalSummary + Belos::TimingDetails +
//      Belos::StatusTestDetails + Belos::Debug );
//  solverParams->get()->set( "Output Stream", outLap2 );
//  lap_problem->setParameters( solverParams->get() );
//
//  lap_problem->solve(u,tempv);
//
//  u->write();
//
//  if(rank==0) {
//    Teuchos::rcp_static_cast<std::ofstream>( outLap1)->close();
//    Teuchos::rcp_static_cast<std::ofstream>( outLap2)->close();
//    Teuchos::rcp_static_cast<std::ofstream>(outSchur)->close();
//  }


  MPI_Finalize();
  return( 0 );
}
