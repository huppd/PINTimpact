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
#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinSolverParameter.hpp"
#include "BelosPimpactAdapter.hpp"

#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Group.hpp"

#include "NOX_Pimpact_StatusTest.hpp"

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


  typedef Pimpact::OperatorBase<MVF> BOp;

  typedef Pimpact::MultiDtHelmholtz<S,O>  DtL;
  typedef Pimpact::MultiHarmonicNonlinear<S,O>  MAdv;
  typedef Pimpact::MultiHarmonicOpWrap< Pimpact::ForcingOp<Pimpact::VectorField<S,O> > > Fo;

  typedef Pimpact::MultiOpWrap< Pimpact::Add3Op<DtL,MAdv,Fo> > Op;

  typedef Pimpact::MultiHarmonicNonlinearJacobian<S,O>  JMAdv;


  typedef NOX::Pimpact::Interface<MVF> Inter;
  typedef NOX::Pimpact::Vector<typename Inter::Field> NV;

  // intialize MPI
  MPI_Init( &argi, &argv );

  //// get problem values form Comand line
  //  Teuchos::CommandLineProcessor my_CLP;

  // physical constants
  S re = 1.e3;
  //  my_CLP.setOption( "re", &re, "Reynolds number" );

  S alpha2 = 4.*(4*std::atan(1));
  //  my_CLP.setOption( "alpha2", &alpha2, "introduced frequency" );


  // domain type
  int domain = 4;
  //  my_CLP.setOption( "domain", &domain,
  //      "Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel" );

  // domain size
  S l1 = 1.;
  //  my_CLP.setOption( "lx", &l1, "length in x-direction" );

  S l2 = 1.;
  //  my_CLP.setOption( "ly", &l2, "length in y-direction" );

  S l3 = 1.;
  //  my_CLP.setOption( "lz", &l3, "length in z-direction" );

  int dim = 2;
  //  my_CLP.setOption( "dim", &dim, "dimension of problem" );

  // grid size
  O n1 = 33;
  //  my_CLP.setOption( "nx", &n1, "amount of grid points in x-direction: a*2**q+1" );

  O n2 = 33;
  //  my_CLP.setOption( "ny", &n2, "amount of grid points in y-direction: a*2**q+1" );

  O n3 = 2.;
  //  my_CLP.setOption( "nz", &n3, "amount of grid points in z-direction: a*2**q+1" );

  O nf = 16.;
  //  my_CLP.setOption( "nf", &nf, "amount of grid points in f-direction" );

  // processor grid size
  O np1 = 2;
  //  my_CLP.setOption( "npx", &np1, "amount of processors in x-direction" );

  O np2 = 2;
  //  my_CLP.setOption( "npy", &np2, "amount of processors in y-direction" );

  O np3 = 1.;
  //  my_CLP.setOption( "npz", &np3, "amount of processors in z-direction" );

  // solver stuff
  std::string linSolName = "GMRES";
  //  my_CLP.setOption( "linSolName", &linSolName, "name of the linear solver" );

  std::string nonLinSolName = "Newton";
  //  my_CLP.setOption( "nonLinSolName", &nonLinSolName , "name of the non linear solver" );

  std::string lineSearchName = "Backtrack";
  //  my_CLP.setOption( "linesearch", &lineSearchName, "name of the line search" );

  int maxI = 30;
  //  my_CLP.setOption( "maxI", &maxI, "maximum iterations" );

  S tol = 1.e-6;
  //  my_CLP.setOption( "tol", &tol, "tolerance for linear solver" );

  S tolNOX = 1.e-1;
  //  my_CLP.setOption( "tolNOX", &tolNOX, "tolerance for non-linear solver" );

  // end of parsing

  // starting with ininializing
  int rank = Pimpact::init_impact_pre();

  // outputs
  Teuchos::RCP<std::ostream> outPar;
  Teuchos::RCP<std::ostream> outLinSolve;
  //  Teuchos::RCP<std::ostream> outLap2;
  //  Teuchos::RCP<std::ostream> outSchur;

  if(rank==0) {
    outPar   = Teuchos::rcp( new std::ofstream("para_case.txt") );
    outLinSolve  = Teuchos::rcp( new std::ofstream("stats_linSolve.txt") );
    //    outLap2  = Teuchos::rcp( new std::ofstream("stats_solvLap2.txt") );
    //    outSchur = Teuchos::rcp( new std::ofstream("stats_solvSchur.txt") );
  } else
    //    outPar = Teuchos::rcp( &blackhole, false) ;
    outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;

  //  *outPar << " \tflow=" << flow << "\n";
  *outPar << " \tdomain=" << domain << "\n";

  auto ds = Pimpact::createDomainSize<S>( dim, re, alpha2, l1, l2, l3 );
  ds->set_Impact();
  ds->print( *outPar );

  //  auto bc = Pimpact::createBC( Pimpact::AllPeriodic );
  auto bc = Pimpact::createBoudaryConditionsGlobal( Pimpact::EDomainType(domain) );
  bc->set_Impact();

  auto gs = Pimpact::createGridSizeGlobal(n1,n2,n3);
  gs->set_Impact();
  gs->print( *outPar );
  *outPar << " \tnf=" << nf << "\n";

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
//  auto fS = Pimpact::createFieldSpace<O>();
//
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto space = Pimpact::createSpace();


  // init vectors
  //  auto p     = Pimpact::createInitMSF<S,O>( fS );
  //  auto temps = Pimpact::createInitMSF<S,O>( fS );
  //  auto fp    = Pimpact::createInitMSF<S,O>( fS );

  auto x    = Pimpact::createMultiField( Pimpact::createMultiHarmonicVectorField<S,O>( space, nf ) );
  auto temp = Pimpact::createMultiField( Pimpact::createMultiHarmonicVectorField<S,O>( space, nf ) );
  auto fu   = Pimpact::createMultiField( Pimpact::createMultiHarmonicVectorField<S,O>( space, nf ) );

  x->init(0);
  //  x->random();
  //  x->scale();
  temp->init(0);

  auto force = x->getConstFieldPtr(0)->getConst0FieldPtr()->clone();
  if( 1==dim )
    force->initField( Pimpact::BoundaryFilter1D );
  else
    force->initField( Pimpact::BoundaryFilter2D );

  //  auto op = Pimpact::createOperatorBase<MVF,Op>(
  //       Pimpact::createMultiOpWrap(
  //           Pimpact::createAdd2Op<MAdv,DtL>(
  //               Pimpact::createMultiHarmonicNonlinear<S,O>( x->getConstFieldPtr(0)->getConst0FieldPtr()->clone() ),
  //               Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
  //               temp->getFieldPtr(0)->clone() ) )
  //  );

  auto op = Pimpact::createOperatorBase<MVF,Op>(
      Pimpact::createMultiOpWrap(
          Pimpact::createAdd3Op<DtL,MAdv,Fo>(
              temp->getFieldPtr(0)->clone(),
              Pimpact::createMultiDtHelmholtz<S,O>( space, alpha2, 1./re ),
              Pimpact::createMultiHarmonicNonlinear<S,O>( /*x->getConstFieldPtr(0)->getConst0FieldPtr()->clone()*/ ),
              Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp( force ) )
          )
      )
  );


  Teuchos::RCP<BOp> jop;
  if(0==rank) std::cout << "\n\t---\titeration matrix(2): full Newton iteration\t---\n";

  typedef Pimpact::MultiOpWrap< Pimpact::Add2Op< Pimpact::Add2Op<JMAdv,DtL>, Fo > > JOp;
  jop = Pimpact::createOperatorBase<MVF,JOp>(
      Pimpact::createMultiOpWrap(
          Pimpact::createAdd2Op<Pimpact::Add2Op<JMAdv,DtL>,Fo>(
              Pimpact::createAdd2Op<JMAdv,DtL>(
                  Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
                      x->getConstFieldPtr(0)->clone() ),
                      Pimpact::createMultiDtHelmholtz<S,O>( space, alpha2, 1./re ),
                      temp->getFieldPtr(0)->clone() ) ,
                      Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp( force ) ) ,
                      temp->getFieldPtr(0)->clone() ) )
  );



  // init Fields, init and rhs
  if( 1==dim ) {
    fu->getFieldPtr(0)->getCFieldPtr(0)->initField( Pimpact::BoundaryFilter1D );
    //    fu->getFieldPtr(0)->getCFieldPtr(0)->scale( *force );
    fu->getFieldPtr(0)->getCFieldPtr(0)->scale( -1 );
    fu->getFieldPtr(0)->get0FieldPtr()->initField( Pimpact::BoundaryFilter1D );
    fu->scale( 0.5 );
  }
  else {
    fu->getFieldPtr(0)->getCFieldPtr(0)->initField( Pimpact::GaussianForcing2D );
    fu->getFieldPtr(0)->getCFieldPtr(0)->scale( *force );
    fu->getFieldPtr(0)->getCFieldPtr(0)->scale( -1 );
    fu->getFieldPtr(0)->get0FieldPtr()->initField( Pimpact::BoundaryFilter2D );
    fu->scale( 0.5 );
  }

  x->init( 0. );
  //  fu->write(100);


  auto para = Pimpact::createLinSolverParameter( linSolName, tol*l1*l2/n1/n2, -1 );
  //  auto para = Teuchos::parameterlist();
  //   para->set( "Num Blocks",          800/4  );
  para->set( "Maximum Iterations", 3000 );
  //  para->set( "Num Recycled Blocks",  20  );
  para->set( "Implicit Residual Scaling", "Norm of RHS" );
  para->set( "Explicit Residual Scaling", "Norm of RHS" );
  //   para->set( "Output Stream", outLinSolve );



  auto lp_ = Pimpact::createLinearProblem<MVF>(
      jop, x->clone(), fu->clone(), para, linSolName );

  auto lp = Pimpact::createInverseOperatorBase<MVF>( lp_ );

  auto inter = NOX::Pimpact::createInterface<MVF>( fu, op, lp );

  auto nx = NOX::Pimpact::createVector(x);

  auto bla = Teuchos::parameterList();

  auto group = NOX::Pimpact::createGroup<Inter>( bla, inter, nx );


  // Set up the status tests
  auto statusTest = NOX::Pimpact::createStatusTest( maxI, tolNOX*l1*l2/n1/n2, tol*l1*l2/n1/n2 );

  // Create the list of solver parameters
  auto solverParametersPtr =
      NOX::Pimpact::createNOXSolverParameter( nonLinSolName, lineSearchName );

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver( group, statusTest, solverParametersPtr);

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
  Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(900);

  MPI_Finalize();
  return( 0 );

}
