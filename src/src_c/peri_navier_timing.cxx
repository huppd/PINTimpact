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
#include "Pimpact_Space.hpp"
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

#include "NOX_Pimpact_StatusTest.hpp"

#include "NOX.H"



int main(int argi, char** argv ) {

  typedef double S;
  typedef int O;
  typedef Pimpact::MultiHarmonicField< Pimpact::VectorField<S,O> > VF;
  typedef Pimpact::MultiHarmonicField< Pimpact::ScalarField<S,O> > SF;
  typedef Pimpact::CompoundField< VF, SF> CF;
  typedef Pimpact::MultiField<CF> MF;



  typedef Pimpact::MultiHarmonicOpWrap< Pimpact::Grad<S,O> > OpS2V;
  typedef Pimpact::MultiHarmonicOpWrap< Pimpact::Div<S,O> >  OpV2S;

  typedef Pimpact::MultiDtHelmholtz<S,O>  DtL;
  typedef Pimpact::MultiHarmonicNonlinear<S,O>  MAdv;
  typedef Pimpact::Add2Op<DtL,MAdv> OpV2V;
  typedef Pimpact::MultiOpWrap< Pimpact::CompoundOpWrap<OpV2V,OpS2V,OpV2S> > Op;

  typedef Pimpact::OperatorBase<MF> BOp;


  typedef NOX::Pimpact::Interface<MF> Inter;
  typedef NOX::Pimpact::Vector<typename Inter::Field> NV;

  // intialize MPI
  MPI_Init( &argi, &argv );

  //// get problem values form Comand line
  Teuchos::CommandLineProcessor my_CLP;

  // physical constants
  S re = 1.e2;

  S alpha2 = 1.e3;

  // flow type
  int flow = 7;

  // domain type
  int domain = 2;

  // domain size
  S l1 = 1.;
  S l2 = 1.;
  S l3 = 1.;

  // grid size
  O n1 = 65;
  O n2 = 65;
  O n3 = 2.;

  O nf = 4.;

  O nfs = 1.;

  O nfe = 1.;

  // processor grid size
  O np1 = 2;
  O np2 = 2;
  O np3 = 1.;

  // solver stuff
  std::string linSolName = "GMRES";
  std::string nonLinSolName = "Newton";

  std::string lineSearchName = "Backtrack";

  int maxIter = 10;

  S tolBelos = 1.e-1;
  S tolNOX   = 1.e-1;
  S tolNF    = 1.e-1;

  // end of parsing

  if( nfs==nfe ) {
    nfs=nf;
    nfe=nf+1;
  }
  else {
    nf=nfe;
  }
  // starting with ininializing
  int rank = Pimpact::init_impact_pre();

  // outputs
  Teuchos::RCP<std::ostream> outPar;
  Teuchos::RCP<std::ostream> outLinSolve;
  Teuchos::RCP<std::ostream> outPrec;
  //  Teuchos::RCP<std::ostream> outLap2;
  //  Teuchos::RCP<std::ostream> outSchur;

  if(rank==0) {
    outPar   = Teuchos::rcp( new std::ofstream("para_case.txt") );
    outLinSolve  = Teuchos::rcp( new std::ofstream("stats_linSolve.txt") );
    outPrec  = Teuchos::rcp( new std::ofstream("stats_solvPrec.txt") );
    //    outSchur = Teuchos::rcp( new std::ofstream("stats_solvSchur.txt") );
  } else
    //    outPar = Teuchos::rcp( &blackhole, false) ;
    outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;

  *outPar << " \tflow=" << flow << "\n";
  *outPar << " \tdomain=" << domain << "\n";

  auto ds = Pimpact::createDomainSize<S>( re, alpha2, l1, l2, l3 );
  ds->set_Impact();
  ds->print( *outPar );

  //  auto bc = Pimpact::createBC( Pimpact::AllPeriodic );
  auto bc = Pimpact::createBoudaryConditionsGlobal( Pimpact::EDomainType(domain) );
  bc->set_Impact();

  auto gs = Pimpact::createGridSizeGlobal( n1, n2, n3 );
  gs->set_Impact();
  gs->print( *outPar );
  *outPar << " \tnf=" << nf << "\n";

  auto pgs = Pimpact::createProcGridSize<O>( np1, np2, np3 );
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
//  auto iS = Pimpact::createScalarIndexSpace<O>();
//  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
//  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto space = Pimpact::createSpace();

  // init vectors
  auto x    = Pimpact::createMultiField( Pimpact::createCompoundField(
      Pimpact::createMultiHarmonicVectorField<S,O>( space, nfs ),
      Pimpact::createMultiHarmonicScalarField<S,O>( space, nfs )) );
  //  auto temp = x->clone();
  auto fu   = x->clone();
  //  auto force = x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone();

  // init Fields, init and rhs
  //  fu->init(1.);
  //  x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->initField( Pimpact::EFlowProfile(flow), re, alpha2/re );
  x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::EFlowProfile(flow), 1 );

  //  x->init( 0. );
  fu->init( 0. );
  //  x->random();



  /******************************************************************************************/
  //  for( nf=nfs; nf<nfe; ++nf) {
  //  for( nf=nfs; nf<nfe; nf+=2 ) {
  for( nf=nfs; nf<nfe; nf*=2) {

    if( nf!=nfs ) {
      S toltemp = x->getConstFieldPtr(0)->getConstVFieldPtr()->getConstFieldPtr(x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes()-1)->norm()/std::sqrt(l1*l2/n1/n2);
      if(0==rank) std::cout << "\n\t--- ||u_Nf||: "<<toltemp<<"\t---\n";
      if( toltemp < tolNF ) {
        if(0==rank) std::cout << "\n\t--- Nf: "<<x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes()<<"\tdof: "<<x->getLength(true)<<"\t---\n";
        break;
      }
      do {
        x->getFieldPtr(0)->getVFieldPtr()->push_back();
        x->getFieldPtr(0)->getSFieldPtr()->push_back();
        fu->getFieldPtr(0)->getVFieldPtr()->push_back();
        fu->getFieldPtr(0)->getSFieldPtr()->push_back();
      } while( x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes() < nf );
      tolNOX /= 10;
    }

    //    if(0==rank) std::cout << "\n\t--- Nf: "<<x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes()<<"\tdof: "<<x->getLength(true)<<"\t---\n";

    auto para = Pimpact::createLinSolverParameter( linSolName, tolBelos*l1*l2/n1/n2*(nfe-1)/nf, -1 );
    para->set( "Maximum Iterations", 3000 );
    para->set( "Implicit Residual Scaling", "Norm of RHS" );
    para->set( "Explicit Residual Scaling", "Norm of RHS" );


    auto opV2V =
        Pimpact::createAdd2Op(
            Pimpact::createMultiDtHelmholtz<S,O>( alpha2/re, 1./re ),
            Pimpact::createMultiHarmonicNonlinear<S,O>( /*x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone()*/ ),
            x->getConstFieldPtr(0)->getConstVFieldPtr()->clone()
        );
    auto opS2V = Pimpact::createMultiHarmonicOpWrap< Pimpact::Grad<S,O> >();
    auto opV2S = Pimpact::createMultiHarmonicOpWrap< Pimpact::Div<S,O> >();

    auto op =
        Pimpact::createMultiOperatorBase<MF>(
            Pimpact::createCompoundOpWrap(
                x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(), opV2V, opS2V, opV2S )
        );


    Teuchos::RCP<BOp> jop;
    jop = Pimpact::createMultiOperatorBase<MF>(
        Pimpact::createCompoundOpWrap(
            x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
            Pimpact::createAdd2Op(
                Pimpact::createMultiDtHelmholtz<S,O>( alpha2/re, 1./re ),
                Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
                    //                    x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(),
                    x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::DeepCopy) ),
                    x->getConstFieldPtr(0)->getConstVFieldPtr()->clone() ),
                    opS2V, opV2S )
    );


    auto lp_ = Pimpact::createLinearProblem<MF>(
        jop, x->clone(), fu->clone(), para, linSolName );
    auto lp = Pimpact::createInverseOperatorBase<MF>( lp_ );

    //    lp->setLeftPrec( lprec );

    auto inter = NOX::Pimpact::createInterface<MF>( fu, op, lp );

    auto nx = NOX::Pimpact::createVector(x);

    auto bla = Teuchos::parameterList();

    auto group = NOX::Pimpact::createGroup<Inter>( bla, inter, nx );

    // Set up the status tests
    auto statusTest = NOX::Pimpact::createStatusTest( maxIter, tolNOX, tolBelos );

    // Create the list of solver parameters
    auto solverParametersPtr =
        NOX::Pimpact::createNOXSolverParameter( nonLinSolName, lineSearchName );

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
        NOX::Solver::buildSolver( group, statusTest, solverParametersPtr);

    // Solve the nonlinear system
    //    NOX::StatusTest::StatusType status =
    solver->solve();

    // Print the parameter list
    //    if( nf==nfs ) {
    //      if(rank==0) std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
    //      if(rank==0) std::cout << "\n" << status << "\n";
    //      if(rank==0) solver->getList().print(std::cout);
    //    }

    // Get the answer
    *group = solver->getSolutionGroup();

    // Print the answer if(rank==0) std::cout << "\n" << "-- Final Solution From Solver --" << "\n";

    //    Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->write(800);
    //    Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(900);

    x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr();
    //    x = Teuchos::rcp_const_cast<NV>( group->getXPtr() )->getFieldPtr();

    //    x->write(nf*100);
  }
  /******************************************************************************************/
  x->write(800);

  MPI_Finalize();
  return( 0 );

}
