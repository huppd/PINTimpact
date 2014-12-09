#include <iostream>
#include <vector>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Fields.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"

#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Group.hpp"
#include "NOX_Pimpact_StatusTest.hpp"

#include "NOX.H"




namespace {


typedef double S;
typedef int O;

const int d = 3;
const int dNC = 2;

typedef Pimpact::Space<S,O,d,dNC> SpaceT;
bool testMpi = true;

S eps = 1.e0;
bool isImpactInit = false;



TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used." );
  clp.setOption(
      "eps", &eps,
      "epsilon used to check test results" );
}



TEUCHOS_UNIT_TEST( NOXPimpact_Group, SimpleNonlinear ) {

  typedef Pimpact::VectorField<SpaceT> VF;

  typedef Pimpact::MultiField<VF> MVF;

//  typedef Pimpact::ConvectionJacobianOp<SpaceT>  JOp1;
//  typedef Pimpact::HelmholtzOp<SpaceT>  JOp2;


  typedef NOX::Pimpact::Interface<MVF> Inter;
  typedef NOX::Pimpact::Vector<typename Inter::Field> NV;

  int rank = 0;

  auto pl = Teuchos::parameterList();
  pl->set( "domain", 1 );

  pl->set("nx", 257 );
  pl->set("ny", 257 );
  pl->set("lx", 2. );
  pl->set("ly", 2. );
//
//  pl->set("nx", 125 ); pl->set("ny", 125 );

//  pl->set("nx", 65 );
//  pl->set("ny", 65 );

  pl->set("Re", 1./eps );

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl, !isImpactInit );


  if( !isImpactInit ) isImpactInit=true;

  auto vel = Pimpact::create<Pimpact::VectorField>( space );

  vel->initField( Pimpact::ZeroFlow );

  auto x = Pimpact::createMultiField( *vel->clone(), 1 );
  auto f = Pimpact::createMultiField( *vel->clone(), 1 );

  auto sop = Pimpact::create<Pimpact::ConvectionDiffusionSOp>( space ) ;


  auto op =
      Pimpact::createOperatorBase(
          Pimpact::createMultiOpWrap(
              Pimpact::create<Pimpact::ConvectionVOp>(
                  Pimpact::create<Pimpact::ConvectionVWrap>(
                      sop
                  )
              )
          )
      );

  auto smoother =
      Pimpact::createOperatorBase(
          Pimpact::createMultiOpWrap(
              Pimpact::create<Pimpact::ConvectionVOp>(
                  Pimpact::create<Pimpact::ConvectionVWrap>(
                      Pimpact::create<Pimpact::ConvectionDiffusionSORSmoother>( sop )
                  )
              )
          )
      );

  auto jop = op;


  // init Fields, init and rhs
  x->getFieldPtr(0)->initField( Pimpact::RankineVortex2D );
  f->getFieldPtr(0)->initField( Pimpact::ZeroFlow );

  x->write( 97 );
  op->assignField(*x);
  op->apply( *x, *f );
  f->write( 98 );

  x->init( 0. );

  auto para = Pimpact::createLinSolverParameter("GMRES",1.e-6);

  auto lp_ = Pimpact::createLinearProblem<MVF>(
      jop, x->clone(), f->clone(), para , "GMRES" );

  lp_->setRightPrec( smoother );
//  lp_->setLeftPrec( smoother );

  auto lp = Pimpact::createInverseOperatorBase<MVF>( lp_ );

  auto inter = NOX::Pimpact::createInterface<MVF>( f, op, lp );

  auto nx = NOX::Pimpact::createVector( x );

  auto bla = Teuchos::parameterList();

  auto group = NOX::Pimpact::createGroup( bla, inter, nx );

  // Set up the status tests
  auto statusTest = NOX::Pimpact::createStatusTest();

  // Create the list of solver parameters
  Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr =
    Teuchos::rcp(new Teuchos::ParameterList);


  // Create the line search parameters sublist
  Teuchos::ParameterList& lineSearchParameters = solverParametersPtr->sublist("Line Search");
  lineSearchParameters.set("Method","Backtrack");
  lineSearchParameters.sublist("Backtrack").set("Recovery Step",1.e-12/2.);


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
  Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->write(99);

  auto sol = Teuchos::rcp_const_cast<VF>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->getConstFieldPtr(0));
  sol->write(888);

  vel->initField( Pimpact::RankineVortex2D );

  auto er = vel->clone();
  er->initField( Pimpact::ZeroFlow );
  er->add( 1., *sol, -1., *vel );
  er->write(100);

  std::cout << " error: " << er->norm() << "\n";
  TEST_EQUALITY( er->norm() < 1.e-6, true );

  Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(999);

}



} // namespace
