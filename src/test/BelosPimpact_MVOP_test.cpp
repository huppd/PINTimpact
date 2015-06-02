#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosMVOPTester.hpp"

#include "Pimpact_Fields.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "BelosPimpactAdapter.hpp"



namespace {


typedef double S;
typedef int O;
const int d = 3;
const int dNC = 4;

bool testMpi = true;
double eps = 1.e-6;
int domain = 0;


typedef typename Pimpact::Space<S,O,d,dNC> SpaceT;

typedef typename Pimpact::ScalarField<SpaceT> SF;
typedef typename Pimpact::VectorField<SpaceT> VF;
typedef typename Pimpact::ModeField<SF>       MSF;
typedef typename Pimpact::ModeField<VF>       MVF;

typedef Pimpact::CompoundField<VF,SF> CF;
typedef Pimpact::CompoundField<MVF,MSF> CMF;


auto pl = Teuchos::parameterList();


TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used." );
  clp.setOption(
      "error-tol-slack", &eps,
      "Slack off of machine epsilon used to check test results" );
  clp.setOption(
      "domain", &domain,
      "domain" );

  pl->set( "domain", domain );
  pl->set( "dim", 3 );
  pl->set( "nx", 24 );
  pl->set( "ny", 17 );
  pl->set( "nz", 9 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, BelosMVTest, FType ) {

  auto space = Pimpact::createSpace<S,O,d,dNC>( pl );

  auto vel = Pimpact::create<FType>( space );

  auto mvec = Pimpact::createMultiField( *vel, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM =
		Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<S,Pimpact::MultiField<FType> >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, BelosMVTest, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, BelosMVTest, VF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, BelosMVTest, MSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, BelosMVTest, MVF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, BelosMVTest, CF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, BelosMVTest, CMF )




TEUCHOS_UNIT_TEST( BelosOperatorMV, HelmholtzMV ) {

  typedef Pimpact::VectorField<SpaceT> VF;
  typedef Pimpact::MultiField<VF> MVF;

  typedef Pimpact::OperatorBase<MVF> OpBase;


  auto space = Pimpact::createSpace( pl );

  auto field = Pimpact::create<Pimpact::VectorField>(space);

  auto mv = Pimpact::createMultiField<VF>( *field, 10 );

  auto opm = Pimpact::createMultiOpWrap( Pimpact::create<Pimpact::HelmholtzOp>(space) );

  auto op = Pimpact::createOperatorBase( opm );

  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp(
      new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestOperatorTraits< S, MVF, OpBase >( MyOM, mv, op );

  TEST_EQUALITY( res, true );
}



TEUCHOS_UNIT_TEST( BelosOperatorMV, DtLapOp ) {

  typedef Pimpact::MultiField<MVF> BVF;

  typedef Pimpact::OperatorBase<BVF> OpBase;

  auto space = Pimpact::createSpace( pl );

  auto field = Pimpact::create<MVF>( space );

  auto mv = Pimpact::createMultiField<MVF>( *field, 10 );

  auto op = Pimpact::createMultiOperatorBase( Pimpact::createDtLapOp(space) );

  Teuchos::RCP<Belos::OutputManager<S> > myOM = Teuchos::rcp(
      new Belos::OutputManager<S>(Belos::Warnings+Belos::TimingDetails,rcp(&out,false)) );

  bool res = Belos::TestOperatorTraits< S, BVF, OpBase >( myOM, mv, op );

  TEST_EQUALITY( res, true );
}



TEUCHOS_UNIT_TEST( BelosOperatorMV, DivGrad ) {

  typedef Pimpact::ScalarField<SpaceT> SF;
  typedef Pimpact::MultiField<SF> BSF;

  typedef Pimpact::OperatorBase<BSF> OpBase;


  auto space = Pimpact::createSpace( pl );

  auto temp = Pimpact::create<Pimpact::VectorField>( space );

  auto p = Pimpact::create<Pimpact::ScalarField>( space );

  auto mv = Pimpact::createMultiField(*p,10);
  auto mv2 = mv->clone();

  auto lap = Pimpact::createOperatorBase(
      Pimpact::createMultiOpWrap(
          Pimpact::createDivGradOp(
              Pimpact::create<Pimpact::DivOp>( space ),
              Pimpact::create<Pimpact::GradOp>(space )
          )
      )
  );

  Teuchos::RCP<Belos::OutputManager<S> > MyOM =
      Teuchos::rcp( new Belos::OutputManager<S>(Belos::Errors + Belos::Warnings + Belos::IterationDetails +
          Belos::OrthoDetails + Belos::FinalSummary + Belos::TimingDetails +
          Belos::StatusTestDetails + Belos::Debug, rcp(&out,false)) );


  bool res =// true;
      Belos::TestOperatorTraits< S, BSF, OpBase > (MyOM,mv,lap);

  TEST_EQUALITY( res, true );

}



} // namespace

