#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosMVOPTester.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_CompoundField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "BelosPimpactAdapter.hpp"



namespace {

bool testMpi = true;
double eps = 1.e-6;
int domain = 1;

typedef double S;
typedef int O;

typedef typename Pimpact::Space<S,O,3,4> SpaceT;
typedef typename Pimpact::ScalarField<SpaceT> SF;
typedef typename Pimpact::VectorField<SpaceT> VF;
typedef typename Pimpact::ModeField<SF>       MSF;
typedef typename Pimpact::ModeField<VF>       MVF;


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
}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldScalar ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl );

  auto p = Pimpact::createScalarField( space );


  auto mvec = Pimpact::createMultiField(*p,5);

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<S,Pimpact::MultiField<Pimpact::ScalarField<SpaceT> > >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);

}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldVector ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl );

  auto vel = Pimpact::create<Pimpact::VectorField>(space);

  auto mvec = Pimpact::createMultiField( *vel, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<S,Pimpact::MultiField<Pimpact::VectorField<SpaceT> > >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);

}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldScalarMode ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl );

//  auto pc = Pimpact::createScalarField(space);
//  auto ps = Pimpact::createScalarField(space);
//
//  auto vel = Pimpact::createModeField( pc, ps );
  auto vel = Pimpact::create<MSF>( space );

  auto mvec = Pimpact::createMultiField( *vel, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM =
      Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res =
      Belos::TestMultiVecTraits<S,Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<SpaceT> > > >(MyOM,mvec);

  TEST_EQUALITY_CONST(res,true);
}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldVectorMode ) {

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl );

  auto vel = Pimpact::create<MVF>( space );

  auto mvec = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<SpaceT> > >( *vel, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<S,Pimpact::MultiField< Pimpact::ModeField<Pimpact::VectorField<SpaceT> > > >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);

}



TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldCompound ) {

  typedef Pimpact::ScalarField<SpaceT> SF;
  typedef Pimpact::VectorField<SpaceT> VF;
  typedef Pimpact::CompoundField<VF,SF> CF;
  typedef Pimpact::MultiField<CF> MV;

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl );

  auto vel = Pimpact::create<Pimpact::VectorField>( space );
  auto sca = Pimpact::createScalarField( space );

  auto q = Pimpact::createCompoundField( vel, sca );

  auto mvec = Pimpact::createMultiField< CF >( *q, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits< S, MV >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);
}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldCompoundMode ) {

  typedef Pimpact::ScalarField<SpaceT> SF;
  typedef Pimpact::VectorField<SpaceT> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::CompoundField<MVF,MSF> CF;
  typedef Pimpact::MultiField<CF> MV;

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl );

  auto vel = Pimpact::create<MVF>( space );

  auto sca = Pimpact::create<MSF>( space );

  auto q = Pimpact::createCompoundField( vel, sca );

  auto mvec = Pimpact::createMultiField<CF>( *q, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<S,MV>(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);
}



TEUCHOS_UNIT_TEST( BelosOperatorMV, HelmholtzMV ) {

  typedef Pimpact::VectorField<SpaceT> VF;
  typedef Pimpact::MultiField<VF> MVF;

  typedef Pimpact::OperatorBase<MVF> OpBase;

  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

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

  //  typedef Pimpact::DtLapOp<S,O> Op;
  //  typedef Pimpact::MultiOpWrap<Op> MuOp;
  typedef Pimpact::OperatorBase<BVF> OpBase;


  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl, true );

  auto field = Pimpact::create<MVF>( space );

  auto mv = Pimpact::createMultiField<MVF>( *field, 10 );

  auto op = Pimpact::createMultiOperatorBase<BVF>( Pimpact::createDtLapOp(space) );

  Teuchos::RCP<Belos::OutputManager<S> > myOM = Teuchos::rcp(
      new Belos::OutputManager<S>(Belos::Warnings+Belos::TimingDetails,rcp(&out,false)) );

  bool res = Belos::TestOperatorTraits< S, BVF, OpBase >( myOM, mv, op );

  TEST_EQUALITY( res, true );
}


TEUCHOS_UNIT_TEST( BelosOperatorMV, DivGrad ) {

  typedef Pimpact::ScalarField<SpaceT> SF;
  typedef Pimpact::MultiField<SF> BSF;

  typedef Pimpact::OperatorBase<BSF> OpBase;


  auto pl = Teuchos::parameterList();

  pl->set( "domain", domain);

  auto space = Pimpact::createSpace( pl );

  auto temp = Pimpact::create<Pimpact::VectorField>( space );

  auto p = Pimpact::create<Pimpact::ScalarField>( space );

  auto mv = Pimpact::createMultiField(*p,10);
  auto mv2 = mv->clone();

  auto lap = Pimpact::createOperatorBase(
      Pimpact::createMultiOpWrap(
          Pimpact::createDivGradOp(
              temp,
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

