// Pimpact_SalarVectorSpace_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include "Teuchos_Range1D.hpp"
#include <Teuchos_CommHelpers.hpp>
#include "BelosOutputManager.hpp"

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_CompoundField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorMV.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "BelosPimpactAdapter.hpp"

#include "BelosMVOPTester.hpp"

#include <iostream>

namespace {

bool testMpi = true;
double errorTolSlack = 1e+1;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used." );
  clp.setOption(
      "error-tol-slack", &errorTolSlack,
      "Slack off of machine epsilon used to check test results" );
}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldScalar ) {
  init_impact(0,0);

  //	auto fS = Pimpact::createFieldSpace<int>();
  //	auto sIS = Pimpact::createScalarIndexSpace<int>();
  auto space = Pimpact::createSpace();

  auto p = Pimpact::createScalarField<double,int>( space );


  auto mvec = Pimpact::createMultiField<Pimpact::ScalarField<double,int> >(*p,5);

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<double,Pimpact::MultiField<Pimpact::ScalarField<double,int> > >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);

}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldVector ) {

  auto space = Pimpact::createSpace();

  auto vel = Pimpact::createVectorField<double,int>(space);

  auto mvec = Pimpact::createMultiField<Pimpact::VectorField<double,int> >(*vel,5);

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<double,Pimpact::MultiField<Pimpact::VectorField<double,int> > >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);

}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldScalarMode ) {

  auto space = Pimpact::createSpace();

  auto pc = Pimpact::createScalarField<double,int>(space);
  auto ps = Pimpact::createScalarField<double,int>(space);

  auto vel = Pimpact::createModeField( pc, ps );

  auto mvec = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > >( *vel, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<double> > MyOM =
      Teuchos::rcp( new Belos::OutputManager<double>(Belos::Warnings,rcp(&out,false)) );

  bool res =
      Belos::TestMultiVecTraits<double,Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > > >(MyOM,mvec);

  TEST_EQUALITY_CONST(res,true);
}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldVectorMode ) {

  auto space = Pimpact::createSpace();

  auto velc = Pimpact::createVectorField<double,int>(space);
  auto vels = Pimpact::createVectorField<double,int>(space);

  auto vel = Pimpact::createModeField( velc, vels );

  auto mvec = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> > >( *vel, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<double,Pimpact::MultiField< Pimpact::ModeField<Pimpact::VectorField<double,int> > > >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);

}



TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldCompound ) {

  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::CompoundField<VF,SF> CF;
  typedef Pimpact::MultiField<CF> MV;

  auto space = Pimpact::createSpace();

  auto vel = Pimpact::createVectorField<S,O>( space );
  auto sca = Pimpact::createScalarField<S,O>( space );

  auto q = Pimpact::createCompoundField( vel, sca );

  auto mvec = Pimpact::createMultiField< CF >( *q, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits< S, MV >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);
}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldCompoundMode ) {

  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<SF> MSF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::CompoundField<MVF,MSF> CF;
  typedef Pimpact::MultiField<CF> MV;

  auto space = Pimpact::createSpace();

  auto velc = Pimpact::createVectorField<double,int>(space);
  auto vels = Pimpact::createVectorField<double,int>(space);

  auto vel = Pimpact::createModeField( velc, vels );

  auto scac = Pimpact::createScalarField<S,O>( space );
  auto scas = Pimpact::createScalarField<S,O>( space );

  auto sca = Pimpact::createModeField( scac, scas );

  auto q = Pimpact::createCompoundField( vel, sca );

  auto mvec = Pimpact::createMultiField<CF>( *q, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<S,MV>(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);
}



TEUCHOS_UNIT_TEST( BelosOperatorMV, HelmholtzMV ) {

  typedef int O;
  typedef double S;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::MultiField<VF> MVF;

  typedef Pimpact::HelmholtzOp<S,O> Op;
  typedef Pimpact::MultiOpWrap<Op> MuOp;
  typedef Pimpact::OperatorBase<MVF> OpBase;

  auto space = Pimpact::createSpace();

  auto field = Pimpact::createVectorField<S,O>(space);

  auto mv = Pimpact::createMultiField<VF>( *field, 10 );

  auto opm = Pimpact::createMultiOpWrap<Op>();

  auto op = Pimpact::createOperatorBase<MVF,MuOp>( opm );

  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp(
      new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestOperatorTraits< S, MVF, OpBase >( MyOM, mv, op );

  TEST_EQUALITY( res, true );
}


TEUCHOS_UNIT_TEST( BelosOperatorMV, DtLapOp ) {

  typedef int O;
  typedef double S;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField<MVF> BVF;

  typedef Pimpact::DtLapOp<S,O> Op;
  typedef Pimpact::MultiOpWrap<Op> MuOp;
  typedef Pimpact::OperatorBase<BVF> OpBase;


  auto space = Pimpact::createSpace();

  auto fieldc = Pimpact::createVectorField<S,O>(space);
  auto fields = fieldc->clone();

  auto field = Pimpact::createModeField( fieldc, fields );

  auto mv = Pimpact::createMultiField<MVF>( *field, 10 );

  auto op = Pimpact::createOperatorBase<BVF,MuOp>();

  Teuchos::RCP<Belos::OutputManager<S> > myOM = Teuchos::rcp(
      new Belos::OutputManager<S>(Belos::Warnings+Belos::TimingDetails,rcp(&out,false)) );

  bool res = Belos::TestOperatorTraits< S, BVF, OpBase >( myOM, mv, op );

  TEST_EQUALITY( res, true );
}


TEUCHOS_UNIT_TEST( BelosOperatorMV, DivGrad ) {
  typedef double S;
  typedef int O;
  typedef Pimpact::ScalarField<S,O> SF;
  typedef Pimpact::MultiField<SF> BSF;

  typedef Pimpact::DivGradOp<S,O> Op;
  //  typedef Pimpact::OperatorMV<Op> BOp;
  typedef Pimpact::MultiOpWrap<Op> MuOp;
  typedef Pimpact::OperatorBase<BSF> OpBase;


	auto space = Pimpact::createSpace();

  auto temp = Pimpact::createVectorField<S,O>( space );

  auto p = Pimpact::createScalarField<double,int>( space );

  auto mv = Pimpact::createMultiField<SF>(*p,10);
  auto mv2 = mv->clone();

  auto lap = Pimpact::createOperatorBase<BSF,MuOp>(
      Pimpact::createMultiOpWrap<Op>(
          Pimpact::createDivGradOp<S,O>(temp) ) );

  Teuchos::RCP<Belos::OutputManager<S> > MyOM =
      Teuchos::rcp( new Belos::OutputManager<S>(Belos::Errors + Belos::Warnings + Belos::IterationDetails +
          Belos::OrthoDetails + Belos::FinalSummary + Belos::TimingDetails +
          Belos::StatusTestDetails + Belos::Debug, rcp(&out,false)) );


  bool res =// true;
      Belos::TestOperatorTraits< S, BSF, OpBase > (MyOM,mv,lap);

  TEST_EQUALITY( res, true );
}


TEUCHOS_UNIT_TEST( BelosOperatorMV, CompoundStokes ) {
  typedef double					                S;
  typedef int   						              O;
  typedef Pimpact::ScalarField<S,O>       SF;
  typedef Pimpact::VectorField<S,O>       VF;
  typedef Pimpact::ModeField<SF>          MSF;
  typedef Pimpact::ModeField<VF>          MVF;
  typedef Pimpact::CompoundField<MVF,MSF> CF;
  typedef Pimpact::MultiField<CF>         MV;
  typedef Pimpact::MultiOpWrap<Pimpact::CompoundStokes<S,O> > Op;
  typedef Pimpact::OperatorBase<MV> OpBase;

  auto space = Pimpact::createSpace();

  auto velc = Pimpact::createVectorField<double,int>(space);
  auto vels = Pimpact::createVectorField<double,int>(space);

  auto vel = Pimpact::createModeField( velc, vels );

  auto scac = Pimpact::createScalarField<S,O>( space );
  auto scas = Pimpact::createScalarField<S,O>( space );

  auto sca = Pimpact::createModeField( scac, scas );

  auto q = Pimpact::createCompoundField( vel, sca );

  auto mv = Pimpact::createMultiField<CF>( *q, 10 );

  auto op =
      Pimpact::createOperatorBase<MV,Op>(
          Pimpact::createMultiOpWrap(
              Pimpact::createCompoundStokes( 1., 0., 1., velc ) ) );

  Teuchos::RCP<Belos::OutputManager<S> > MyOM =
      Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res =
      Belos::TestOperatorTraits< S, MV, OpBase >( MyOM, mv, op );

  TEST_EQUALITY( res, true );
}

} // namespace

