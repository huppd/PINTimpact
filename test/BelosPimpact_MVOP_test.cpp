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

	TEUCHOS_STATIC_SETUP()
	{
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

	auto fS = Pimpact::createFieldSpace<int>();
//		auto sIS = Pimpact::createScalarIndexSpace<int>();
	auto p = Pimpact::createScalarField<double,int>(fS);


	auto mvec = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,5);

 // Create an output manager to handle the I/O from the solver
	Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>(Belos::Warnings,rcp(&out,false)) );

	bool res = Belos::TestMultiVecTraits<double,Pimpact::MultiField<Pimpact::ScalarField<double,int> > >(MyOM,mvec);
	TEST_EQUALITY_CONST(res,true);

}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldVector ) {

	auto fS = Pimpact::createFieldSpace<int>();

	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto mvec = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,5);

 // Create an output manager to handle the I/O from the solver
	Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>(Belos::Warnings,rcp(&out,false)) );

	bool res = Belos::TestMultiVecTraits<double,Pimpact::MultiField<Pimpact::VectorField<double,int> > >(MyOM,mvec);
	TEST_EQUALITY_CONST(res,true);

}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldScalarMode ) {
	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mvec = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>( *vel, 5 );

 // Create an output manager to handle the I/O from the solver
	Teuchos::RCP<Belos::OutputManager<double> > MyOM =
			Teuchos::rcp( new Belos::OutputManager<double>(Belos::Warnings,rcp(&out,false)) );

	bool res =
			Belos::TestMultiVecTraits<double,Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> > > >(MyOM,mvec);

	TEST_EQUALITY_CONST(res,true);
}


TEUCHOS_UNIT_TEST( AelosPimpactMV, MultiFieldVectorMode ) {

	auto fS = Pimpact::createFieldSpace<int>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto mvec = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<double,int> >,double,int>( *vel, 5 );

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

	auto fS = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto vel = Pimpact::createVectorField<S,O>( fS,iIS,fIS );
	auto sca = Pimpact::createScalarField<S,O>( fS );

	auto q = Pimpact::createCompoundField( vel, sca );

	auto mvec = Pimpact::createMultiField< CF, S, O >( *q, 5 );

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

	auto fS = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto scac = Pimpact::createScalarField<S,O>( fS );
	auto scas = Pimpact::createScalarField<S,O>( fS );

	auto sca = Pimpact::createModeField( scac, scas );

	auto q = Pimpact::createCompoundField( vel, sca );

	auto mvec = Pimpact::createMultiField<CF,S,O>( *q, 5 );

 // Create an output manager to handle the I/O from the solver
	Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

	bool res = Belos::TestMultiVecTraits<S,MV>(MyOM,mvec);
	TEST_EQUALITY_CONST(res,true);
}



TEUCHOS_UNIT_TEST( BelosOperatorMV, HelmholtzMV ) {

	auto fS = Pimpact::createFieldSpace<int>();

	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto mv = Pimpact::createMultiField<Pimpact::VectorField<double,int>,double,int>(*vel,10);

	auto helm = Pimpact::createOperatorMV<Pimpact::Helmholtz<double,int> >();

	Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>(Belos::Warnings,rcp(&out,false)) );

	bool res = Belos::TestOperatorTraits< double, Pimpact::MultiField<Pimpact::VectorField<double,int> >,Pimpact::OperatorMV<Pimpact::Helmholtz<double,int> > >(MyOM,mv,helm);

	TEST_EQUALITY( res, true );
}


TEUCHOS_UNIT_TEST( BelosOperatorMV, HelmholtzMV2 ) {

  typedef int O;
  typedef double S;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::MultiField<VF> MVF;

  typedef Pimpact::Helmholtz<S,O> Op;
  typedef Pimpact::OperatorMV<Op> OpMV;
  typedef Pimpact::OperatorBase<MVF> OpBase;
  typedef Pimpact::OperatorPimpl<MVF,Op> OpPimpl;


  auto fS = Pimpact::createFieldSpace<O>();

  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto field = Pimpact::createVectorField<S,O>(fS,iIS,fIS);

  auto mv = Pimpact::createMultiField<VF,S,O>( *field, 10 );

  auto opm = Pimpact::createOperatorMV<Op>();

  auto opPimpl = Teuchos::rcp( new OpPimpl(opm) );
  auto opBase = Teuchos::rcp_dynamic_cast<OpBase>( opPimpl );

  auto op  = Teuchos::rcp( &opBase, false );

  Teuchos::RCP<Belos::OutputManager<S> > MyOM = Teuchos::rcp(
      new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestOperatorTraits< S, MVF, Teuchos::RCP<OpBase> >( MyOM, mv, op );

  TEST_EQUALITY( res, true );
}


TEUCHOS_UNIT_TEST( BelosOperatorMV, DtL ) {

  typedef int O;
  typedef double S;
  typedef Pimpact::VectorField<S,O> VF;
  typedef Pimpact::ModeField<VF> MVF;
  typedef Pimpact::MultiField<MVF> BVF;

  typedef Pimpact::DtL<S,O> Op;
  typedef Pimpact::OperatorMV<Op> OpMV;
  typedef Pimpact::OperatorBase<BVF> OpBase;
  typedef Pimpact::OperatorPimpl<BVF,Op> OpPimpl;


  auto fS = Pimpact::createFieldSpace<O>();

  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto fieldc = Pimpact::createVectorField<S,O>(fS,iIS,fIS);
  auto fields = Pimpact::createVectorField<S,O>(fS,iIS,fIS);

  auto field = Pimpact::createModeField( fieldc, fields );

  auto mv = Pimpact::createMultiField<MVF,S,O>( *field, 10 );

  auto op = Pimpact::createOperatorBase<BVF,Op>();
  auto opp = Teuchos::rcp( &op, false );

  Teuchos::RCP<Belos::OutputManager<S> > myOM = Teuchos::rcp(
      new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestOperatorTraits< S, BVF, Teuchos::RCP<OpBase> >( myOM, mv, opp );

  TEST_EQUALITY( res, true );
}


TEUCHOS_UNIT_TEST( BelosOperatorMV, DivGrad ) {

	auto fS = Pimpact::createFieldSpace<int>();

	auto p = Pimpact::createScalarField<double,int>(fS);

	auto mv = Pimpact::createMultiField<Pimpact::ScalarField<double,int>,double,int>(*p,10);

	auto lap = Pimpact::createOperatorMV<Pimpact::Div_Grad<double,int> >();

	Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>(Belos::Warnings,rcp(&out,false)) );

	bool res =
			Belos::TestOperatorTraits<
				double,
				Pimpact::MultiField<Pimpact::ScalarField<double,int> >,
				Pimpact::OperatorMV<Pimpact::Div_Grad<double,int> >
			> (MyOM,mv,lap);

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
	typedef Pimpact::OperatorMV<Pimpact::CompoundStokes<S,O> > OP;

	auto fS = Pimpact::createFieldSpace<O>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto scac = Pimpact::createScalarField<S,O>( fS );
	auto scas = Pimpact::createScalarField<S,O>( fS );

	auto sca = Pimpact::createModeField( scac, scas );

	auto q = Pimpact::createCompoundField( vel, sca );

	auto mv = Pimpact::createMultiField<CF,S,O>( *q, 10 );

	auto op = Pimpact::createCompoundStokes( 1., 0., 1., velc );

	Teuchos::RCP<Belos::OutputManager<S> > MyOM =
			Teuchos::rcp( new Belos::OutputManager<S>(Belos::Warnings,rcp(&out,false)) );

	bool res =
			Belos::TestOperatorTraits< S, MV, OP >( MyOM, mv, op );

	TEST_EQUALITY( res, true );
}

} // namespace

