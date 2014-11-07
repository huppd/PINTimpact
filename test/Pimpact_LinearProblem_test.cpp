#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "BelosPimpactAdapter.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"




namespace {

typedef double  S;
typedef int     O;

typedef Pimpact::Space<S,O,3> SpaceT;

bool testMpi = true;
double eps = 1e+1;


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
}



TEUCHOS_UNIT_TEST( BelosSolver, HelmholtzMV ) {

	typedef Pimpact::MultiField< Pimpact::VectorField<SpaceT> > MF;



	auto space = Pimpact::createSpace();

	auto vel = Pimpact::createVectorField(space);

	auto x = Pimpact::createMultiField(*vel,1);
	auto b = Pimpact::createMultiField(*vel,1);

	x->init(0.);
	b->init(1.);

	auto A = Pimpact::createMultiOperatorBase< MF, Pimpact::HelmholtzOp<SpaceT> >( Pimpact::createHelmholtzOp(space) );


	auto param = Pimpact::createLinSolverParameter("GMRES",1.e-4);

	auto linprob = Pimpact::createLinearProblem<MF>(A,x,b,param,"GMRES");

	linprob->solve(x,b);
	x->write();

  TEST_EQUALITY( linprob->getSolver()->achievedTol()<1.e-4, true );

}


} // end of namespace

