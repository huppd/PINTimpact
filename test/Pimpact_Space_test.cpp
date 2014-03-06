// Pimpact_SalarVectorSpace_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"

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


// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST( FieldSpace, local_consistency ) {
  // init impact
  init_impact(0,0);

  auto sVS = Pimpact::createFieldSpace<int>();
  auto sIS = Pimpact::createScalarIndexSpace<int>();

//  sVS->print();
//  sIS->print();

	TEST_EQUALITY( sVS->nLoc_[0], sIS->eInd_[0]-sIS->sInd_[0]+1 );
	TEST_EQUALITY( sVS->nLoc_[1], sIS->eInd_[1]-sIS->sInd_[1]+1 );
	TEST_EQUALITY( sVS->nLoc_[2], sIS->eInd_[2]-sIS->sInd_[2]+1 );
}


// shows that local start/end  indexs are consisten wich nGlo
TEUCHOS_UNIT_TEST( FieldSpace, global_consistency ) {
  auto sVS = Pimpact::createFieldSpace<int>();
  auto sIS = Pimpact::createScalarIndexSpace<int>();
  auto fIIS = Pimpact::createInnerFieldIndexSpaces<int>();
  auto fFIS = Pimpact::createFullFieldIndexSpaces<int>();

//	sVS->print();
//	sIS->print();
	const int dim = 2;

//	for(int i=0; i<dim; ++i) {
//	  fIIS[i]->print();
//		fFIS[i]->print();
//	}

	Teuchos::Tuple<int,dim> nloc;
	for(int i=0; i<dim; ++i)
	  nloc[i] = sIS->eInd_[i]-sIS->sInd_[i]+1;

	std::cout << "nloc: " << nloc << "\n";
	Teuchos::Tuple<int,dim> nglo;

	MPI_Allreduce( nloc.getRawPtr(), nglo.getRawPtr(), dim,
	    MPI_INT, MPI_SUM, sVS->comm_);

	for(int i=0; i<dim; ++i)
	  TEST_EQUALITY( sVS->nGlo_[i], nglo[i] );
}

} // endo of namespace

