// Pimpact_SalarVectorSpace_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "pimpact.hpp"

#include "Pimpact_Space.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"

#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_Grid.hpp"
#include "Pimpact_ProcGrid.hpp"

#include <iostream>


namespace {

bool testMpi = true;
double errorTolSlack = 1e+1;

bool isImpactInit = false;

typedef int O;


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
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto sVS = Pimpact::createFieldSpace<O>();
  auto sIS = Pimpact::createScalarIndexSpace<O>();

  //  sVS->print();
  //  sIS->print();

  TEST_EQUALITY( sVS->nLoc_[0], sIS->eInd_[0]-sIS->sInd_[0]+1 );
  TEST_EQUALITY( sVS->nLoc_[1], sIS->eInd_[1]-sIS->sInd_[1]+1 );
  TEST_EQUALITY( sVS->nLoc_[2], sIS->eInd_[2]-sIS->sInd_[2]+1 );
}


// shows that local start/end  indexs are consisten wich nGlo
TEUCHOS_UNIT_TEST( FieldSpace, global_consistency ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto sVS = Pimpact::createFieldSpace<O>();
  auto sIS = Pimpact::createScalarIndexSpace<O>();
  auto fIIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fFIS = Pimpact::createFullFieldIndexSpaces<O>();

  //	sVS->print();
  //	sIS->print();
  const int dim = 2;

  //	for(int i=0; i<dim; ++i) {
  //	  fIIS[i]->print();
  //		fFIS[i]->print();
  //	}

  Teuchos::Tuple<O,dim> nloc;
  for(int i=0; i<dim; ++i)
    nloc[i] = sIS->eInd_[i]-sIS->sInd_[i]+1;

  std::cout << "nloc: " << nloc << "\n";
  Teuchos::Tuple<O,dim> nglo;

  MPI_Allreduce( nloc.getRawPtr(), nglo.getRawPtr(), dim,
      MPI_INT, MPI_SUM, sVS->comm_);

  for(int i=0; i<dim; ++i)
    TEST_EQUALITY( sVS->nGlo_[i], nglo[i] );
}


TEUCHOS_UNIT_TEST( ProcGrid, initialization ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto sVS = Pimpact::createFieldSpace<O>();
  auto bcg = Pimpact::createBoudaryConditionsGlobal();
  auto pgs = Pimpact::createProcGridSize<O>();

  auto pg = Pimpact::createProcGrid<O>( sVS, bcg, pgs );


}


} // end of namespace
