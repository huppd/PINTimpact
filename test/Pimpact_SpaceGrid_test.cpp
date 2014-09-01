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

//#include "Pimpact_BoundaryConditionsGlobal.hpp"
//#include "Pimpact_BoundaryConditionsLocal.hpp"
//#include "Pimpact_DomainSize.hpp"
#include "Pimpact_Domain.hpp"
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

  sVS->print();
  sIS->print();

  //  for( int i=0; i<3; ++i )
  //    TEST_EQUALITY( sVS->nLoc_[i], sIS->eInd_[i]-sIS->sInd_[i]+1 );

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

  sVS->print();
  sIS->print();
  const int dim = 2;

  for(int i=0; i<dim; ++i) {
    fIIS[i]->print();
    fFIS[i]->print();
  }

  Teuchos::Tuple<O,dim> nloc;
  for(int i=0; i<dim; ++i)
    nloc[i] = sIS->eInd_[i]-sIS->sInd_[i];

  std::cout << "nloc: " << nloc << "\n";
  Teuchos::Tuple<O,dim> nglo;

//  MPI_Allreduce( nloc.getRawPtr(), nglo.getRawPtr(), dim,
//      MPI_INT, MPI_SUM, sVS->comm_);

  //  for(int i=0; i<dim; ++i)
  //    TEST_EQUALITY( sVS->nGlo_[i], nglo[i] );
}



TEUCHOS_UNIT_TEST( ProcGrid, initialization3D ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

//  auto sVS = Pimpact::createFieldSpace<O>();
  auto bcg = Pimpact::createBoudaryConditionsGlobal();
  auto gsl = Pimpact::createGridSizeLocal();
  auto pgs = Pimpact::createProcGridSize<O>(2,2,1);

  auto pg = Pimpact::createProcGrid<O>( gsl, bcg, pgs );

}



TEUCHOS_UNIT_TEST( ProcGrid, initialization4D ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }

  auto gsl = Pimpact::createGridSizeLocal<O,4>();
  auto bcg = Pimpact::createBoudaryConditionsGlobal();
  auto pgs = Pimpact::createProcGridSize<O,4>(2,1,1,2);

  auto pg = Pimpact::createProcGrid<O,4>( gsl, bcg, pgs );

}

TEUCHOS_UNIT_TEST( BoundaryConditionsGlobal, createFromImpact ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }
  auto bc = Pimpact::createBoudaryConditionsGlobal();

  bc->print();

}

TEUCHOS_UNIT_TEST( BoundaryConditionsLocal, createFromImpact ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }
  auto bc = Pimpact::createBoudaryConditionsLocal();

  bc->print();

}




TEUCHOS_UNIT_TEST( Space, create ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }
  auto space = Pimpact::createSpace();

  space->print();

}

TEUCHOS_UNIT_TEST( Space, GlobalGridCoordinates ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }
  auto space = Pimpact::createSpace();

//  auto coord =
//      Pimpact::createGridCoordinatesGlobal(space->getGridSizeGlobal(),space->getDomain()->getDomainSize());
  auto coord = space->getCoordinatesGlobal();

  coord->print();

}


TEUCHOS_UNIT_TEST( Space, LocalGridCoordinates ) {
  // init impact
  if( !isImpactInit ) {
    init_impact(0,0);
    isImpactInit=true;
  }
  auto space = Pimpact::createSpace();

//  auto coord = Pimpact::createGridCoordinatesLocal(
//      space->getFieldSpace(),
//      space->getDomain()->getDomainSize(),
//      space->getGridSizeGlobal(),
//      space->getGridSizeLocal(),
//      space->getDomain()->getBCGlobal(),
//      space->getDomain()->getBCLocal(),
//      space->getProcGrid(),
//      space->getCoordinatesGlobal() );

  auto coord = space->getCoordinatesLocal();

  if( space->getProcGrid()->getRank()==3 )
    coord->print();

}

} // end of namespace
