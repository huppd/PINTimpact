#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_VectorFieldOpWrap.hpp"
#include "Pimpact_MultiGrid.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Operator.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"

#include "Pimpact_InterpolationTimeOp.hpp"
#include "Pimpact_RestrictionTimeOp.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"


#include "Pimpact_IntResCompoundOp.hpp"
#include "Pimpact_TransferCompoundOp.hpp"
#include "Pimpact_TransferTimeOp.hpp"
#include "Pimpact_TimeStokesBSmoother.hpp"

namespace {

typedef double S;
typedef int O;

typedef Pimpact::Space<S,O,3,4> FSpace3T;
typedef Pimpact::Space<S,O,4,4> FSpace4T;

typedef Pimpact::Space<S,O,3,2> CSpace3T;
typedef Pimpact::Space<S,O,4,2> CSpace4T;

typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS3L;
typedef Pimpact::CoarsenStrategyGlobal<FSpace3T,CSpace3T,5> CS3G;

typedef Pimpact::CoarsenStrategy<FSpace4T,CSpace4T> CS4L;
typedef Pimpact::CoarsenStrategyGlobal<FSpace4T,CSpace4T,5> CS4G;

template<class ST> using BSF = Pimpact::MultiField< Pimpact::ScalarField<ST> >;
//template<class T> using BVF = Pimpact::MultiField< Pimpact::VectorField<T> >;

template<class ST> using TSF = Pimpact::TimeField< Pimpact::ScalarField<ST> >;
template<class ST> using TOPF = Pimpact::TimeOpWrap< Pimpact::DivGradOp<ST> >;
template<class ST> using TOPC = Pimpact::TimeOpWrap< Pimpact::DivGradO2Op<ST> >;


template<class ST> using BOPF = Pimpact::MultiOpWrap< Pimpact::DivGradOp<ST> >;
template<class ST> using BOPC = Pimpact::MultiOpWrap< Pimpact::DivGradO2Op<ST> >;
template<class ST> using BSM = Pimpact::MultiOpWrap< Pimpact::DivGradO2JSmoother<ST> >;

template<class T> using ConvDiffOpT = Pimpact::ConvectionVOp<Pimpact::ConvectionDiffusionSOp<T> >;

//template<class T> using ConvDiffOpT = Pimpact::ConvectionVOp<Pimpact::ConvectionDiffusionSOp<T> >;

template<class T> using ConvDiffSORT = Pimpact::ConvectionVSmoother<T,Pimpact::ConvectionDiffusionSORSmoother >;

template<class T> using ConvDiffJT = Pimpact::ConvectionVSmoother<T,Pimpact::ConvectionDiffusionJSmoother >;

bool testMpi = true;
double eps = 1e-6;

int domain = 0;
int ftype = 0;

int fs = 0;
int fe = 4;

int npx = 1;
int npy = 1;
int npz = 1;
int npf = 1;

int nx = 32;
int ny = 32;
int nz = 32;
int nf = 32;

int rankbla = -1;

int maxGrids = 10;

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
	    "Slack off of machine epsilon used to check test results" );
	clp.setOption(
	    "ftype", &ftype,
	    "Slack off of machine epsilon used to check test results" );
	clp.setOption(
	    "fs", &fs,
	    "Slack off of machine epsilon used to check test results" );
	clp.setOption(
	    "fe", &fe,
	    "Slack off of machine epsilon used to check test results" );
	clp.setOption(
	    "npx", &npx,
	    "" );
	clp.setOption(
	    "npy", &npy,
	    "" );
	clp.setOption(
	    "npz", &npz,
	    "" );
	clp.setOption(
	    "npf", &npf,
	    "" );
	clp.setOption(
	    "nx", &nx,
	    "" );
	clp.setOption(
	    "ny", &ny,
	    "" );
	clp.setOption(
	    "nz", &nz,
	    "" );
	clp.setOption(
	    "nf", &nf,
	    "" );
	clp.setOption(
	    "rank", &rankbla,
	    "" );
	clp.setOption(
	    "maxGrids", &maxGrids,
	    "" );

  pl->set( "dim", 3 );
  pl->set( "domain", domain );

  pl->set( "lx", 8. );
  pl->set( "ly", 2. );
  pl->set( "lz", 4. );

  pl->set<S>( "Re", 1 );

}

/*
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSpaces, constructor4D, CS ) {

//  typedef Pimpact::CoarsenStrategy<FSpace4T,CSpace4T> CS;

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );


  auto space = Pimpact::createSpace<S,O,4>( pl );

  space->print();

  auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, maxGrids );
  std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";

	if( 0==space->rankST() )
		mgSpaces->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSpaces, constructor4D, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSpaces, constructor4D, CS4G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGFields, SF_constructor4D, CS ) {

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,4>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, maxGrids );
  std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
  if( space->rankST()==0 )
    mgSpaces->print();

  auto mgFields = Pimpact::createMGFields<TSF>( mgSpaces );

  auto field = mgFields->get( -1 );
	if( mgSpaces->participating(-1) ){

    field->init(5.);
    TEST_FLOATING_EQUALITY( std::sqrt( std::pow(5.,2)*field->getLength() ), field->norm(Belos::TwoNorm), eps );

		field->write(0);
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, SF_constructor4D, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, SF_constructor4D, CS4G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGFields, VF_constructor4D, CS ) {

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,4>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, maxGrids );
  std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
  if( space->rankST()==0 )
    mgSpaces->print();

  auto mgFields = Pimpact::createMGFields<TSF>( mgSpaces );

  auto field = mgFields->get( -1 );
	if( mgSpaces->participating(-1) ){

    field->init(5.);
    TEST_FLOATING_EQUALITY( std::sqrt( std::pow(5.,2)*field->getLength() ), field->norm(Belos::TwoNorm), eps );

		field->write(0);
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, VF_constructor4D, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, VF_constructor4D, CS4G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGOperators, SF_constructor4D, CS ) {

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,4>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, maxGrids );

  //auto OpsT = Pimpact::createTimeOpWrap< Pimpact::DivGradOp<SpaceT>, true >(
    //  Pimpact::create<Pimpact::DivGradOp>( space ) );

  //auto Ops2T = Pimpact::createTimeOpWrap< Pimpact::DivGradO2Op<SpaceT>, true >(
    //  Pimpact::create<Pimpact::DivGradO2Op>( space ) );

  auto mgOps2 = Pimpact::createMGOperators<TOPF,TOPC>( mgSpaces );

  auto op = mgOps2->get( -1 );

	if( mgSpaces->participating(-1) )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, SF_constructor4D, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, SF_constructor4D, CS4G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGOperators, VF_constructor4D, CS ) {

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,4>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, maxGrids );

  auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>( mgSpaces );


  auto op = mgOps->get( -1 );

	if( mgSpaces->participating(-1) )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, VF_constructor4D, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, VF_constructor4D, CS4G )
*/

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, Restrictor4D, CS ) {

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

	auto space = Pimpact::createSpace<S,O,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, maxGrids );

 	auto op = Pimpact::createRestrictionTimeOp<Pimpact::RestrictionOp<CSpace4T> >(mgSpaces->get(-2),mgSpaces->get(-1));

       	auto op1 = Pimpact::createRestrictionTimeOp<Pimpact::RestrictionOp<CSpace4T> >(mgSpaces->get(-3),mgSpaces->get(-2));

	auto fieldff = Pimpact::createTimeField<Pimpact::ScalarField<CSpace4T>,CSpace4T>( mgSpaces->get( -3 ));
	auto fieldf = Pimpact::createTimeField<Pimpact::ScalarField<CSpace4T>,CSpace4T>( mgSpaces->get( -2 ));
	auto fieldc = Pimpact::createTimeField<Pimpact::ScalarField<CSpace4T>,CSpace4T>( mgSpaces->get( -1 ));
		
	auto sol = fieldc->clone();
	auto er = fieldc->clone();

		// the zero test
		fieldff->init( 0. );
		fieldf->init( 1. );
		fieldc->init( 1. );

		if( mgSpaces->participating(-3) )
			op1->apply( *fieldff, *fieldf );
		if( mgSpaces->participating(-2) )
			op->apply( *fieldf, *fieldc );

		if( mgSpaces->participating(-2) )
			TEST_FLOATING_EQUALITY( 0., fieldf->norm(), eps );

		if( mgSpaces->participating(-1) )
			TEST_FLOATING_EQUALITY( 0., fieldc->norm(), eps );
                
		//fieldff->write(100);
                //fieldf->write(200);
                //fieldc->write(300);



		// the random test
		fieldff->random();

		if( mgSpaces->participating(-3) )
			op1->apply( *fieldff, *fieldf );
		if( mgSpaces->participating(-2) )
			op->apply( *fieldf, *fieldc );
		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );


		// the strong test
		fieldff->init(1.);
		fieldf->init(0.);
		fieldc->init(0.);
		sol->init(1.);
		er->random();

		if( mgSpaces->participating(-3) )
			op1->apply( *fieldff, *fieldf );
		if( mgSpaces->participating(-2) )
			op->apply( *fieldf, *fieldc );

		if( mgSpaces->participating(-1) ) {

			er->add( 1., *sol, -1., *fieldc );
			double bla = er->norm(Belos::InfNorm);
			if( 0==space->rankST() )
				std::cout << "error Const: " << bla << "\n";
			TEST_EQUALITY( er->norm(Belos::InfNorm)<eps, true ); // boundaries?
			if( bla>=eps )
				er->write(0);
		}
//		fieldff->write(100);
//		fieldf->write(200);
//		fieldc->write(300);

/*
		// the hard test in X
		fieldff->initField( Pimpact::Grad2D_inX, 1. );
		fieldf->initField( Pimpact::ConstField, 0. );

		fieldc->initField( Pimpact::ConstField, 0. );
		sol->initField( Pimpact::Grad2D_inX, 1. );
		er->initField( Pimpact::ConstField, 0. );
		er->random();

		if( mgSpaces->participating(-3) )
			op1->apply( *fieldff, *fieldf );
		if( mgSpaces->participating(-2) )
			op2->apply( *fieldf, *fieldc );


//		if( mgSpaces->participating(-2) ) fieldf->write( 0 );
//		if( mgSpaces->participating(-1) )fieldc->write( 1 );

		if( mgSpaces->participating(-1) ) {
			er->add( 1., *sol, -1., *fieldc );
			double bla = er->norm(Belos::InfNorm);
			if( 0==space->rankST() )
				std::cout << "error GradX: " << bla << "\n";
//			if( i>0 )
			TEST_EQUALITY( bla<eps, true ); // boundaries?
			if( bla>=eps ) {
				er->write(1);
				fieldc->write(10);
				sol->write(100);
			}
		}

		// the hard test in Y
		fieldff->initField( Pimpact::Grad2D_inY, 1.);
		fieldf->initField( Pimpact::ConstField, 0. );
		fieldc->initField( Pimpact::ConstField, 0. );
		sol->initField( Pimpact::Grad2D_inY, 1. );
		er->initField( Pimpact::ConstField, 0. );
		er->random();

		if( mgSpaces->participating(-3) )
			op1->apply( *fieldff, *fieldf );
		if( mgSpaces->participating(-2) )
			op2->apply( *fieldf, *fieldc );


		if( mgSpaces->participating(-1) ) {
			er->add( 1., *sol, -1., *fieldc );
			double bla = er->norm(Belos::InfNorm);
			if( 0==space->rankST() )
				std::cout << "error GradY: " << er->norm(Belos::InfNorm) << "\n";
//			if( i>0 )
			TEST_EQUALITY( er->norm(Belos::InfNorm)<eps, true ); // boundaries?
			if( bla>=eps ) {
				er->write(2);
				fieldc->write(20);
				sol->write(200);
			}
		}

		// the hard test in Z
		fieldff->initField( Pimpact::Grad2D_inZ, 1.);
		fieldf->initField( Pimpact::ConstField, 0. );
		fieldc->initField( Pimpact::ConstField, 0. );

		sol->initField( Pimpact::Grad2D_inZ, 1. );
		er->initField( Pimpact::ConstField, 0. );
		er->random();

		if( mgSpaces->participating(-3) )
			op1->apply( *fieldff, *fieldf );
		if( mgSpaces->participating(-2) )
			op2->apply( *fieldf, *fieldc );

////		fieldf->write( 6 );
////		fieldc->write( 7 );

		if( mgSpaces->participating(-1) ) {
			er->add( 1., *sol, -1., *fieldc );
			double bla = er->norm(Belos::InfNorm);
			if( 0==space->rankST() )
				std::cout << "error GradZ: " << er->norm(Belos::InfNorm) << "\n";
//			if( i>0 )
			TEST_EQUALITY( er->norm(Belos::InfNorm)<eps, true ); // boundaries?
			if( bla>=eps ) {
				er->write(3);
				fieldc->write(30);
				sol->write(300);
			}
		}
	

*/

} // end of TEUCHOS_UNIT_TEST_TEMPLATE


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Restrictor4D, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Restrictor4D, CS4G )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, Interpolator4D, CS ) {
	
	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// processor grid size
	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

  auto space = Pimpact::createSpace<S,O,4>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, maxGrids );

  	//auto mgTransfers = Pimpact::createMGTransfers<Pimpact::TransferOp,Pimpact::RestrictionOp,Pimpact::InterpolationOp>( mgSpaces );

  auto op = Pimpact::createInterpolationTimeOp<Pimpact::InterpolationOp<CSpace4T> >(mgSpaces->get(-1),mgSpaces->get(-2));

  auto op1 = Pimpact::createInterpolationTimeOp<Pimpact::InterpolationOp<CSpace4T> >(mgSpaces->get(-2),mgSpaces->get(-3));

    //auto fieldf = Pimpact::createScalarField( mgSpaces->get( -2 ), type[i] );
    auto fieldf = Pimpact::createTimeField<Pimpact::ScalarField<CSpace4T>,CSpace4T>( mgSpaces->get( -2 ));
    auto fieldff = Pimpact::createTimeField<Pimpact::ScalarField<CSpace4T>,CSpace4T>( mgSpaces->get( -3 ));
    auto fieldc = Pimpact::createTimeField<Pimpact::ScalarField<CSpace4T>,CSpace4T>( mgSpaces->get( -1 ));
    auto sol = fieldff->clone();
    auto er = fieldff->clone();


	// the zero test
	fieldff->init(1.);
	fieldf->init(1.);
	fieldc->init(0.); //initField( Pimpact::ConstField, 0. );// this will be a loop on the time field
		
                //fieldc->write();
                //fieldf->write(1000);
		//fieldf->print();

		if( mgSpaces->participating(-2) )
				op->apply( *fieldc, *fieldf );
		
		if( mgSpaces->participating(-3) )
                                op1->apply( *fieldf, *fieldff );
		
		if( mgSpaces->participating(-3) )
			TEST_EQUALITY( eps>fieldff->norm(Belos::InfNorm), true );
                        
		if( fieldff->norm(Belos::InfNorm)>=eps )
			std::cout << "error = " << fieldff->norm(Belos::InfNorm) << std::endl;
		if( mgSpaces->participating(-1) )
			TEST_EQUALITY( eps>fieldc->norm(), true );
		
		//fieldc->write();
		//fieldf->write(1000);
		

	// the random test
	fieldc->random();
	fieldf->init(0.);
	fieldff->init(0.);

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );

		if( mgSpaces->participating(-2) )
			op->apply( *fieldc, *fieldf );
		        
	        if( mgSpaces->participating(-3) )
                        op1->apply( *fieldf, *fieldff );

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );


    // the stronger init test
    fieldc->init(1.);
    fieldf->init(0.);
    fieldff->init(0.);
    sol->init(1.);
    er->random();

		if( mgSpaces->participating(-2) )
			op->apply( *fieldc, *fieldf );
	
	       if( mgSpaces->participating(-3) )
                        op1->apply( *fieldf, *fieldff );

		if( mgSpaces->participating(-3) ) {
			er->add( 1., *sol, -1., *fieldff );
			double bla = er->norm(Belos::InfNorm);
			if( 0==space->rankST() )
			std::cout << "error Const: " << bla << "\n";
			TEST_EQUALITY( bla < eps, true  );
			if( bla>=eps )
				er->write(0);
		}

    // gradient tes
    //fieldc->init(1.);
    //fieldf->init(0.);
    //sol->init(1.);
    //er->random();

      //          if( mgSpaces->participating(-2) )
        //                op->apply( *fieldc, *fieldf );

	
	// compounf field test
	//
        //auto fn = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<CSpace4T> >( mgSpaces->get( -2 ) ),
//                                               Pimpact::createTimeField< Pimpact::ScalarField<CSpace4T> >( mgSpaces->get( -2 ) ));
        
	//auto cr = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<CSpace4T> >( mgSpaces->get( -1 ) ),
  //                                             Pimpact::createTimeField< Pimpact::ScalarField<CSpace4T> >( mgSpaces->get( -1 ) ));


        // RHS
        //Pimpact::initVectorTimeField( cr->getVFieldPtr(), Pimpact::ConstVel_inX, 10);
	                
	//if( mgSpaces->participating(-2) )
          //              op->apply( *cr, *fn );
	
	// hardcore test init test in X


        //auto fieldf2 = Pimpact::createTimeField<Pimpact::VectorField<CSpace4T>>( mgSpaces->get( -2 ));
    	//auto fieldc2 = Pimpact::createTimeField<Pimpact::VectorField<CSpace4T>>( mgSpaces->get( -1 ));
	
	//Pimpact::initVectorTimeField( fieldc2, Pimpact::Grad2D_inX );
	//Pimpact::initVectorTimeField( fieldf2, Pimpact::Zero2DFlow );
		
	//fieldc->initField( Pimpact::Grad2D_inX );
	//fieldf->initField( Pimpact::ConstField, 0. );
	//auto sol2 = fieldf2->clone();
        //auto er2 = fieldf2->clone();
	
	//	sol->initField( Pimpact::Grad2D_inX );
	//	er->random();
	//	if( mgSpaces->participating(-1) )
	//		fieldc->write(1001);
//	//	else
//	//		fieldc->print();

	//	if( mgSpaces->participating(-2) )
	//		op->apply( *fieldc, *fieldf );

	// er->add( 1., *sol, -1., *fieldf );

	//	if( mgSpaces->participating(-1) )
	//		fieldc->write(1000);
//	//	else
//	//		fieldc->print();

	//	if( mgSpaces->participating(-2) ){
//	//		if( rankbla==space->rankST() ) {
//	//			er->print();
//	//			std::cout << "rank: " << space->rankST() << " procCord: " << space->	procCoordinate()[0] << ", " << space->	procCoordinate()[1] << ", "<< space->	procCoordinate()[2] << "\n";
//	//		}
//	//		if( rankbla==space->rankST() )
	//		{
	//			auto bla = er->norm( Belos::InfNorm, false );
	//			std::cout << "rank: " << space->rankST() << " procCord: " << space->procCoordinate()[0] << ", " << space->procCoordinate()[1] << ", "<< space->procCoordinate()[2] << ", local error:" << bla << "\n";
//	//			if( bla>1.e-12 ) 
	//			if( rankbla==space->rankST() ) {
	//				er->print();
	//				fieldf->print();
	//				sol->print();
	//				mgSpaces->get(-2)->getIndexSpace()->print();
	//				op->print();
	//			}
	//		}

	//		double bla = er->norm(Belos::InfNorm);
	//		if( 0==space->rankST() )
	//			std::cout << "error GradX: " << bla << "\n";
//	//	if( i>0 ) // corners for scalar are "wrong"
	//			TEST_EQUALITY( bla<eps, true  );
	//		if( bla>=eps ) {
	//			er->write(1);
	//			fieldf->write(10);
	//			sol->write(100);
	//		}
	//	}


	//	// hardcore test init test in Y
	//	fieldc->initField( Pimpact::Grad2D_inY );
	//	fieldf->initField( Pimpact::ConstField, 0. );
	//	sol->initField( Pimpact::Grad2D_inY );
	//	er->random();

	//	if( mgSpaces->participating(-2) )
	//		op->apply( *fieldc, *fieldf );

	//	er->add( 1., *sol, -1., *fieldf );

	//	if( mgSpaces->participating(-2) ) {
	//		double bla = er->norm(Belos::InfNorm);
	//		if( 0==space->rankST() )
	//			std::cout << "error GradY: " <<  bla << "\n";
//	//	if( i>0 ) // corners for scalar are "wrong"
	//		TEST_EQUALITY( bla < eps, true  );
	//		if( bla>=eps ) {
	//			er->write(2);
	//			fieldf->write(20);
	//		}
	//	}

	// // hardcore test init test in Z
	//	fieldc->initField( Pimpact::Grad2D_inZ );
	//	fieldf->initField( Pimpact::ConstField, 0. );
	//	sol->initField( Pimpact::Grad2D_inZ );
	//	er->random();

	//	if( mgSpaces->participating(-2) )
	//		op->apply( *fieldc, *fieldf );

	//	er->add( 1., *sol, -1., *fieldf );


	//	if( mgSpaces->participating(-2) ) {
	//		double bla = er->norm(Belos::InfNorm);
	//		if( 0==space->rankST() )
	//			std::cout << "error GradZ: " << bla << "\n";
//	//	if( i>0 ) // corners for scalar are "wrong"
	//		TEST_EQUALITY( bla < eps, true  );
	//		if( bla>=eps ) {
	//			er->write(3);
	//			fieldf->write(30);
	//		}
	//	}
//	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Interpolator4D, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Interpolator4D, CS4G )


template<class SpaceT> using CVF = Pimpact::CompoundField<Pimpact::TimeField<Pimpact::VectorField<SpaceT> >,
		              			          Pimpact::TimeField<Pimpact::ScalarField<SpaceT> > >;

template<class SpaceT> using INT = Pimpact::IntResCompoundOp<
																				Pimpact::InterpolationTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<SpaceT> > >,
																				Pimpact::InterpolationTimeOp<                           Pimpact::InterpolationOp<SpaceT> > >;

template<class SpaceT> using RES = Pimpact::IntResCompoundOp<
                                        Pimpact::RestrictionTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::RestrictionOp<SpaceT> > >,
                                        Pimpact::RestrictionTimeOp<                           Pimpact::RestrictionOp<SpaceT> > >;

template<class SpaceT1, class SpaceT2> using TCO = Pimpact::TransferCompoundOp<
																				Pimpact::TransferTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<SpaceT1, SpaceT2> > >,
																				Pimpact::TransferTimeOp<                           Pimpact::TransferOp<SpaceT1, SpaceT2> > >; 

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, MG, CS ) {

	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

	auto space = Pimpact::createSpace<S,O,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, 3 );


//	auto pl = Teuchos::parameterList();
	auto mg = Pimpact::createMultiGrid<
									CVF,
									TCO,
									RES,
									INT,
									Pimpact::TimeStokesOp,
									Pimpact::TimeStokesOp,															       Pimpact::TimeStokesBSmoother,
									Pimpact::TimeStokesBSmoother > (mgSpaces);

	mg->print();

auto x = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<FSpace4T> >( space ),
                                       Pimpact::createTimeField< Pimpact::ScalarField<FSpace4T> >( space ));


typedef Pimpact::TimeStokesOp<FSpace4T> OpT;

auto op = Pimpact::create<OpT>( space );

double p = 10;
double alpha = std::sqrt(pl->get<double>("alpha2"));

auto b = x->clone();
auto true_sol = x->clone();
auto err = x->clone();
err->init(1);

Pimpact::initVectorTimeField( true_sol->getVFieldPtr(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha );

x=true_sol->clone();

op->apply(*true_sol,*b);

for( int i=0; i<2; ++i ) {
	mg->apply( *b, *x );
//	x->write(i*100);

	err->add( -1, *x, 1., *true_sol );
	std::cout << "err: " << err->norm() << "\n";
}

TEST_EQUALITY( err->norm()<1.e-5, true );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, MG, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, MG, CS4G )

} // end of namespae
