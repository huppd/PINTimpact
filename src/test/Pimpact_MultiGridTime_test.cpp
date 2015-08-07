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

int domain = 1;
int ftype = 0;

int fs = 0;
int fe = 4;

int npx = 1;
int npy = 1;
int npz = 1;
int npf = 1;

int nx = 17;
int ny = 17;
int nz = 17;
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

  pl->set( "lx", 2. );
  pl->set( "ly", 2. );
  pl->set( "lz", 2. );

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );
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
/*
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
		auto solf = fieldf->clone();
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

    // strong test with cosinus

std::cout << "----------- strng test ------------" << std::endl;

  auto op_v = Pimpact::createRestrictionTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::RestrictionOp<CSpace4T> > >(mgSpaces->get(-5),mgSpaces->get(-4));

  //auto op1_v = Pimpact::createInterpolationTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<CSpace4T> > >(mgSpaces->get(-2),mgSpaces->get(-3));


        auto fieldf_v = Pimpact::createTimeField<Pimpact::VectorField<CSpace4T>,CSpace4T>( mgSpaces->get( -5 ));
        auto fieldc_v = Pimpact::createTimeField<Pimpact::VectorField<CSpace4T>,CSpace4T>( mgSpaces->get( -4));


fieldc_v->init(2);

        Pimpact::initVectorTimeField( fieldf_v, Pimpact::ConstVel_inX, 10);

        if( mgSpaces->participating(-5) )
               op_v->apply( *fieldf_v, *fieldc_v);
      
      fieldc_v->write();
      fieldf_v->write(100); 

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

    // strong test with cosinus

  auto op_v = Pimpact::createInterpolationTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<CSpace4T> > >(mgSpaces->get(-3),mgSpaces->get(-4));

  //auto op1_v = Pimpact::createInterpolationTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<CSpace4T> > >(mgSpaces->get(-2),mgSpaces->get(-3));


	auto fieldf_v = Pimpact::createTimeField<Pimpact::VectorField<CSpace4T>,CSpace4T>( mgSpaces->get( -4 ));
	auto fieldc_v = Pimpact::createTimeField<Pimpact::VectorField<CSpace4T>,CSpace4T>( mgSpaces->get( -3));
        
        Pimpact::initVectorTimeField( fieldc_v, Pimpact::ConstVel_inX, 10);
	                
	if( mgSpaces->participating(-4) )
               op_v->apply( *fieldc_v, *fieldf_v);
//	fieldc_v->write();
//	fieldf_v->write(100);	
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Interpolator4D, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Interpolator4D, CS4G )
*/

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

template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, MG, CS ) {

	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

	auto space = Pimpact::createSpace<S,O,4,4>( pl ); 

	auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, maxGrids );


	auto mgPL = Teuchos::parameterList();
	mgPL->sublist("Smoother").set<std::string>("Solver name", "GMRES" );
	mgPL->sublist("Smoother").set( "Solver",
			*Pimpact::createLinSolverParameter( "GMRES", 1.e-16, -1,
				Teuchos::rcp<std::ostream>( new Teuchos::oblackholestream() ), 4 ) );
	
        mgPL->sublist("Smoother").set<int>( "numIters", 4 );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<int>( "Maximum Iterations", 1000 );
	mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "GMRES" );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver" );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<S>("Convergence Tolerance" , 9.e-3 );

	auto mg = Pimpact::createMultiGrid<
									CVF,
									TCO,
									RES,
									INT,
									Pimpact::TimeStokesOp,
									Pimpact::TimeStokesOp,								
									Pimpact::TimeStokesBSmoother,
						//			Pimpact::TimeStokesBSmoother
									MOP
										> ( mgSpaces, mgPL );

	mg->print();

auto x = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<FSpace4T> >( space ),
                                       Pimpact::createTimeField< Pimpact::ScalarField<FSpace4T> >( space ));


//typedef Pimpact::TimeStokesOp<FSpace4T> OpT;

//auto op = Pimpact::create<OpT>( space );

//double p = 10;
//double alpha = std::sqrt(pl->get<double>("alpha2"));

auto b = x->clone();
auto true_sol = x->clone();
auto err = x->clone();
err->init(1);

//Pimpact::initVectorTimeField( true_sol->getVFieldPtr(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha );

//x=true_sol->clone();

//op->apply(*true_sol,*b);

x->random();
x->scale(10);

err->add( -1, *x, 1., *true_sol );
std::cout << "err: " << err->norm() << "\n";
err->write();

for( int i=0; i<4; ++i ) {
	mg->apply( *b, *x );
	
	x->level();

	err->add( -1, *x, 1., *true_sol );
	err->write((i+1)*100);
	std::cout << "err: " << err->norm() << "\n";
}

TEST_EQUALITY( err->norm()<1.e-5, true );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, MG, CS4L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, MG, CS4G )

} // end of namespae
