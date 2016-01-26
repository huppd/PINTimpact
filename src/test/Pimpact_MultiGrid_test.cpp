#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_Operator.hpp"




namespace {



using S = double;
using O = int;

const int dimension = 4;

using FSpace3T = Pimpact::Space<S,O,dimension,4>;
using CSpace3T = Pimpact::Space<S,O,dimension,2>;

using CS3L = Pimpact::CoarsenStrategy<FSpace3T,CSpace3T>;
using CS3G = Pimpact::CoarsenStrategyGlobal<FSpace3T,CSpace3T>;

template<class ST> using BSF = Pimpact::MultiField< Pimpact::ScalarField<ST> >;

template<class ST> using BOPF = Pimpact::MultiOpWrap< Pimpact::DivGradOp<ST> >;
template<class ST> using BOPC = Pimpact::MultiOpWrap< Pimpact::DivGradO2Op<ST> >;
template<class ST> using BSM = Pimpact::MultiOpWrap< Pimpact::DivGradO2JSmoother<ST> >;

template<class T> using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;

template<class T> using ConvDiffSORT = Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionSORSmoother >;

template<class T> using ConvDiffJT = Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionJSmoother >;

template<class T1,class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionHWOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

bool testMpi = true;
double eps = 1e-6;

int dim = 3;
int domain = 0;

int ftype = 0;

int fs = 0;
int fe = 4;

int level = -1;

S lx = 30.;
S ly = 20.;
S lz = 20.;

int npx = 1;
int npy = 1;
int npz = 1;
int npf = 1;

O nx = 65;
O ny = 49;
O nz = 97;
O nf = 1;

int rankbla = -1;

int maxGrids = 10;

auto pl = Teuchos::parameterList();



TEUCHOS_STATIC_SETUP() {

	Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
	clp.addOutputSetupOptions(true);
	clp.setOption( "test-mpi", "test-serial", &testMpi,
			"Test MPI (if available) or force test of serial.  In a serial build,"
			" this option is ignored and a serial comm is always used." );
	clp.setOption( "error-tol-slack", &eps,
			"Slack off of machine epsilon used to check test results" );

	clp.setOption( "domain", &domain, "domain" );
	clp.setOption( "dim", &dim, "dim" );

	clp.setOption( "ftype", &ftype,
			"Slack off of machine epsilon used to check test results" );
	clp.setOption( "fs", &fs,
			"Slack off of machine epsilon used to check test results" );
	clp.setOption( "fe", &fe,
			"Slack off of machine epsilon used to check test results" );
	clp.setOption( "level", &level, "" );

	clp.setOption( "lx", &lx, "" );
	clp.setOption( "ly", &ly, "" );
	clp.setOption( "lz", &lz, "" );

	clp.setOption( "npx", &npx, "" );
	clp.setOption( "npy", &npy, "" );
	clp.setOption( "npz", &npz, "" );
	clp.setOption( "npf", &npf, "" );

	clp.setOption( "nx", &nx, "" );
	clp.setOption( "ny", &ny, "" );
	clp.setOption( "nz", &nz, "" );
	clp.setOption( "nf", &nf, "" );
	clp.setOption( "rank", &rankbla, "" );
	clp.setOption( "maxGrids", &maxGrids, "" );

	////pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "para" );
	//pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "cos" );
	//pl->sublist("Stretching in X").set<S>( "N metr L", nx );
	//pl->sublist("Stretching in X").set<S>( "N metr U", nx  );
	//pl->sublist("Stretching in X").set<S>( "x0 L", 0.05 );
	////pl->sublist("Stretching in X").set<S>( "x0 L", 0. );
	//pl->sublist("Stretching in X").set<S>( "x0 U", 0. );

	////pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "cos" );
	////pl->sublist("Stretching in Y").set<S>( "N metr L", 1 );
	////pl->sublist("Stretching in Y").set<S>( "N metr U", ny  );
	////pl->sublist("Stretching in Y").set<S>( "x0 L", 0. );
	////pl->sublist("Stretching in Y").set<S>( "x0 U", 0. );

	////pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "cos" );
	////pl->sublist("Stretching in Z").set<S>( "N metr L", 1 );
	////pl->sublist("Stretching in Z").set<S>( "N metr U", nz  );
	////pl->sublist("Stretching in Z").set<S>( "x0 L", 0. );
	////pl->sublist("Stretching in Z").set<S>( "x0 U", 0. );

	pl->set<S>( "Re", 400 );

}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSpaces, constructor3D, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	int rank = space->rankST();
	space->print();

	{
		Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "coord_rank_"+std::to_string(rank)+".txt", 0 );

		space->getCoordinatesGlobal()->print( *fstream );
	}

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids);

	std::cout << "rank: " << space->rankST() << "\tnGridLevels: " << mgSpaces->getNGrids() << "\n";

	if( space->rankST()==0 )
		mgSpaces->print();

	for( int level=0; level<mgSpaces->getNGrids(); ++level ) {
		Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "coord_l"+std::to_string(level)+"rank_"+std::to_string(rank)+".txt",0 );
		mgSpaces->get(level)->getCoordinatesGlobal()->print( *fstream );
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSpaces, constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSpaces, constructor3D, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGFields, SF_constructor3D, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set( "nx", nx );
	pl->set( "ny", ny );
	pl->set( "nz", nz );
	pl->set( "nf", nf );

	// processor grid size
	pl->set( "npx", npx );
	pl->set( "npy", npy );
	pl->set( "npz", npz );
	pl->set( "npf", npf );

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );
	std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
	if( space->rankST()==0 )
		mgSpaces->print();

	auto mgFields = Pimpact::createMGFields<Pimpact::ScalarField>( mgSpaces );

	auto field = mgFields->get( -1 );
	if( mgSpaces->participating(-1) ){

		field->init(5.);
		TEST_FLOATING_EQUALITY( std::sqrt( std::pow(5.,2)*field->getLength() ), field->norm(Belos::TwoNorm), eps );

		field->write(0);
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, SF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, SF_constructor3D, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGFields, VF_constructor3D, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );
	std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
	if( space->rankST()==0 )
		mgSpaces->print();

	auto mgFields = Pimpact::createMGFields<Pimpact::VectorField>( mgSpaces );

	auto field = mgFields->get( -1 );
	if( mgSpaces->participating(-1) ){

		field->init(5.);
		TEST_FLOATING_EQUALITY( std::sqrt( std::pow(5.,2)*field->getLength() ), field->norm(Belos::TwoNorm), eps );

		field->write(0);
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, VF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, VF_constructor3D, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGOperators, SF_constructor3D, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );

	auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradO2Op>( mgSpaces );

	auto mgOps2 = Pimpact::createMGOperators<Pimpact::DivGradOp,Pimpact::DivGradO2Op>( mgSpaces );

	for( int i=0; i<3; ++i ) {
		auto op = mgOps->get( i );
		if( 0==space->rankST() )
			std::cout << "\n\n --- level: " << i << " ---\n\n";

		if( mgSpaces->participating(i) )
			op->print();
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, SF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, SF_constructor3D, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGOperators, VF_constructor3D, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );

	auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>( mgSpaces );


	auto op = mgOps->get( -1 );

	if( mgSpaces->participating(-1) )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, VF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, VF_constructor3D, CS3G )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSmoothers, SF_constructor3D, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );


	auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradOp,Pimpact::DivGradO2Op>( mgSpaces );

	auto mgSmoother = Pimpact::createMGSmoothers<Pimpact::DivGradO2JSmoother>( mgOps );

	auto op = mgSmoother->get( -1 );

	if( mgSpaces->participating(-1) )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, SF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, SF_constructor3D, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSmoothers, VF_constructor3D, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );

	auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>( mgSpaces );

	auto mgSmoother = Pimpact::createMGSmoothers<ConvDiffSORT>( mgOps );

	auto op = mgSmoother->get( -1 );

	if( mgSpaces->participating(-1) )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, VF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, VF_constructor3D, CS3G )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGTransfers, Restrictor, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );

	for( int level=1; level<mgSpaces->getNGrids(); ++level ) {
		if( 0==space->rankST() ) {
			std::cout << "\n\n\t--- level: " << level-1 << "---\n";
		}
		Teuchos::RCP<const Pimpact::RestrictionHWOp<CSpace3T> > op = 
			Teuchos::rcp( new Pimpact::RestrictionHWOp<CSpace3T>(
						mgSpaces->get(level-1),
						mgSpaces->get(level),
						mgSpaces->get()->getProcGrid()->getNP() ) );

		//if( mgSpaces->participating(level-1) ) op->print();

		Teuchos::Tuple<Pimpact::EField,4> type =
			Teuchos::tuple(
					Pimpact::EField::S,
					Pimpact::EField::U,
					Pimpact::EField::V,
					Pimpact::EField::W );

		for( int i=fs; i<fe; ++i ) {
			if( 2==dim && i==Pimpact::EField::W ) break;
			if( 0==space->rankST() )
				std::cout << " --- ftype: " << i << " ---\n";

			Teuchos::RCP< Pimpact::ScalarField<CSpace3T> > fieldf;
			Teuchos::RCP< Pimpact::ScalarField<CSpace3T> > fieldc;
			Teuchos::RCP< Pimpact::ScalarField<CSpace3T> > sol;
			Teuchos::RCP< Pimpact::ScalarField<CSpace3T> > er;

			fieldf = Pimpact::createScalarField( mgSpaces->get( level-1 ), type[i] );

			fieldc = Pimpact::createScalarField( mgSpaces->get( level ), type[i] );
			sol = fieldc->clone();
			er = fieldc->clone();

			// the zero test
			fieldf->initField( Pimpact::ConstField, 0. );
			fieldc->init( 1. );

			if( mgSpaces->participating(level-1) ) op->apply( *fieldf, *fieldc );

			if( mgSpaces->participating(level-1) )
				TEST_FLOATING_EQUALITY( 0., fieldf->norm(), eps );

			if( mgSpaces->participating(level) ) 
				TEST_FLOATING_EQUALITY( 0., fieldc->norm(), eps );

			// the random test
			fieldf->random();

			if( mgSpaces->participating(level-1) ) op->apply( *fieldf, *fieldc );

			if( mgSpaces->participating(level) )
				TEST_INEQUALITY( 0., fieldc->norm() );

			// the const test
			fieldf->initField( Pimpact::ConstField, 1. );
			fieldc->initField( Pimpact::ConstField, 0. );
			sol->initField( Pimpact::ConstField, 1. );

			if( mgSpaces->participating(level-1) ) op->apply( *fieldf, *fieldc );

			if( mgSpaces->participating(level) ) {

				er->add( 1., *sol, -1., *fieldc );
				S errInf = er->norm(Belos::InfNorm);
				if( 0==space->rankST() )
					std::cout << "error Const: " << errInf << " ("<< op->getDD() << ")\n";
				//if( i>0 )
				TEST_EQUALITY( errInf<eps, true ); // boundaries?
				if( errInf>=eps )
					er->write(0);
			}

			// the hard test
			for( int dir=3; dir<6; ++ dir ) {
				fieldf->initField( static_cast<Pimpact::EScalarField>(dir) );
				fieldc->initField( Pimpact::ConstField, 0. );
				sol->initField( static_cast<Pimpact::EScalarField>(dir) );

				if( mgSpaces->participating(level-1) ) op->apply( *fieldf, *fieldc );

				if( mgSpaces->participating(level) ) {
					er->add( 1., *sol, -1., *fieldc );
					double errInf = er->norm(Belos::InfNorm);
					if( 0==space->rankST() )
						std::cout << "error Grad in "<< dir-3<< "-dir: " << errInf << " ("<< op->getDD() << ")\n";
					TEST_EQUALITY( errInf<eps, true );
					if( errInf>=eps ) {
						er->print();
						er->write(1*(dir-2));
						fieldc->write(10*(dir-2));
						sol->write(100*(dir-2));
					}
				}
			}
		}
	} // end of for( )

} // end of TEUCHOS_UNIT_TEST_TEMPLATE


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, Restrictor, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, Restrictor, CS3G )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGTransfers, Interpolator, CS ) {

	pl->set( "domain", domain );
	pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	Teuchos::RCP<const Pimpact::Space<S,O,dimension,4> > space =
		Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );

	std::cout << "\nrank: " << space->rankST() << "\tnGridLevels: " << mgSpaces->getNGrids() << "\n";

	auto mgTransfers = Pimpact::createMGTransfers<Pimpact::TransferOp,Pimpact::RestrictionHWOp,Pimpact::InterpolationOp>( mgSpaces );

	for( int level=0; level<mgSpaces->getNGrids()-1; ++level ) {

		Teuchos::RCP< const Pimpact::InterpolationOp<CSpace3T> > op = 
			Teuchos::rcp( new Pimpact::InterpolationOp<CSpace3T>(
						mgSpaces->get(level+1),
						mgSpaces->get(level),
						mgSpaces->get()->getProcGrid()->getNP() ) );
		if( 0==space->rankST() ) {
			std::cout << "\n\n\t--- level: " << level << "---\n";
		}
		if( space->rankST()==rankbla ) {
			op->print();
		}

		Teuchos::Tuple<Pimpact::EField,4> type =
			Teuchos::tuple(
					Pimpact::EField::S,
					Pimpact::EField::U,
					Pimpact::EField::V,
					Pimpact::EField::W );

		for( int i=fs; i<fe; ++i ) {
			if( 0==space->rankST() )
				std::cout << "field type: " << i << "\n";

			auto fieldf = Pimpact::createScalarField( mgSpaces->get( level   ), type[i] );
			auto fieldc = Pimpact::createScalarField( mgSpaces->get( level+1 ), type[i] );
			auto sol = fieldf->clone();
			auto er = fieldf->clone();


			// the zero test
			fieldf->init( 1. );
			fieldc->initField( Pimpact::ConstField, 0. );

			if( mgSpaces->participating(level) )
				op->apply( *fieldc, *fieldf );

			if( mgSpaces->participating(level) )
				TEST_EQUALITY( eps>fieldf->norm(Belos::InfNorm), true );
			if( mgSpaces->participating(level+1) )
				TEST_EQUALITY( eps>fieldc->norm(Belos::InfNorm), true );

			// the random test
			fieldc->random();
			fieldf->init(0.);

			if( mgSpaces->participating(level+1) )
				TEST_INEQUALITY( 0., fieldc->norm() );

			if( mgSpaces->participating(level) )
				op->apply( *fieldc, *fieldf );

			if( mgSpaces->participating(level+1) )
				TEST_INEQUALITY( 0., fieldc->norm() );


			// the stronger init test
			fieldc->initField( Pimpact::ConstField, 1. );
			fieldf->initField( Pimpact::ConstField, 0. );
			sol->initField( Pimpact::ConstField, 1. );
			er->random();

			if( mgSpaces->participating(level) )
				op->apply( *fieldc, *fieldf );


			er->add( 1., *sol, -1., *fieldf );
			if( mgSpaces->participating(level) ) {
				S errInf = er->norm( Belos::InfNorm );
				if( 0==space->rankST() )
					std::cout << "error Const: " << errInf << "\n";
				TEST_EQUALITY( errInf < eps, true  );
				if( errInf>=eps )
					er->write(0);
			}
			if( er->norm(Belos::InfNorm, false )>=eps && space->rankST()==rankbla ){
				//			std::cout << "rank: " << space->rankST() << "\n";
				er->print();
			}


			// --- hardcore test ---
			for( int dir=3; dir<6; ++dir ) {
				fieldc->initField( static_cast<Pimpact::EScalarField>(dir) );
				fieldf->initField( Pimpact::ConstField, 0. );
				sol->initField( static_cast<Pimpact::EScalarField>(dir) );
				er->random();

				if( mgSpaces->participating(level) ){
					op->apply( *fieldc, *fieldf );

					er->add( 1., *sol, -1., *fieldf );

					S errInf = er->norm( Belos::InfNorm );
					if( 0==space->rankST() )
						std::cout << "error Grad in " << dir-3 << "-dir: " << errInf << "\n";
					TEST_EQUALITY( errInf<eps, true  );
					if( errInf>=eps ) {
						int i = dir;
						er->write( i );
						fieldf->write( 10*i );
						sol->write( 100*i );
						er->print();
					}
				}
			}
		}
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, Interpolator, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, Interpolator, CS3G )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGTransfers, MGTransfersSF, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );

	auto mgTransfers = Pimpact::createMGTransfers<
		Pimpact::TransferOp,Pimpact::RestrictionHWOp,Pimpact::InterpolationOp>(
				mgSpaces );

	auto x = Pimpact::createMGFields<Pimpact::ScalarField>( mgSpaces );

	auto fieldf = x->get();
	auto fieldc = x->get(-1);

	S errInf = 0.;

	// interpolation 
	{
		auto sol = fieldf->clone( Pimpact::ShallowCopy );
		auto er = fieldf->clone( Pimpact::ShallowCopy );

		// the zero test
		fieldf->init( 1. );
		fieldc->initField();

		mgTransfers->interpolation( x );

		if( mgSpaces->participating(-1) ) {
			errInf = fieldc->norm( Belos::InfNorm );
			TEST_EQUALITY( errInf<eps, true );
		}
		if( mgSpaces->participating(0) ) {
			errInf = fieldf->norm( Belos::InfNorm );
			TEST_EQUALITY( errInf<eps, true );
			if( 0==space->rankST() )
				std::cout << "\ninterpolation error zero: " << errInf << "\n";
		}

		// the random test
		fieldc->random();
		fieldf->init(0.);

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm( Belos::InfNorm ) );

		mgTransfers->interpolation( x );

		if( mgSpaces->participating(0) )
			TEST_INEQUALITY( 0., fieldf->norm( Belos::InfNorm ) );

		// the Const test
		fieldc->initField( Pimpact::ConstField, 1. );
		fieldf->initField();
		sol->initField( Pimpact::ConstField, 1. );

		mgTransfers->interpolation( x );

		er->add( 1., *sol, -1., *fieldf );

		if( mgSpaces->participating(0) ) {
			errInf = er->norm( Belos::InfNorm );
			if( 0==space->rankST() )
					std::cout << "interpolation error Const: " << errInf << "\n";
			TEST_EQUALITY( errInf<eps, true  );
			if( errInf>=eps ) {
				int i = 0;
				er->write( i );
				std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
				er->print( *Pimpact::createOstream( "int_error_c_r"+r+".txt" ) );
				//er->print(  );
			}
		}


		// hardcore test
		for( int dir=3; dir<6; ++dir ) {
			fieldc->initField( static_cast<Pimpact::EScalarField>(dir) );
			fieldf->initField();
			sol->initField( static_cast<Pimpact::EScalarField>(dir) );

			mgTransfers->interpolation( x );

			er->add( 1., *sol, -1., *fieldf );

			errInf = er->norm( Belos::InfNorm );
			if( 0==space->rankST() )
				std::cout << "interpolation error Grad in " << dir-3 << ": " << errInf << "\n";
			TEST_EQUALITY( errInf<eps, true  );
			if( errInf>=eps ) {
				int i = dir-2;
				er->write( i );
				std::string d = std::to_string( static_cast<long long>( dir-3 ) ); // long long needed on brutus(intel)
				std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
				er->print( *Pimpact::createOstream( "int_error_g"+d+"_r"+r+".txt" ) );
				//er->print(  );
			}
		}

	}

	// restriction 
	{
		auto sol = fieldc->clone( Pimpact::ShallowCopy );
		auto er = fieldc->clone( Pimpact::ShallowCopy );


		// the zero test
		fieldf->initField();
		fieldc->init( 1. );

		mgTransfers->restriction( x );

		if( mgSpaces->participating(0) )
			TEST_EQUALITY( fieldf->norm()<eps, true );
		if( mgSpaces->participating(-1) )
			TEST_EQUALITY( fieldc->norm()<eps, true );

		// the random test
		fieldf->random();
		fieldc->init(0.);

		if( mgSpaces->participating(0) )
			TEST_INEQUALITY( 0., fieldf->norm() );

		mgTransfers->restriction( x );

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );

		// the const test
		fieldf->initField( Pimpact::ConstField, 1. );
		fieldc->initField();
		sol->initField( Pimpact::ConstField, 1. );

		mgTransfers->restriction( x );

		er->add( 1., *sol, -1., *fieldc );

		if( mgSpaces->participating(-1) ) {
			errInf = er->norm( Belos::InfNorm );
			TEST_EQUALITY( errInf<eps, true  );
		}

		if( 0==space->rankST() )
			std::cout << "\nrestriction error Const: " << errInf << "\n";


		// hardcore Grad test
		for( int dir=3; dir<6; ++ dir ) {
			fieldf->initField( static_cast<Pimpact::EScalarField>(dir) );
			fieldc->initField();

			sol->initField( static_cast<Pimpact::EScalarField>(dir) );

			mgTransfers->restriction( x );

			er->add( 1., *sol, -1., *fieldc );
			if( mgSpaces->participating(-1) ) {

				errInf = er->norm( Belos::InfNorm );
				TEST_EQUALITY( errInf < eps, true  );
				if( errInf>=eps )
					er->write( 1*(dir-2) );
			}
			if( 0==space->rankST() )
				std::cout << "restriction error Grad in " << dir-3 << ": " << errInf << "\n";
		}

	}

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, MGTransfersSF, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, MGTransfersSF, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGTransfers, MGTransfersVF, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

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

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );

	auto mgTransfers = Pimpact::createMGTransfers<
		TransVF,RestrVF,InterVF>( mgSpaces );

	auto x = Pimpact::createMGFields<Pimpact::VectorField>( mgSpaces );

	auto fieldf = x->get();
	auto fieldc = x->get(-1);

	// interpolation 
	{
		auto sol = fieldf->clone( Pimpact::ShallowCopy );
		auto er = fieldf->clone( Pimpact::ShallowCopy );


		// the zero test
		fieldf->init( 1. );
		fieldc->initField( Pimpact::ConstFlow, 0., 0., 0. );

		mgTransfers->interpolation( x );

		if( mgSpaces->participating(0) )
			TEST_EQUALITY( fieldf->norm()<eps, true );
		if( mgSpaces->participating(-1) )
			TEST_EQUALITY( fieldc->norm()<eps, true );


		// the random test
		fieldc->random();
		fieldf->init(0.);

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );

		mgTransfers->interpolation( x );

		if( mgSpaces->participating(0) )
			TEST_INEQUALITY( 0., fieldf->norm() );

		// the stronger init test
		fieldc->initField( Pimpact::ConstFlow, 1., 1., 1. );
		fieldf->init(0.);
		sol->initField( Pimpact::ConstFlow, 1., 1., 1. );

		mgTransfers->interpolation( x );

		er->add( 1., *sol, -1., *fieldf );

		if( mgSpaces->participating(0) ) {
			S rel_error = er->norm( Belos::InfNorm );
			if( 0==space->rankST() )
				std::cout << "\nint. error Const: " << rel_error << "\n";
			TEST_EQUALITY( rel_error < eps, true  );
			if( rel_error>=eps || isnan(rel_error) ) {
				er->write(0);
			}
		}


		//// hardcore test init test in X
		//fieldc->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
		//fieldc->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX );
		//fieldc->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX );
		//fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );

		//sol->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
		//sol->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX );
		//sol->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX );

		//mgTransfers->interpolation( x );

		//er->add( 1., *sol, -1., *fieldf );

		//if( mgSpaces->participating(0) ) {
			//S rel_error = er->norm( Belos::InfNorm );
			//if( 0==space->rankST() )
				//std::cout << "int. error GradX: " << rel_error << "\n";
			//TEST_EQUALITY( rel_error<eps, true  );
			//if( rel_error>=eps || isnan(rel_error) ) {
				//er->write(1);
				//fieldf->write(10);
				//sol->write(100);
			//}
		//}

		//// hardcore test init test in Y
		//fieldc->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY );
		//fieldc->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
		//fieldc->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inY );
		//fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );

		//sol->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY );
		//sol->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
		//sol->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inY );

		//mgTransfers->interpolation( x );

		//er->add( 1., *sol, -1., *fieldf );
		//if( mgSpaces->participating(0) ) {
			//S rel_error = er->norm( Belos::InfNorm );
			//if( 0==space->rankST() )
				//std::cout << "int. error GradY: " << rel_error << "\n";
			//TEST_EQUALITY( rel_error<eps, true  );
			//if( rel_error>=eps || isnan(rel_error) ) {
				//er->write(2);
				//fieldf->write(20);
				//sol->write(200);
			//}
		//}

		// hardcore Grad test
		for( int dir=3; dir<6; ++ dir ) {
			fieldc->getFieldPtr( Pimpact::U )->initField(static_cast<Pimpact::EScalarField>(dir) );
			fieldc->getFieldPtr( Pimpact::V )->initField(static_cast<Pimpact::EScalarField>(dir) );
			fieldc->getFieldPtr( Pimpact::W )->initField(static_cast<Pimpact::EScalarField>(dir) );
			fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );

			sol->getFieldPtr( Pimpact::U )->initField( static_cast<Pimpact::EScalarField>(dir) );
			sol->getFieldPtr( Pimpact::V )->initField( static_cast<Pimpact::EScalarField>(dir) );
			sol->getFieldPtr( Pimpact::W )->initField( static_cast<Pimpact::EScalarField>(dir) );

			mgTransfers->interpolation( x );

			er->add( 1., *sol, -1., *fieldf );

			if( mgSpaces->participating(0) ) {
				S rel_error = er->norm( Belos::InfNorm );
				if( 0==space->rankST() )
					std::cout << "interpolatio. error in "<<dir-3<<"-dir: " << rel_error << "\n";
				TEST_EQUALITY( rel_error < eps, true  );
				if( rel_error>=eps || isnan(rel_error) ) {
					er->write(1*(dir-2));
					fieldf->write(10*(dir-2));
					sol->write(100*(dir-2));
				}
			}
		}

	}
	// restriction 
	{
		auto sol = fieldc->clone( Pimpact::ShallowCopy );
		auto er = fieldc->clone( Pimpact::ShallowCopy );

		// the zero test
		fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );
		fieldc->init( 1. );

		mgTransfers->restriction( x );

		if( mgSpaces->participating(0) )
			TEST_EQUALITY( fieldf->norm()<eps, true );
		if( mgSpaces->participating(-1) )
			TEST_EQUALITY( fieldc->norm()<eps, true );


		// the random test
		fieldf->random();
		fieldc->init(0.);

		if( mgSpaces->participating(0) )
			TEST_INEQUALITY( 0., fieldf->norm() );

		mgTransfers->restriction( x );

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );


		// the stronger init test
		fieldf->initField( Pimpact::ConstFlow, 1., 1., 1. );
		fieldc->init(0.);
		sol->initField( Pimpact::ConstFlow, 1., 1., 1. );

		mgTransfers->restriction( x );

		er->add( 1., *sol, -1., *fieldc );
		if( mgSpaces->participating(-1) ) {
			S rel_error = er->norm()/std::sqrt( (S)er->getLength() );
			std::cout << "res. error Const: " << rel_error << "\n";
			TEST_EQUALITY( rel_error<eps, true  );
			if( rel_error>eps ) {
				er->write(0);
			}
		}

		// hardcore grad test
		for( int dir=3; dir<6; ++ dir ) {
			fieldf->getFieldPtr( Pimpact::U )->initField( static_cast<Pimpact::EScalarField>(dir) );
			fieldf->getFieldPtr( Pimpact::V )->initField( static_cast<Pimpact::EScalarField>(dir) );
			fieldf->getFieldPtr( Pimpact::W )->initField( static_cast<Pimpact::EScalarField>(dir) );
			fieldc->initField( Pimpact::ConstFlow, 0., 0., 0. );

			sol->getFieldPtr( Pimpact::U )->initField( static_cast<Pimpact::EScalarField>(dir) );
			sol->getFieldPtr( Pimpact::V )->initField( static_cast<Pimpact::EScalarField>(dir) );
			sol->getFieldPtr( Pimpact::W )->initField( static_cast<Pimpact::EScalarField>(dir) );

			mgTransfers->restriction( x );

			er->add( 1., *sol, -1., *fieldc );
			if( mgSpaces->participating(-1) ) {
				S rel_error = er->norm()/std::sqrt( (S)er->getLength() );
				std::cout << "restriction grad  error in " << dir-3 << "-dir: " << rel_error << "\n";
				TEST_EQUALITY( rel_error < eps, true  );
				if( rel_error>eps ) {
					er->write(dir-2);
					fieldf->write(10*(dir-2));
					sol->write(100*(dir-2));
				}
			}
		}

	}

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, MGTransfersVF, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, MGTransfersVF, CS3G )


template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;
template<class T> using POP = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;
template<class T> using POP2 = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2SORSmoother >;




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, DivGradOp, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set( "nx", nx );
	pl->set( "ny", ny );
	pl->set( "nz", nz );
	pl->set( "nf", nf );

	// processor grid size
	pl->set( "npx", npx );
	pl->set( "npy", npy );
	pl->set( "npz", npz );
	pl->set( "npf", npf );

	Teuchos::RCP<const Pimpact::Space<S,O,dimension,4> > space = Pimpact::createSpace<S,O,dimension,4>( pl );

	int nMax = 10;

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );
	if( space->rankST()==0 ) mgSpaces->print();

	auto mgPL = Teuchos::parameterList();
	mgPL->set<int>( "numCycles", 1 );
	mgPL->set<bool>( "defect correction", false );
	//mgPL->set<bool>( "defect correction", true );
	mgPL->set<bool>( "init zero", false );

	//mgPL->sublist("Smoother").set( "omega", 1. );
	mgPL->sublist("Smoother").set<int>( "numIters", 5 );
	//mgPL->sublist("Smoother").set<int>( "RBGS mode", 0 );
	//mgPL->sublist("Smoother").set<bool>( "level", false );

	//mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "Fixed Point" );
	//mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "CG" );
	//mgPL->sublist("Coarse Grid Solver").set<std::string>( "Solver name", "Flexible GMRES" );
	mgPL->sublist("Coarse Grid Solver").set<std::string>( "Solver name", "GMRES" );
	//mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "TFQMR" );
	//mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver" );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<S>( "Convergence Tolerance" , 2.e-1 );

	//mgPL->sublist("Coarse Grid Solver").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
			//"Output Stream", Teuchos::rcp( &std::cout, false ) );
	//mgPL->sublist("Coarse Grid Solver").sublist("Solver").set("Verbosity",
			//Belos::Errors +
			//Belos::Warnings +
			////Belos::IterationDetails +
			////Belos::OrthoDetails +
			//Belos::FinalSummary +
			////Belos::TimingDetails +
			////Belos::StatusTestDetails +
			//Belos::Debug
			//);

	//mgPL->sublist("Coarse Grid Solver").set<S>( "omega", 1. );
	//mgPL->sublist("Coarse Grid Solver").set<int>( "RBGS mode", 2 );
	//mgPL->sublist("Coarse Grid Solver").set<int>( "numIters", 12 );
	//mgPL->sublist("Coarse Grid Solver").set<bool>( "level", false );
	
	mgPL->sublist("Coarse Grid Solver").sublist("Preconditioner").set<S>( "omega", 1. );
	mgPL->sublist("Coarse Grid Solver").sublist("Preconditioner").set<int>( "RBGS mode", 0 );
	mgPL->sublist("Coarse Grid Solver").sublist("Preconditioner").set<int>( "numIters", 1 );
	mgPL->sublist("Coarse Grid Solver").sublist("Preconditioner").set<bool>( "level", false );
	//mgPL->sublist("Coarse Grid Solver").sublist("Preconditioner").set<bool>( "level", true );

	auto mg =
		Pimpact::createMultiGrid<
			Pimpact::ScalarField,
			Pimpact::TransferOp,
			Pimpact::RestrictionHWOp,
			Pimpact::InterpolationOp,
			Pimpact::DivGradOp,
			Pimpact::DivGradO2Op,
			Pimpact::DivGradO2JSmoother,
			//Pimpact::DivGradO2SORSmoother,
			POP2 >( mgSpaces, mgPL );
			//POP >( mgSpaces, mgPL );
			//MOP >( mgSpaces, mgPL );
			//Pimpact::DivGradO2SORSmoother >( mgSpaces, mgPL );

	auto x = Pimpact::create<Pimpact::ScalarField>( space );
	auto b = Pimpact::create<Pimpact::ScalarField>( space );
	auto res = x->clone( Pimpact::ShallowCopy );
	auto sol = x->clone( Pimpact::ShallowCopy );

	auto xm = Pimpact::createMultiField( x );
	auto bm = Pimpact::createMultiField( b );

	auto op   = Pimpact::create<Pimpact::DivGradOp>( space );
	auto opO2 = Pimpact::create<Pimpact::DivGradO2Op>( space );

	auto prec = Pimpact::createMultiOperatorBase( mg );

	//auto solvName = "Block GMRES";
	 //auto solvName = "TFQMR";
	auto solvName = "GMRES";
	//auto solvName = "Fixed Point";

	auto param = Pimpact::createLinSolverParameter( solvName, 1.e-6 );
	param->set( "Output Frequency", 10 );
	param->set( "Maximum Iterations", 10 );
	//param->set( "Flexible Gmres", false );

	auto bop = Pimpact::createMultiOperatorBase( op );

	auto linprob = Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::ScalarField<FSpace3T> > >( bop, xm, bm, param, solvName );
	linprob->setRightPrec(prec);
	//linprob->setLeftPrec(prec);
	
	// --- zero rhs test ---
	b->init(0.);
	x->random();
	x->setCornersZero();
	
	S error0 = x->norm();
	S errorp = error0;

	opO2->apply( *x, *res );
	S res0 = res->norm();
	S resP = res0;

	if( space()->rankST()==0 ) {
		std::cout << "\n\n\t\t\t--- zero rhs test ---\n";
		std::cout << "\tresidual:\trate:\t\t\terror:\t\trate: \n";
		std::cout <<  std::scientific;
		std::cout << "\t" << 1. << "\t\t\t\t"  << 1.  << "\n";
	}

	for( int i=0; i<nMax; ++i ) {
		mg->apply( *b, *x );
		x->level();

		S error = x->norm();

		opO2->apply( *x, *res );
		S residual = res->norm();

		if( space()->rankST()==0 )
			std::cout << "\t" << residual/res0 << "\t" << residual/resP << "\t\t" << error/error0 << "\t" <<  error/errorp << "\n";

		//if( error>= errorp )
			//break;
		//else
		errorp = error;
		resP = residual;
	}
	
	//x->print();
	x->write();
	TEST_EQUALITY( x->norm()/std::sqrt( static_cast<S>(x->getLength()) )<1.e-3, true );
	
	//bm->init( 0. );
	//xm->random();
	//linprob->solve( xm, bm );

	// --- grad test ---
	auto e = x->clone();
	for( int dir=3; dir<6; ++dir ) {
		x->initField( static_cast<Pimpact::EScalarField>(dir) );
		sol->initField( static_cast<Pimpact::EScalarField>(dir) );
		sol->level();

		opO2->apply( *x, *b );
		//b->write(1);

		//x->random();
		x->init( 0 );

		// residual
		opO2->apply( *x, *res );
		res->add( -1, *b, 1., *res );
		S res0 = res->norm();
		S resP = res0;
		//error
		res->add( 1., *sol, -1., *x );
		S err0 = res->norm();
		S errP = err0;
		//S err0 = 1.;
		//S errP =1.;

		if( space()->rankST()==0 ) {
			std::cout << "\n\n\t\t\t--- " << dir << " test ---\n";
			std::cout << "\tresidual:\trate:\t\t\terror:\t\trate:\n";
			std::cout << "\t"  << 1.<< "\t\t\t\t" << 1.  << "\n";
		}
		for( int i=0; i<nMax; ++i ) {
			mg->apply( *b, *x );
			x->level();

			opO2->apply( *x, *res );
			res->add( -1, *b, 1., *res );
			S residual = res->norm();
			res->add( 1., *sol, -1., *x );
			S err = res->norm();

			if( space()->rankST()==0 )
				std::cout << "\t" << residual/res0 << "\t" <<  residual/resP << "\t\t" << err/err0 << "\t" <<  err/errP  << "\n";
			resP = residual;
			errP = err;
		}
		TEST_EQUALITY( res->norm()/std::sqrt( static_cast<S>( res->getLength() ) )<1.e-3, true );

		////bm->get( 0 )-init( 0. );
		//xm->init( 0. );
		//linprob->solve( xm, bm );
	}


}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, DivGradOp, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, DivGradOp, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, ConvDiffSOR, CS ) {

	pl->set( "domain", domain );
	pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set( "nx", nx );
	pl->set( "ny", ny );
	pl->set( "nz", nz );
	pl->set( "nf", nf );

	// processor grid size
	pl->set( "npx", npx );
	pl->set( "npy", npy );
	pl->set( "npz", npz );
	pl->set( "npf", npf );

	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );

	auto mgPL = Teuchos::parameterList();
	mgPL->set<int>("numCycles", 1);
	mgPL->sublist("Smoother").set( "omega", 1. );
	mgPL->sublist("Smoother").set( "numIters", 2 );
	mgPL->sublist("Smoother").set( "Ordering", 1 );
	mgPL->sublist("Smoother").set<short int>( "dir X", -1 );
	mgPL->sublist("Smoother").set<short int>( "dir Y", -1 );
	mgPL->sublist("Smoother").set<short int>( "dir Z", -1 );

	mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "GMRES" );

//	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<Teuchos::RCP<std::ostream> >( "Output Stream", Teuchos::rcp( &std::cout, false ) );
//	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set("Verbosity",
//			Belos::Errors + Belos::Warnings + Belos::IterationDetails +
//			Belos::OrthoDetails + Belos::FinalSummary + Belos::TimingDetails +
//			Belos::StatusTestDetails + Belos::Debug );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver" );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<S>("Convergence Tolerance"
			, 0.1 );

	auto mg =
		Pimpact::createMultiGrid<
		Pimpact::VectorField,
		TransVF,
		RestrVF,
		InterVF,
		ConvDiffOpT,
		ConvDiffOpT,
		ConvDiffSORT,
		MOP
			>( mgSpaces, mgPL );

	auto op = Pimpact::create< ConvDiffOpT >( space );

	auto x = Pimpact::create<Pimpact::VectorField>( space );
	auto b = Pimpact::create<Pimpact::VectorField>( space );
	auto temp = x->clone();

	{
		auto wind = x->clone();
		wind->initField( Pimpact::ConstFlow, 0., 0., 0. );
//		wind->initField( Pimpact::ConstFlow, 1., 1., 1. );
		op->assignField( *wind );
		mg->assignField( *wind );
	}

	std::ofstream ofs;
	if( space()->rankST()==0 )
		ofs.open("MG2.txt", std::ofstream::out);

	// 
	x->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inX );
	x->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inY );
	x->getFieldPtr(Pimpact::W)->initField( Pimpact::Grad2D_inZ );
	auto sol = x->clone( Pimpact::DeepCopy );

	op->apply(*x,*b);
	{
		x->init(0);
		auto bc = x->clone( Pimpact::ShallowCopy );
		op->apply( *x, *bc );
		b->add( 1., *b, -1., *bc );
	}
	b->write(1);

	x->initField( Pimpact::ConstFlow, 0., 0., 0. );
	//x->random();

	temp->add( -1, *x, 1., *sol );
	S res = temp->norm();
	S res_0 = res;
	S res_p = res;

	if( space()->rankST()==0 ) {
		std::cout << "\n\n\tresidual:\trate: \n";
		std::cout << "\t"  << 1.  << "\n";
		ofs << res << "\n";
	}

	for( int i=0; i<20; ++i ) {
		mg->apply( *b, *x );
		//x->write(i+10);

		temp->add( -1, *x, 1., *sol );
		S res = temp->norm();

		if( space()->rankST()==0 ) {
			std::cout << "\t" << res/res_0 << "\t" <<  res/res_p << "\n";
			ofs << res << "\n";
		}
		res_p = res;
	}

	if( space()->rankST()==0 ) {
		std::cout << "\n";
	}

	TEST_EQUALITY( temp->norm()<0.5, true );

	x->write(2);

	if( space()->rankST()==0 )
		ofs.close();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffSOR, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffSOR, CS3G )



//TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiGrid, ConvDiffJ, CS, SmootherT ) {
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, ConvDiffJ, CS ) {

  pl->set( "domain", domain );
  pl->set( "dim", dim );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	pl->set<S>("Re", 100. );

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

	//auto space = Pimpact::createSpace<S,O,3,2>( pl );
	auto space = Pimpact::createSpace<S,O,dimension,4>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, maxGrids );

	auto mgPL = Teuchos::parameterList();
	mgPL->set<int>("numCycles", 1 );
//	mgPL->sublist("Smoother").set( "omega", 0.6 );
//	mgPL->sublist("Smoother").set( "numIters", 10 );

	mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "GMRES" );

//	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<Teuchos::RCP<std::ostream> >( "Output Stream", Teuchos::rcp( &std::cout, false ) );
//	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set("Verbosity", Belos::Errors + Belos::Warnings +
//			Belos::IterationDetails + Belos::OrthoDetails +
//			Belos::FinalSummary + Belos::TimingDetails +
//			Belos::StatusTestDetails + Belos::Debug );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver" );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<S>("Convergence Tolerance" , 1.e-12 );

	auto mg =
		Pimpact::createMultiGrid<
		Pimpact::VectorField,
		TransVF,
		RestrVF,
		InterVF,
		ConvDiffOpT,
		ConvDiffOpT,
		ConvDiffJT,
//		ConvDiffSORT,
		MOP
			>( mgSpaces, mgPL );

	mg->print();

	auto op = Pimpact::create< ConvDiffOpT >( space );

	auto x = Pimpact::create<Pimpact::VectorField>( space );
	auto b = Pimpact::create<Pimpact::VectorField>( space );
	auto temp = x->clone();

	{
		auto wind = x->clone();
		wind->initField( Pimpact::ConstFlow, 0., 0., 0. );
//		wind->initField( Pimpact::ConstFlow, 1., 1., 1. );
//		wind->random();
		wind->getFieldPtr(Pimpact::U)->initField( Pimpact::Poiseuille2D_inX );
		wind->getFieldPtr(Pimpact::V)->random();
		wind->getFieldPtr(Pimpact::V)->scale(0.1);
		wind->getFieldPtr(Pimpact::W)->random();
		wind->getFieldPtr(Pimpact::W)->scale(0.1);
		op->assignField( *wind );
		mg->assignField( *wind );
	}

	std::ofstream ofs;
	if( space()->rankST()==0 )
		ofs.open("MG2.txt", std::ofstream::out);

	// 
//	x->getFieldPtr(Pimpact::V)->initField( Pimpact::Poiseuille2D_inY );
//	x->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inX );
//	x->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inY );
	x->getFieldPtr(Pimpact::W)->initField( Pimpact::Grad2D_inY );
	auto sol = x->clone( Pimpact::DeepCopy );

//	op->apply(*x,*b);
//	{
//		x->init(0);
//		auto bc = x->clone( Pimpact::ShallowCopy );
//		op->apply( *x, *bc );
//		b->add( 1., *b, -1., *bc );
//	}
//	b->write(1);

	x->initField( Pimpact::ConstFlow, 0., 0., 0. );
	sol->initField( Pimpact::ConstFlow, 0., 0., 0. );
	x->random();

	temp->add( -1, *x, 1., *sol );
	S res = temp->norm();
	S res_0 = res;
	S res_p = res;

	if( space()->rankST()==0 ) {
		std::cout << "\n\n\tresidual:\trate: \n";
		std::cout << "\t"  << 1.  << "\n";
		ofs << res << "\n";
	}

	for( int i=0; i<20; ++i ) {
		mg->apply( *b, *x );
		//x->write(i+10);

		temp->add( -1, *x, 1., *sol );
		S res = temp->norm();

		if( space()->rankST()==0 ) {
			std::cout << "\t" << res/res_0 << "\t" <<  res/res_p << "\n";
			ofs << res << "\n";
		}
		res_p = res;
	}

	if( space()->rankST()==0 ) {
		std::cout << "\n";
	}

	TEST_EQUALITY( temp->norm()<0.5, true );

	x->write(2);

	if( space()->rankST()==0 )
		ofs.close();

}

//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, ConvDiffJ, CS3L, ConvDiffJT )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffJ, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffJ, CS3G )
//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, ConvDiffJ, CS3G, ConvDiffJT )
//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, ConvDiffJ, CS3L, ConvDiffSORT )
//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, ConvDiffJ, CS3G, ConvDiffSORT )



} // end of namespace
