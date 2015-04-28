#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_VectorFieldOpWrap.hpp"
#include "Pimpact_MultiGrid.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Operator.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"



namespace {

typedef double S;
typedef int O;

typedef Pimpact::Space<S,O,3,4> FSpace3T;
typedef Pimpact::Space<S,O,4,4> FSpace4T;

typedef Pimpact::Space<S,O,3,2> CSpace3T;
typedef Pimpact::Space<S,O,4,2> CSpace4T;

typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS3L;
typedef Pimpact::CoarsenStrategyGlobal<FSpace3T,CSpace3T,5> CS3G;

template<class ST> using BSF = Pimpact::MultiField< Pimpact::ScalarField<ST> >;
//template<class T> using BVF = Pimpact::MultiField< Pimpact::VectorField<T> >;

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

  pl->set( "dim", 3 );
  pl->set( "domain", 0 );

	pl->set( "nx", 33 );
	pl->set( "ny", 33 );
	pl->set( "nz", 33 );
//	pl->set( "nz", 33 );
//	pl->set( "nx", 65 );
//	pl->set( "ny", 65 );
//	pl->set( "nz", 65 );
//	pl->set( "nx", 49 );
//	pl->set( "ny", 49 );
//	pl->set( "nz", 49 );
//	pl->set( "nx", 129 );
//	pl->set( "ny", 129 );
//	pl->set( "nz", 129 );
//	pl->set( "nx", 1025 );
//	pl->set( "ny", 1025 );
	
// processor grid size
	pl->set("npx", 2 );
	pl->set("npy", 2 );
	pl->set("npz", 2 );

}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSpaces, constructor3D, CS ) {

  auto space = Pimpact::createSpace( pl );

 auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );
 std::cout << "rank: " << space->rankST() << "\tnGridLevels: " << mgSpaces->getNGrids() << "\n";
 if( space->rankST()==0 )
	 mgSpaces->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSpaces, constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSpaces, constructor3D, CS3G )



TEUCHOS_UNIT_TEST( MGSpaces, constructor4D ) {

  typedef Pimpact::CoarsenStrategy<FSpace4T,CSpace4T> CS;

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );

  pl->set( "lx", 1. );
  pl->set( "ly", 1. );
  pl->set( "lz", 1. );

//  pl->set( "dim", 2 );

//  pl->set("nx", 33 );
//  pl->set("ny", 17 );
//  pl->set("nz", 9 );

  pl->set("nf", 65 );

  auto space = Pimpact::createSpace<S,O,4>( pl );

  space->print();

  auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, 10 );
  std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";

 if( 0==space->rankST() )
	mgSpaces->print();

}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGFields, SF_constructor3D, CS ) {

//  typedef Pimpact::CoarsenStrategyGlobal<FSpace3T,CSpace3T> CS;

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );
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

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );
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

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );

  auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradO2Op>( mgSpaces );

  auto mgOps2 = Pimpact::createMGOperators<Pimpact::DivGradOp,Pimpact::DivGradO2Op>( mgSpaces );

  auto op = mgOps->get( -1 );

	if( mgSpaces->participating(-1) )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, SF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, SF_constructor3D, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGOperators, VF_constructor3D, CS ) {

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );

  auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>( mgSpaces );


  auto op = mgOps->get( -1 );

	if( mgSpaces->participating(-1) )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, VF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, VF_constructor3D, CS3G )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSmoothers, SF_constructor3D, CS ) {

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );


  auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradOp,Pimpact::DivGradO2Op>( mgSpaces );

  auto mgSmoother = Pimpact::createMGSmoothers<Pimpact::DivGradO2JSmoother>( mgOps );

  auto op = mgSmoother->get( -1 );

	if( mgSpaces->participating(-1) )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, SF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, SF_constructor3D, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSmoothers, VF_constructor3D, CS ) {

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );


  auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>( mgSpaces );

  auto mgSmoother = Pimpact::createMGSmoothers<ConvDiffSORT>( mgOps );

  auto op = mgSmoother->get( -1 );

	if( mgSpaces->participating(-1) )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, VF_constructor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, VF_constructor3D, CS3G )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, Restrictor3D, CS ) {

	auto space = Pimpact::createSpace<S,O,3>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );

	auto mgTransfers = Pimpact::createMGTransfers<Pimpact::TransferOp,Pimpact::RestrictionOp,Pimpact::InterpolationOp>( mgSpaces );

//	mgTransfers->print();

	Pimpact::EField type[] = {Pimpact::EField::S, Pimpact::EField::U, Pimpact::EField::V, Pimpact::EField::W };

	for( int i=fs; i<fe; ++i ) {
		std::cout << "type: " << i << "\n";

		Teuchos::RCP< Pimpact::ScalarField<CSpace3T> > fieldff;
		Teuchos::RCP< Pimpact::ScalarField<CSpace3T> > fieldf;
		Teuchos::RCP< Pimpact::ScalarField<CSpace3T> > fieldc;
		Teuchos::RCP< Pimpact::ScalarField<CSpace3T> > sol;
		Teuchos::RCP< Pimpact::ScalarField<CSpace3T> > er;

		fieldff = Pimpact::createScalarField( mgSpaces->get( -3 ), type[i] );
		fieldf = Pimpact::createScalarField( mgSpaces->get( -2 ), type[i] );

		fieldc = Pimpact::createScalarField( mgSpaces->get( -1 ), type[i] );
		sol = fieldc->clone();
		er = fieldc->clone();


		Teuchos::RCP<const Pimpact::RestrictionOp<CSpace3T> > op1;
		Teuchos::RCP<const Pimpact::RestrictionOp<CSpace3T> > op2;

		if( mgSpaces->participating(-3) ) {
			std::cout << "\nrank: " << space->rankST() << "\top-3\n";
			op1 = mgTransfers->getRestrictionOp( -3 );
//			op1->print();
		}
		if( mgSpaces->participating(-2) ) {
			std::cout << "\nrank: " << space->rankST() << "\top: 2\n";
			op2 = mgTransfers->getRestrictionOp( -2 );
			if( 1==space->rankST() )
				op2->print();
		}

		// the zero test
		if( mgSpaces->participating(-3) )
			fieldff->init( 0. );
		if( mgSpaces->participating(-2) ) {
			fieldf->initField( Pimpact::ConstField, 0. );
			fieldf->init( 1. );
		}
		if( mgSpaces->participating(-1) )
			fieldc->init( 1. );

		if( mgSpaces->participating(-3) )
			op1->apply( *fieldff, *fieldf );
		if( mgSpaces->participating(-2) )
			op2->apply( *fieldf, *fieldc );

		if( mgSpaces->participating(-2) )
			TEST_FLOATING_EQUALITY( 0., fieldf->norm(), eps );

		if( mgSpaces->participating(-1) ) 
			TEST_FLOATING_EQUALITY( 0., fieldc->norm(), eps );


		// the random test
		fieldff->random();

		if( mgSpaces->participating(-3) )
			op1->apply( *fieldff, *fieldf );
		if( mgSpaces->participating(-2) )
			op2->apply( *fieldf, *fieldc );

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );

		// the strong test
		fieldff->initField( Pimpact::ConstField, 1. );
		fieldf->initField( Pimpact::ConstField, 0. );
		fieldc->initField( Pimpact::ConstField, 0. );
		sol->initField( Pimpact::ConstField, 1. );
		er->initField( Pimpact::ConstField, 0. );
		er->random();

		if( mgSpaces->participating(-3) )
			op1->apply( *fieldff, *fieldf );
		if( mgSpaces->participating(-2) )
			op2->apply( *fieldf, *fieldc );

		if( mgSpaces->participating(-1) ) {

			er->add( 1., *sol, -1., *fieldc );
			std::cout << "const error: " << er->norm() << "\n";
			TEST_EQUALITY( er->norm()<eps, true ); // boundaries?

		}

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
			fieldc->write(1);
			sol->write(2);
			er->add( 1., *sol, -1., *fieldc );
			er->write(4);

			std::cout << "x error: " << er->norm() << "\n";
			TEST_EQUALITY( er->norm()<eps, true ); // boundaries?
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

//		fieldf->write( 4 );
//		fieldc->write( 5 );

		if( mgSpaces->participating(-1) ) {
			er->add( 1., *sol, -1., *fieldc );
			//		er->write(6);

			std::cout << "Y error: " << er->norm() << "\n";
			TEST_EQUALITY( er->norm()<eps, true ); // boundaries?
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

			std::cout << "Z error: " << er->norm() << "\n";
			TEST_EQUALITY( er->norm()<eps, true ); // boundaries?
			er->write(6);
		}

	} // end of for( )

} // end of TEUCHOS_UNIT_TEST_TEMPLATE


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Restrictor3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Restrictor3D, CS3G )




/// \todo remove corners for test(Scalar case)
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, Interpolator3D, CS ) {


  auto space = Pimpact::createSpace<S,O,3>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );

  auto mgTransfers = Pimpact::createMGTransfers<Pimpact::TransferOp,Pimpact::RestrictionOp,Pimpact::InterpolationOp>( mgSpaces );

	auto op = mgTransfers->getInterpolationOp( -2 );
	if( 0==space->rankST() ) op->print();

	Pimpact::EField type[] = { Pimpact::EField::S, Pimpact::EField::U, Pimpact::EField::V, Pimpact::EField::W };

  for( int i=fs; i<fe; ++i ) {
    std::cout << "type: " << i << "\n";

    auto fieldf = Pimpact::createScalarField( mgSpaces->get( -2 ), type[i] );
    auto fieldc = Pimpact::createScalarField( mgSpaces->get( -1 ), type[i] );
    auto sol = fieldf->clone();
    auto er = fieldf->clone();


	 // the zero test
	 fieldf->init( 1. );
	 fieldc->initField( Pimpact::ConstField, 0. );

		if( mgSpaces->participating(-2) )
				op->apply( *fieldc, *fieldf );

		if( mgSpaces->participating(-2) )
			TEST_EQUALITY( eps>fieldf->norm(), true );
		if( mgSpaces->participating(-1) )
			TEST_EQUALITY( eps>fieldc->norm(), true );


	 // the random test
	 fieldc->random();
	 fieldf->init(0.);

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );

		if( mgSpaces->participating(-2) )
			op->apply( *fieldc, *fieldf );

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );


    // the stronger init test
		if( mgSpaces->participating(-1) )
			fieldc->initField( Pimpact::ConstField, 1. );
		else
			fieldc->initField();
    fieldf->initField( Pimpact::ConstField, 0. );
    sol->initField( Pimpact::ConstField, 1. );
		er->random();

		if( mgSpaces->participating(-2) )
			op->apply( *fieldc, *fieldf );

		er->add( 1., *sol, -1., *fieldf );
		if( mgSpaces->participating(-2) ) {
			er->write(0);
			if( 1==space->rankST() ) fieldf->print();
		}

		if( mgSpaces->participating(-2) )
			std::cout << "error Const: " << er->norm()/std::sqrt( (S)er->getLength() ) << "\n";

		if( mgSpaces->participating(-2) )
			TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() ) < eps, true  );


	 // hardcore test init test in X

		if( mgSpaces->participating(-1) )
			fieldc->initField( Pimpact::Grad2D_inX );
		else
			fieldc->initField();
	 fieldf->initField( Pimpact::ConstField, 0. );
	 sol->initField( Pimpact::Grad2D_inX );
		er->random();

		if( mgSpaces->participating(-1) )
			fieldc->write(1001);
//		else
//			fieldc->print();

		if( mgSpaces->participating(-2) )
			op->apply( *fieldc, *fieldf );

	 er->add( 1., *sol, -1., *fieldf );

		if( mgSpaces->participating(-1) )
			fieldc->write(1000);
//		else
//			fieldc->print();

		if( mgSpaces->participating(-2) ){
			er->write(1);
//			if( 1==space->rankST() ) er->print();
			fieldf->write(10);
			sol->write(100);
		}

		if( mgSpaces->participating(-2) )
			std::cout << "error GradX: " << er->norm() << "\n";

		if( i>0 ) // corners for scalar are "wrong"
			if( mgSpaces->participating(-2) )
				TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() )< eps, true  );


	 // hardcore test init test in Y
		if( mgSpaces->participating(-1) )
			fieldc->initField( Pimpact::Grad2D_inY );
		else
			fieldc->initField();
	 fieldf->initField( Pimpact::ConstField, 0. );
	 sol->initField( Pimpact::Grad2D_inY );
		er->random();

		if( mgSpaces->participating(-2) )
			op->apply( *fieldc, *fieldf );

		er->add( 1., *sol, -1., *fieldf );

		if( mgSpaces->participating(-2) ) {
			er->write(2);
			fieldf->write(20);
		}

		if( mgSpaces->participating(-2) )
			std::cout << "error GradY: " << er->norm()/std::sqrt( (S)er->getLength() ) << "\n";

		if( i>0 ) // corners for scalar are "wrong"
			if( mgSpaces->participating(-2) )
				TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() ) < eps, true  );


	 // hardcore test init test in Z
		if( mgSpaces->participating(-1) )
			fieldc->initField( Pimpact::Grad2D_inZ );
		else
			fieldc->initField();
	 fieldf->initField( Pimpact::ConstField, 0. );
	 sol->initField( Pimpact::Grad2D_inZ );
		er->random();

		if( mgSpaces->participating(-2) )
			op->apply( *fieldc, *fieldf );

		er->add( 1., *sol, -1., *fieldf );

		if( mgSpaces->participating(-2) ) er->write(3);
		if( mgSpaces->participating(-2) ) fieldf->write(30);

		if( mgSpaces->participating(-2) )
			std::cout << "error GradZ: " << er->norm()/std::sqrt( (S)er->getLength() ) << "\n";

		if( i>0 ) // corners for scalar are "wrong"
			if( mgSpaces->participating(-2) )
				TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() ) < eps, true  );
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Interpolator3D, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, Interpolator3D, CS3G )



template<class T1,class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, MGTransfersVF, CS ) {


//  pl->set( "Re", 1. );
//  pl->set( "alpha2", 1. );
//
//  pl->set( "dim", 2 );
//
//  pl->set("nf", 0 );
  //  pl->set("nfs", 0 );
  //  pl->set("nfe", 0 );

	auto space = Pimpact::createSpace<S,O,3>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 2 );

	auto mgTransfers = Pimpact::createMGTransfers<
      TransVF,RestrVF,InterVF>( mgSpaces );

	auto fieldf = Pimpact::create<Pimpact::VectorField>( mgSpaces->get( 0 ) );
	auto fieldc = Pimpact::create<Pimpact::VectorField>( mgSpaces->get( 1 ) );

	auto sol = fieldf->clone( Pimpact::ShallowCopy );
	auto er = fieldf->clone( Pimpact::ShallowCopy );

	auto op = mgTransfers->getInterpolationOp( 0 );

	if( space->rankST()==1 ) op->print();


	// the zero test

	fieldf->init( 1. );

	fieldc->initField( Pimpact::ConstFlow, 0., 0., 0. );

	op->apply( *fieldc, *fieldf );


	TEST_EQUALITY( fieldf->norm()<eps, true );
	TEST_EQUALITY( fieldc->norm()<eps, true );


	// the random test
	fieldc->random();
	fieldf->init(0.);

	TEST_INEQUALITY( 0., fieldc->norm() );

	op->apply( *fieldc, *fieldf );

	TEST_INEQUALITY( 0., fieldc->norm() );


	// the stronger init test
	fieldc->initField( Pimpact::ConstFlow, 1., 1., 1. );
	fieldf->init(0.);
	sol->initField( Pimpact::ConstFlow, 1., 1., 1. );

	op->apply( *fieldc, *fieldf );

	er->add( 1., *sol, -1., *fieldf );
	er->write(0);

	std::cout << "error Const: " << er->norm() << "\n";
	TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() )< eps, true  );


	// hardcore test init test in X
	fieldc->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
	fieldc->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX );
	fieldc->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX );
	fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );

	sol->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
	sol->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX );
	sol->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX );

	op->apply( *fieldc, *fieldf );

	er->add( 1., *sol, -1., *fieldf );
	er->write(1);

	std::cout << "error GradX: " << er->norm() << "\n";

	TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() )< eps, true  );

	// hardcore test init test in Y
	fieldc->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY );
	fieldc->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
	fieldc->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inY );
	fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );

	sol->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY );
	sol->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
	sol->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inY );

	op->apply( *fieldc, *fieldf );

	er->add( 1., *sol, -1., *fieldf );
	er->write(2);

	std::cout << "error GradY: " << er->norm()/std::sqrt( (S)er->getLength() ) << "\n";

	TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() ) < eps, true  );

	// hardcore test init test in Z
	fieldc->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inZ );
	fieldc->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inZ );
	fieldc->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );
	fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );

	sol->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inZ );
	sol->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inZ );
	sol->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );

	op->apply( *fieldc, *fieldf );

	er->add( 1., *sol, -1., *fieldf );
	er->write(2);

	std::cout << "error GradZ: " << er->norm()/std::sqrt( (S)er->getLength() ) << "\n";

	TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() ) < eps, true  );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, MGTransfersVF, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, MGTransfersVF, CS3G )


template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, DivGradOp, CS ) {

	auto space = Pimpact::createSpace( pl );

	auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );
	if( space->rankST()==0 ) mgSpaces->print();


	auto mg =
		Pimpact::createMultiGrid<
		Pimpact::ScalarField,
		Pimpact::TransferOp,
		Pimpact::RestrictionOp,
		Pimpact::InterpolationOp,
		Pimpact::DivGradOp,
		Pimpact::DivGradO2Op,
		Pimpact::DivGradO2JSmoother,
		MOP>( mgSpaces );

	auto x = Pimpact::create<Pimpact::ScalarField>( space );
	auto b = Pimpact::create<Pimpact::ScalarField>( space );
	auto op = Pimpact::create<Pimpact::DivGradOp>( space );

	std::ofstream ofs;
	if( space()->rankST()==0 )
		ofs.open("MG.txt", std::ofstream::out);

	//	Grad in x
	x->initField( Pimpact::Grad2D_inZ );
	x->write(0);

	op->apply(*x,*b);
	b->write(1);

	x->init( 0. );
	auto e = x->clone();

	e->init( 1. );
	auto sol = x->clone();
	sol->initField( Pimpact::Grad2D_inZ );

	S bla = b->dot(*e)/x->getLength();
	std::cout<< " rhs nullspace: " << bla << "\n";

	for( int i=0; i<10; ++i ) {
		mg->apply( *b, *x );
		x->level();
//		x->write(i+10);

		op->apply(*x,*sol);
		sol->add( -1, *b, 1., *sol );
		S res = sol->norm();

		if( space()->rankST()==0 ) {
			std::cout << "res: " << res << "\n";
			ofs << res << "\n";
		}
	}
	TEST_EQUALITY( sol->norm()<1.e-3, true );

//	x->write(2);

	if( space()->rankST()==0 )
		ofs.close();

	auto xm = Pimpact::createMultiField( x );
	auto bm = Pimpact::createMultiField( b );

	xm->init(0.);
	bm->level();


	auto prec = Pimpact::createMultiOperatorBase( mg );

	auto solvName = "Block GMRES";
//	auto solvName = "GMRES";

	auto param = Pimpact::createLinSolverParameter(solvName,1.e-8);
	param->set( "Output Frequency", 1);
	param->set( "Maximum Iterations", 50 );
	param->set( "Flexible Gmres", true );

	auto bop = Pimpact::createMultiOperatorBase( op );

	auto linprob = Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::ScalarField<FSpace3T> > >( bop, xm, bm, param, solvName );
	linprob->setRightPrec(prec);

	linprob->solve(xm,bm);
	xm->write();
//   auto bop = Pimpact::createMultiOperatorBase( mg );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, DivGradOp, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, DivGradOp, CS3G )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, ConvDiffOp, CS ) {

	pl->set<S>( "Re", 1000 );
	//pl->set<S>( "Re", 100. );
	//pl->set<S>( "Re", 10. );
	//pl->set<S>( "Re", 1. );
	//pl->set<S>( "Re", 0.01 );

  //pl->set( "lx", 2. );
  //pl->set( "ly", 2. );
//	pl->set( "nx", 65 );
//	pl->set( "ny", 65 );


  //auto space = Pimpact::createSpace<S,O,3,2>( pl );
  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );

	auto mgPL = Teuchos::parameterList();
	mgPL->sublist("Smoother").set( "omega", 1. );
	mgPL->sublist("Smoother").set( "numIters", 4 );
	mgPL->sublist("Smoother").set( "Ordering", 1 );
	mgPL->sublist("Smoother").set<short int>( "dir X", -1 );
	mgPL->sublist("Smoother").set<short int>( "dir Y", -1 );


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

	std::cout << "blabla\n";
  auto op = Pimpact::create< ConvDiffOpT >( space );

	std::cout << "blupblup\n";
  auto x = Pimpact::create<Pimpact::VectorField>( space );
  auto b = Pimpact::create<Pimpact::VectorField>( space );
  auto temp = x->clone();

  {
    auto wind = x->clone();
		//wind->initField( Pimpact::ConstFlow, 0., 0., 0. );
		wind->initField( Pimpact::ConstFlow, 1., 1., 1. );
    op->assignField( *wind );
    mg->assignField( *wind );
  }

  std::ofstream ofs;
  if( space()->rankST()==0 )
    ofs.open("MG2.txt", std::ofstream::out);

  // 
	x->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inY );
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

   for( int i=0; i<20; ++i ) {
     mg->apply( *b, *x );
     //x->write(i+10);

     temp->add( -1, *x, 1., *sol );
     S res = temp->norm();

     if( space()->rankST()==0 ) {
			 std::cout << "res: " << res << "\n";
       ofs << res << "\n";
		 }
   }

   TEST_EQUALITY( temp->norm()<0.5, true );

   x->write(2);

   if( space()->rankST()==0 )
     ofs.close();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffOp, CS3L )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffOp, CS3G )




} // end of namespace
