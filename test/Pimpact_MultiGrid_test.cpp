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

int domain = 1;
int ftype = 0;

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

  pl->set( "domain", 1 );
	pl->set( "nx", 65 );
	pl->set( "ny", 65 );
	//pl->set( "nx", 1025 );
	//pl->set( "ny", 1025 );
	
	// processor grid size
  pl->set("npx", 1 );
  pl->set("npy", 1 );

}



TEUCHOS_UNIT_TEST( MGSpaces, constructor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );
  std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
  if( space->rankST()==0 )
    mgSpaces->print();

}



TEUCHOS_UNIT_TEST( MGSpaces, constructor4D ) {

  typedef Pimpact::CoarsenStrategy<FSpace4T,CSpace4T> CS;

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );

  pl->set( "lx", 1. );
  pl->set( "ly", 1. );
  pl->set( "lz", 1. );

  pl->set( "dim", 2 );

  pl->set("nx", 129 );
  pl->set("ny", 65 );
  pl->set("nz", 2 );

  pl->set("nf", 256 );
  //  pl->set("nfs", 0 );
  //  pl->set("nfe", 0 );

  auto space = Pimpact::createSpace<S,O,4>( pl );

  space->print();

  auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, 10 );
  std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
  if( 0==space->rankST() )
    mgSpaces->print();

}



TEUCHOS_UNIT_TEST( MGFields, SF_constructor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );
  std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
  if( space->rankST()==0 )
    mgSpaces->print();

  auto mgFields = Pimpact::createMGFields<Pimpact::ScalarField>( mgSpaces );

  auto field = mgFields->get( 2 );
  field->random();
  field->write(0);

}



TEUCHOS_UNIT_TEST( MGFields, VF_constructor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );
  std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
  if( space->rankST()==0 )
    mgSpaces->print();

  auto mgFields = Pimpact::createMGFields<Pimpact::VectorField>( mgSpaces );

  auto field = mgFields->get( -1 );
  field->random();
  field->write(0);

}



TEUCHOS_UNIT_TEST( MGOperators, SF_constructor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );

  auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradO2Op>( mgSpaces );

  auto mgOps2 = Pimpact::createMGOperators<Pimpact::DivGradOp,Pimpact::DivGradO2Op>( mgSpaces );

  auto op = mgOps->get( 2 );

  op->print();
  //  field->random();
  //  field->write(0);

}




TEUCHOS_UNIT_TEST( MGOperators, VF_constructor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );

  auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>( mgSpaces );


  auto op = mgOps->get( 2 );

  op->print();
  //  field->random();
  //  field->write(0);

}



TEUCHOS_UNIT_TEST( MGSmoothers, SF_constructor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;


  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );


  auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradOp,Pimpact::DivGradO2Op>( mgSpaces );

  auto mgSmoother = Pimpact::createMGSmoothers<Pimpact::DivGradO2JSmoother>( mgOps );

  auto op = mgSmoother->get( 2 );
  op->print();

}



TEUCHOS_UNIT_TEST( MGSmoothers, VF_constructor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;


  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );


  auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>( mgSpaces );

  auto mgSmoother = Pimpact::createMGSmoothers<ConvDiffSORT>( mgOps );

  auto op = mgSmoother->get( 2 );
  op->print();

}



TEUCHOS_UNIT_TEST( MultiGrid, Restrictor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  auto space = Pimpact::createSpace<S,O,3>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 2 );

  auto mgTransfers = Pimpact::createMGTransfers<Pimpact::TransferOp,Pimpact::RestrictionOp,Pimpact::InterpolationOp>( mgSpaces );

  Pimpact::EField type[] = {Pimpact::EField::S, Pimpact::EField::U, Pimpact::EField::V };

  for( int i=0; i<3; ++i ) {
    std::cout << "type: " << i << "\n";

    auto fieldf = Pimpact::createScalarField( mgSpaces->get( 0 ), type[i] );
    auto fieldc = Pimpact::createScalarField( mgSpaces->get( 1 ), type[i] );

    auto op = mgTransfers->getRestrictionOp( 0 );

	op->print();

    // the zero test
    fieldf->init( 0. );
    fieldc->init( 1. );

    op->apply( *fieldf, *fieldc );

    TEST_FLOATING_EQUALITY( 0., fieldf->norm(), eps );
    TEST_FLOATING_EQUALITY( 0., fieldc->norm(), eps );

    // the random test
    fieldf->random();

    TEST_INEQUALITY( 0., fieldf->norm() );

    op->apply( *fieldf, *fieldc );

    TEST_INEQUALITY( 0., fieldc->norm() );

    // the strong test
    fieldf->initField( Pimpact::ConstField, 1. );

    TEST_FLOATING_EQUALITY( 1., fieldf->norm(Belos::InfNorm), eps );

    TEST_FLOATING_EQUALITY( (S)fieldf->getLength(), fieldf->norm(Belos::OneNorm), eps  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldf->getLength() ), fieldf->norm(Belos::TwoNorm), eps  );

    fieldc->init( 0. );

    op->apply( *fieldf, *fieldc );

    fieldf->write( 0 );
    fieldc->write( 1 );

    TEST_FLOATING_EQUALITY( 1., fieldc->norm(Belos::InfNorm), eps );

    TEST_FLOATING_EQUALITY( (S)fieldc->getLength(), fieldc->norm(Belos::OneNorm), eps  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldc->getLength() ), fieldc->norm(Belos::TwoNorm), eps  );

    // the hard test in X
    fieldf->initField( Pimpact::Grad2D_inX, 1. );

    fieldc->initField( Pimpact::ConstField, 0. );
    auto sol = fieldc->clone();
		auto er = fieldc->clone();

    sol->initField( Pimpact::Grad2D_inX, 1. );
    er->initField( Pimpact::ConstField, 0. );

    op->apply( *fieldf, *fieldc );

    fieldf->write( 0 );
    fieldc->write( 1 );

    er->add( 1., *sol, -1., *fieldc );
    er->write(3);

		TEST_EQUALITY( er->norm()<eps, true ); // boundaries?

    // the hard test in Y
    fieldf->initField( Pimpact::Grad2D_inY, 1. );

    fieldc->initField( Pimpact::ConstField, 0. );
    sol = fieldc->clone();
	er = fieldc->clone();

    sol->initField( Pimpact::Grad2D_inY, 1. );
    er->initField( Pimpact::ConstField, 0. );

    op->apply( *fieldf, *fieldc );

    fieldf->write( 4 );
    fieldc->write( 5 );

    er->add( 1., *sol, -1., *fieldc );
    er->write(6);

	TEST_EQUALITY( er->norm()<eps, true ); // boundaries?
		

  }

}



/// \todo remove corners for test(Scalar case)
TEUCHOS_UNIT_TEST( MultiGrid, Interpolator3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );
  pl->set( "domain", domain );

  pl->set( "dim", 2 );

  pl->set("nf", 0 );
  //  pl->set("nfs", 0 );
  //  pl->set("nfe", 0 );

  auto space = Pimpact::createSpace<S,O,3>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 2 );

  auto mgTransfers = Pimpact::createMGTransfers<Pimpact::TransferOp,Pimpact::RestrictionOp,Pimpact::InterpolationOp>( mgSpaces );

  Pimpact::EField type[] = { Pimpact::EField::S, Pimpact::EField::U, Pimpact::EField::V };

  for( int i=1; i<3; ++i ) {
    std::cout << "type: " << i << "\n";

    auto fieldf = Pimpact::createScalarField( mgSpaces->get( 0 ), type[i] );
    auto fieldc = Pimpact::createScalarField( mgSpaces->get( 1 ), type[i] );

    auto op = mgTransfers->getInterpolationOp( 0 );

    if( space->rankST()==0 )
      op->print();

    // the zero test

    fieldf->init( 1. );

    fieldc->initField( Pimpact::ConstField, 0. );

    op->apply( *fieldc, *fieldf );

    //    fieldf->write(0);
    //    fieldc->write(1);

    TEST_EQUALITY( eps>fieldf->norm(), true );
    TEST_EQUALITY( eps>fieldc->norm(), true );

    //   fieldc->print();
    //   fieldf->print();

    // the random test
    fieldc->random();
    fieldf->init(0.);

    TEST_INEQUALITY( 0., fieldc->norm() );

    op->apply( *fieldc, *fieldf );

    fieldf->write(0);
    fieldc->write(1);

    TEST_INEQUALITY( 0., fieldc->norm() );


    // the stronger init test
    fieldc->initField( Pimpact::ConstField, 1. );
    fieldf->init(0.);

    op->apply( *fieldc, *fieldf );

    fieldf->write(2);
    fieldc->write(3);

    TEST_FLOATING_EQUALITY( 1., fieldc->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( (S)fieldc->getLength(), fieldc->norm(Belos::OneNorm), eps  );
    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldc->getLength() ), fieldc->norm(Belos::TwoNorm), eps  );

    TEST_FLOATING_EQUALITY( 1., fieldf->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( (S)fieldf->getLength(), fieldf->norm(Belos::OneNorm), eps  );
    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldf->getLength() ), fieldf->norm(Belos::TwoNorm), eps  );

    // hardcore test init test in X
    fieldc->initField( Pimpact::Grad2D_inX );
    fieldf->initField( Pimpact::ConstField, 0. );
    auto sol = fieldf->clone();
    auto er = fieldf->clone();
    sol->initField( Pimpact::Grad2D_inX );

    op->apply( *fieldc, *fieldf );
    fieldf->write(4);
    fieldc->write(5);

    er->add( 1., *sol, -1., *fieldf );
    er->write(90);

    std::cout << "error GradX: " << er->norm() << "\n";

    TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() )< eps, true  );



    // hardcore test init test in Y
    fieldc->initField( Pimpact::Grad2D_inY );
    fieldf->initField( Pimpact::ConstField, 0. );
    sol->initField( Pimpact::Grad2D_inY );
    sol->write(81);

    op->apply( *fieldc, *fieldf );

    fieldf->write(6);
    fieldc->write(7);

    er->add( 1., *sol, -1., *fieldf );
    er->write(91);

    std::cout << "error GradY: " << er->norm()/std::sqrt( (S)er->getLength() ) << "\n";

    TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() ) < eps, true  );

  }

}




template<class T1,class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

TEUCHOS_UNIT_TEST( MultiGrid, MGTransfersVF ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  pl->set( "Re", 1. );
  pl->set( "alpha2", 1. );
  pl->set( "domain", domain );

  pl->set( "dim", 2 );

  pl->set("nf", 0 );
  //  pl->set("nfs", 0 );
  //  pl->set("nfe", 0 );

  auto space = Pimpact::createSpace<S,O,3>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 2 );

  auto mgTransfers = Pimpact::createMGTransfers<
      TransVF,RestrVF,InterVF>( mgSpaces );

    auto fieldf = Pimpact::create<Pimpact::VectorField>( mgSpaces->get( 0 ) );
    auto fieldc = Pimpact::create<Pimpact::VectorField>( mgSpaces->get( 1 ) );

    auto op = mgTransfers->getInterpolationOp( 0 );

    if( space->rankST()==0 )
      op->print();

    // the zero test

    fieldf->init( 1. );

    fieldc->initField( Pimpact::ConstFlow, 0., 0., 0. );

    op->apply( *fieldc, *fieldf );

    //    fieldf->write(0);
    //    fieldc->write(1);

    TEST_EQUALITY( fieldf->norm()<eps, true );
    TEST_EQUALITY( fieldc->norm()<eps, true );

    //   fieldc->print();
    //   fieldf->print();

    // the random test
    fieldc->random();
    fieldf->init(0.);

    TEST_INEQUALITY( 0., fieldc->norm() );

    op->apply( *fieldc, *fieldf );

//    fieldf->write(0);
//    fieldc->write(1);

    TEST_INEQUALITY( 0., fieldc->norm() );


    // the stronger init test
    fieldc->initField( Pimpact::ConstFlow, 1., 1., 1. );
    fieldf->init(0.);

    op->apply( *fieldc, *fieldf );

    //    fieldf->write(2);
    //    fieldc->write(3);

    TEST_FLOATING_EQUALITY( 1., fieldc->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( (S)fieldc->getLength(), fieldc->norm(Belos::OneNorm), eps  );
    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldc->getLength() ), fieldc->norm(Belos::TwoNorm), eps  );

    TEST_FLOATING_EQUALITY( 1., fieldf->norm(Belos::InfNorm), eps );
    TEST_FLOATING_EQUALITY( (S)fieldf->getLength(), fieldf->norm(Belos::OneNorm), eps  );
    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldf->getLength() ), fieldf->norm(Belos::TwoNorm), eps  );

    // hardcore test init test in X
    fieldc->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
    fieldc->getFieldPtr(1)->initField( Pimpact::Grad2D_inX );
    //fieldc->getFieldPtr(1)->initField( Pimpact::ConstField, 0. );
    fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );

    auto sol = fieldf->clone( Pimpact::ShallowCopy );
    auto er = fieldf->clone( Pimpact::ShallowCopy );

    sol->getFieldPtr(0)->initField( Pimpact::Grad2D_inX );
    sol->getFieldPtr(1)->initField( Pimpact::Grad2D_inX );

    op->apply( *fieldc, *fieldf );
    //    fieldf->write(4);
    //    fieldc->write(5);

    er->add( 1., *sol, -1., *fieldf );
    er->write(90);

    std::cout << "error GradX: " << er->norm() << "\n";

    TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() )< eps, true  );

    // hardcore test init test in Y
    fieldc->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
    fieldc->getFieldPtr(1)->initField( Pimpact::Grad2D_inY );
    //fieldc->getFieldPtr(0)->initField( Pimpact::ConstField, 0. );

    sol->getFieldPtr(0)->initField( Pimpact::Grad2D_inY );
    sol->getFieldPtr(1)->initField( Pimpact::Grad2D_inY );

    fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );


    //    sol->write(81);

    op->apply( *fieldc, *fieldf );

    //    fieldf->write(6);
    //    fieldc->write(7);

    er->add( 1., *sol, -1., *fieldf );
    er->write(91);

    std::cout << "error GradY: " << er->norm()/std::sqrt( (S)er->getLength() ) << "\n";

    TEST_EQUALITY( er->norm()/std::sqrt( (S)er->getLength() ) < eps, true  );

}



template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;


TEUCHOS_UNIT_TEST( MultiGrid, DivGradOp ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  pl->set("nx", 33 );
  pl->set("ny", 65 );
  pl->set("nz", 2 );

  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 5 );


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

  // Grad in x
   x->initField( Pimpact::Grad2D_inX );
   x->write(0);

   op->apply(*x,*b);
   b->write(1);

   x->init( 0. );
//   x->random();
   auto e = x->clone();

   e->init( 1. );
   auto sol = x->clone();
//   auto  = x->clone();
   sol->initField( Pimpact::Grad2D_inX );

   S bla = b->dot(*e)/x->getLength();
   std::cout<< " rhs nullspace: " << bla << "\n";
   for( int i=0; i<40; ++i ) {
     mg->apply( *b, *x );
     x->level();
     x->write(i+10);

     op->apply(*x,*sol);
     sol->add( -1, *b, 1., *sol );
     double res = sol->norm();
     std::cout << "res: " << res << "\n";

     if( space()->rankST()==0 )
       ofs << res << "\n";
   }
   TEST_EQUALITY( sol->norm()<0.5, true );

   x->write(2);

   if( space()->rankST()==0 )
     ofs.close();

   auto xm = Pimpact::createMultiField( x );
   auto bm = Pimpact::createMultiField( b );

   xm->init(0.);


   auto prec = Pimpact::createMultiOperatorBase( mg );

//   auto solvName = "Block GMRES";
   auto solvName = "GMRES";

   auto param = Pimpact::createLinSolverParameter(solvName,1.e-6);
   param->set( "Output Frequency", 1);
//   param->set( "Flexible Gmres", true );

   auto bop = Pimpact::createMultiOperatorBase( op );

   auto linprob = Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::ScalarField<FSpace3T> > >( bop, xm, bm, param, solvName );
   linprob->setRightPrec(prec);

   linprob->solve(xm,bm);
   xm->write();
//   auto bop = Pimpact::createMultiOperatorBase( mg );

}




TEUCHOS_UNIT_TEST( MultiGrid, ConvDiffOp ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

	pl->set<S>( "Re", 1000 );
	//pl->set<S>( "Re", 100. );
	//pl->set<S>( "Re", 10. );
	//pl->set<S>( "Re", 1. );
	//pl->set<S>( "Re", 0.01 );

	//pl->set( "nx", 256 );
	//pl->set( "ny", 256 );
	pl->set( "nx", 129 );
	pl->set( "ny", 129 );
  //pl->set( "nx", 65 );
  //pl->set( "ny", 65 );
  //pl->set( "nx", 33 );
  //pl->set( "ny", 33 );

  //pl->set( "lx", 2. );
  //pl->set( "ly", 2. );

	pl->set("npx", 1 );
	pl->set("npy", 1 );
	pl->set("npz", 1 );

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
				//ConvDiffJT,
				//ConvDiffSORT
				MOP
		> ( mgSpaces, mgPL );

  auto op = Pimpact::create< ConvDiffOpT >( space );

  auto x = Pimpact::create<Pimpact::VectorField>( space );
  auto b = Pimpact::create<Pimpact::VectorField>( space );
  auto temp = x->clone();

  {
    auto wind = x->clone();
		//wind->initField( Pimpact::ConstFlow, 0., 0., 0. );
		wind->initField( Pimpact::ConstFlow, 1., 1., 0. );
    op->assignField( *wind );
    mg->assignField( *wind );
  }

  std::ofstream ofs;
  if( space()->rankST()==0 )
    ofs.open("MG2.txt", std::ofstream::out);

  // Grad in x
	x->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inY );
	x->getFieldPtr(Pimpact::V)->initField( Pimpact::ConstField, 0. );
	//x->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inX );
	//x->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inY );
	//x->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inX );
	//x->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inY );
  auto sol = x->clone( Pimpact::DeepCopy );
  x->write(0);
  //sol->write(3);
	
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

   for( int i=0; i<40; ++i ) {
     mg->apply( *b, *x );
     //x->write(i+10);

     temp->add( -1, *x, 1., *sol );
     S res = temp->norm();
     std::cout << "res: " << res << "\n";

     if( space()->rankST()==0 )
       ofs << res << "\n";
   }

   TEST_EQUALITY( temp->norm()<0.5, true );

   x->write(2);

   if( space()->rankST()==0 )
     ofs.close();

}


} // end of namespace
