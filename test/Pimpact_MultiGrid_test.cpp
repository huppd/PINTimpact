#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"

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
  //  // processor grid size
//    pl->set( "nx", 65 );
//    pl->set( "ny", 65 );
    pl->set( "nx", 1025 );
    pl->set( "ny", 1025 );

    pl->set("npx", 2 );
    pl->set("npy", 2 );

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
  pl->set( "domain", 1 );

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



TEUCHOS_UNIT_TEST( MGFields, constructor3D ) {

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



TEUCHOS_UNIT_TEST( MGOperators, constructor3D ) {

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



TEUCHOS_UNIT_TEST( MGSmoothers, constructor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;


  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 10 );


  auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradOp,Pimpact::DivGradO2Op>( mgSpaces );

  auto mgSmoother = Pimpact::createMGSmoothers<Pimpact::DivGradO2JSmoother>( mgOps );

  auto op = mgSmoother->get( 2 );
  op->print();

}



TEUCHOS_UNIT_TEST( MultiGrid, Restrictor3D ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;

  auto space = Pimpact::createSpace<S,O,3>( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 2 );

  auto mgTransfers = Pimpact::createMGTransfers( mgSpaces );

  Pimpact::EField type[] = {Pimpact::EField::S, Pimpact::EField::U, Pimpact::EField::V };

  for( int i=0; i<3; ++i ) {
    std::cout << "type: " << i << "\n";

    auto fieldf = Pimpact::createScalarField( mgSpaces->get( 0 ), type[i] );
    auto fieldc = Pimpact::createScalarField( mgSpaces->get( 1 ), type[i] );

    auto op = mgTransfers->getRestrictionOp( 0 );

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
    fieldf->initField( Pimpact::ConstField,1. );

    TEST_FLOATING_EQUALITY( 1., fieldf->norm(Belos::InfNorm), eps );

    TEST_FLOATING_EQUALITY( (S)fieldf->getLength(), fieldf->norm(Belos::OneNorm), eps  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldf->getLength() ), fieldf->norm(Belos::TwoNorm), eps  );

    fieldc->init(0.);

    op->apply( *fieldf, *fieldc );

    fieldf->write( 0 );
    fieldc->write( 1 );

    TEST_FLOATING_EQUALITY( 1., fieldc->norm(Belos::InfNorm), eps );

    TEST_FLOATING_EQUALITY( (S)fieldc->getLength(), fieldc->norm(Belos::OneNorm), eps  );

    TEST_FLOATING_EQUALITY( std::sqrt( (S)fieldc->getLength() ), fieldc->norm(Belos::TwoNorm), eps  );

    // the hard test
    fieldf->initField( Pimpact::Grad2D_inX, 1. );

    fieldc->init(0.);
    auto sol = fieldc->clone();
    sol->initField( Pimpact::Grad2D_inX, 1. );


    op->apply( *fieldf, *fieldc );

    fieldf->write( 0 );
    fieldc->write( 1 );

    sol->add( 1., *sol, -1., *fieldc );
    sol->write(3);

//    TEST_FLOATING_EQUALITY( 0., sol->norm(), eps ); // boundaries?

  }

}



/// \todo remove corners for test
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

  auto mgTransfers = Pimpact::createMGTransfers( mgSpaces );

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

    TEST_FLOATING_EQUALITY( 0., fieldf->norm(), eps );
    TEST_FLOATING_EQUALITY( 0., fieldc->norm(), eps );

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


template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;


TEUCHOS_UNIT_TEST( MultiGrid, DivGradOp ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;


  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 5 );


  auto mg =
      Pimpact::createMultiGrid<
        Pimpact::ScalarField,
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
   x->random();
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
     TEST_EQUALITY( sol->norm()<0.5, true );
   }

   x->write(2);

   if( space()->rankST()==0 )
     ofs.close();

   auto xm = Pimpact::createMultiField( x );
   auto bm = Pimpact::createMultiField( b );

   xm->init(0.);

   auto opc = Pimpact::create<Pimpact::DivGradOp>( space );

   auto prec = Pimpact::createMultiOperatorBase( mg );

   auto param = Pimpact::createLinSolverParameter("GMRES",1.e-9);
   auto bop = Pimpact::createMultiOperatorBase( opc );

   auto linprob = Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::ScalarField<FSpace3T> > >( bop, xm, bm, param, "GMRES" );
   linprob->setRightPrec(prec);

   linprob->solve(xm,bm);
   xm->write();
//   auto bop = Pimpact::createMultiOperatorBase( mg );

}



template<class ST> using ConvDiffOp = Pimpact::ConvectionVOp< Pimpact::ConvectionVWrap< Pimpact::ConvectionDiffusionSOp<ST> > >;


TEUCHOS_UNIT_TEST( MultiGrid, ConvDiffOp ) {

  typedef Pimpact::CoarsenStrategy<FSpace3T,CSpace3T> CS;


  auto space = Pimpact::createSpace( pl );

  auto mgSpaces = Pimpact::createMGSpaces<FSpace3T,CSpace3T,CS>( space, 5 );


//  auto mg =
//      Pimpact::createMultiGrid<
//        Pimpact::ScalarField,
//        ConvDiffOp,
//        ConvDiffOp,
//        Pimpact::DivGradO2JSmoother,
//        MOP>( mgSpaces );

  auto x = Pimpact::create<Pimpact::ScalarField>( space );
  auto b = Pimpact::create<Pimpact::ScalarField>( space );
//  auto op = Pimpact::create< ConvDiffOp<FSpace3T> >( space );

//  std::ofstream ofs;
//  if( space()->rankST()==0 )
//    ofs.open("MG.txt", std::ofstream::out);
//
//  // Grad in x
//   x->initField( Pimpact::Grad2D_inX );
//   x->write(0);
//
//   op->apply(*x,*b);
//   b->write(1);
//
//   x->init( 0. );
//   x->random();
//   auto e = x->clone();
//
//   e->init( 1. );
//   auto sol = x->clone();
////   auto  = x->clone();
//   sol->initField( Pimpact::Grad2D_inX );
//
//   S bla = b->dot(*e)/x->getLength();
//   std::cout<< " rhs nullspace: " << bla << "\n";
//   for( int i=0; i<40; ++i ) {
//     mg->apply( *b, *x );
//     x->level();
//     x->write(i+10);
//
//     op->apply(*x,*sol);
//     sol->add( -1, *b, 1., *sol );
//     double res = sol->norm();
//     std::cout << "res: " << res << "\n";
//
//     if( space()->rankST()==0 )
//       ofs << res << "\n";
//     TEST_EQUALITY( sol->norm()<0.5, true );
//   }
//
//   x->write(2);
//
//   if( space()->rankST()==0 )
//     ofs.close();
//
//   auto xm = Pimpact::createMultiField( x );
//   auto bm = Pimpact::createMultiField( b );
//
//   xm->init(0.);
//
//   auto opc = Pimpact::create<Pimpact::DivGradOp>( space );
//
//   auto prec = Pimpact::createMultiOperatorBase( mg );
//
//   auto param = Pimpact::createLinSolverParameter("GMRES",1.e-9);
//   auto bop = Pimpact::createMultiOperatorBase( opc );
//
//   auto linprob = Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::ScalarField<FSpace3T> > >( bop, xm, bm, param, "GMRES" );
//   linprob->setRightPrec(prec);
//
//   linprob->solve(xm,bm);
//   xm->write();
////   auto bop = Pimpact::createMultiOperatorBase( mg );
}


} // end of namespace
