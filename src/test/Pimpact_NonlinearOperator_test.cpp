#include <iostream>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Fields.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_LinSolverParameter.hpp"

#include "Pimpact_Test.hpp"



namespace {


const int sd = 3;

using SpaceT = Pimpact::Space<ST,OT,sd,d,dNC>;

using SF = typename Pimpact::ScalarField<SpaceT>;
using VF = typename Pimpact::VectorField<SpaceT>;
using MSF = typename Pimpact::ModeField<SF>;
using MVF = typename Pimpact::ModeField<VF>;

template<class T>
using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;
template<class T>
using ConvDiffSORT = Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionSORSmoother >;



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, ConvectionSOp, SpaceT ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

	auto u =
		Teuchos::tuple(
				Teuchos::tuple(
					Pimpact::createScalarField( space, Pimpact::U ),
					Pimpact::createScalarField( space, Pimpact::U ),
					Pimpact::createScalarField( space, Pimpact::U )
					),
				Teuchos::tuple(
					Pimpact::createScalarField( space, Pimpact::V ),
					Pimpact::createScalarField( space, Pimpact::V ),
					Pimpact::createScalarField( space, Pimpact::V )
					),
				Teuchos::tuple(
					Pimpact::createScalarField( space, Pimpact::W ),
					Pimpact::createScalarField( space, Pimpact::W ),
					Pimpact::createScalarField( space, Pimpact::W )
					)
				);

	auto x    = Pimpact::create<Pimpact::VectorField>( space );
	auto y    = Pimpact::create<Pimpact::VectorField>( space );
	auto solv = Pimpact::create<Pimpact::VectorField>( space );

	auto op = Pimpact::create<Pimpact::ConvectionSOp>( space ) ;

	if( print ) op->print();

	std::cout << "\n";
	for( int dir=-1; dir<2; dir+= 2 ) {

		for( int i=0; i<SpaceT::sdim; ++i )
			for( int j=0; j<SpaceT::sdim; ++j )
				u[i][j]->initField( Pimpact::ConstField, 2.*dir );

		// dx test
		solv->getField(Pimpact::U).initField( Pimpact::ConstField, 2.*dir );
		solv->getField(Pimpact::V).initField( Pimpact::ConstField, 2.*dir );
		if( 3==SpaceT::sdim )
			solv->getField(Pimpact::W).initField( Pimpact::ConstField, 2.*dir );

		x->getField(Pimpact::U).initField( Pimpact::Grad2D_inX );
		x->getField(Pimpact::V).initField( Pimpact::Grad2D_inX );
		if( 3==SpaceT::sdim )
			x->getField(Pimpact::W).initField( Pimpact::Grad2D_inX );

		y->random();

		for( int i=0; i<SpaceT::sdim; ++i ) {
			op->apply( u[i], x->getConstField(i), y->getField(i) );
		}

		for( int i=0; i<SpaceT::sdim; ++i ) {
			auto sol = solv->getField( i );
			sol.add( 1., sol, -1., y->getField(i) );
			ST errorInf = sol.norm( Belos::InfNorm );
			std::cout << "error in "<< Pimpact::toString( static_cast<Pimpact::EField>(i) ) << " (gradX): " << errorInf << "\n";
			TEST_EQUALITY( errorInf<eps, true );
			if( write && errorInf>=eps )
				sol.write();
		}

		// dy test
		solv->getField(Pimpact::U).initField( Pimpact::ConstField, 2.*dir );
		solv->getField(Pimpact::V).initField( Pimpact::ConstField, 2.*dir );
		if( 3==SpaceT::sdim )
			solv->getField(Pimpact::W).initField( Pimpact::ConstField, 2.*dir );

		x->getField(Pimpact::U).initField( Pimpact::Grad2D_inY );
		x->getField(Pimpact::V).initField( Pimpact::Grad2D_inY );
		if( 3==SpaceT::sdim )
			x->getField(Pimpact::W).initField( Pimpact::Grad2D_inY );

		y->random();

		for( int i=0; i<SpaceT::sdim; ++i ) {
			op->apply( u[i], x->getConstField(i), y->getField(i) );
		}

		for( int i=0; i<SpaceT::sdim; ++i ){
			auto sol = solv->getField( i );
			sol.add( 1., sol, -1., y->getField(i) );
			std::cout << "error in "<< Pimpact::toString( static_cast<Pimpact::EField>(i) ) << " (gradY): " << sol.norm( Belos::InfNorm ) << "\n";
			TEST_EQUALITY( sol.norm( Belos::InfNorm )<eps, true );
		}
		x->getField(Pimpact::U).initField( Pimpact::Grad2D_inY );
		x->getField(Pimpact::V).initField( Pimpact::Grad2D_inY );
		if( 3==SpaceT::sdim )
			x->getField(Pimpact::W).initField( Pimpact::Grad2D_inY );

		// dz test
		if( 3==SpaceT::sdim ) {
			solv->getField(Pimpact::U).initField( Pimpact::ConstField, 2.*dir );
			solv->getField(Pimpact::V).initField( Pimpact::ConstField, 2.*dir );
			solv->getField(Pimpact::W).initField( Pimpact::ConstField, 2.*dir );

			x->getField(Pimpact::U).initField( Pimpact::Grad2D_inZ );
			x->getField(Pimpact::V).initField( Pimpact::Grad2D_inZ );
			x->getField(Pimpact::W).initField( Pimpact::Grad2D_inZ );

			y->random();

			for( int i=0; i<SpaceT::sdim; ++i ) {
				op->apply( u[i], x->getConstField(i), y->getField(i) );
			}

			for( int i=0; i<SpaceT::sdim; ++i ) {
				auto sol = solv->getField( i );
				sol.add( 1., sol, -1., y->getField(i) );
				ST errorInf = sol.norm( Belos::InfNorm );
				std::cout << "error in "<< Pimpact::toString( static_cast<Pimpact::EField>(i) ) << " (gradZ): " << errorInf << "\n";
				TEST_EQUALITY( errorInf<eps, true );
				if( write && errorInf>=eps )
					//sol->write();
					sol.print();
				//y->getField(i).print();
			}
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionSOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionSOp, D3 )



template<class T> using ConvOpT = Pimpact::NonlinearOp<Pimpact::ConvectionSOp<T> >;


TEUCHOS_UNIT_TEST( BasicOperator, NonlinearOp ) {

	setParameter( SpaceT::sdim );

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

  auto x = Pimpact::create<Pimpact::VectorField>( space );
  auto y = Pimpact::create<Pimpact::VectorField>( space );
  auto z = Pimpact::create<Pimpact::VectorField>( space );
  auto z2 = Pimpact::create<Pimpact::VectorField>( space );


  auto op = Pimpact::create<ConvOpT>( space  ) ;

	for(int dir=-1; dir<2; dir+=2 ) {
		// x test
		x->getField( Pimpact::U ).initField( Pimpact::ConstField, 2.*dir );
		x->getField( Pimpact::V ).initField( Pimpact::ConstField, 2.*dir );
		x->getField( Pimpact::W ).initField( Pimpact::ConstField, 2.*dir );
		y->getField( Pimpact::U ).initField( Pimpact::Grad2D_inX );
		y->getField( Pimpact::V ).initField( Pimpact::Grad2D_inX );
		y->getField( Pimpact::W ).initField( Pimpact::Grad2D_inX );
		z->random();
		z2->getField( Pimpact::U ).initField( Pimpact::ConstField, 2.*dir );
		z2->getField( Pimpact::V ).initField( Pimpact::ConstField, 2.*dir );
		z2->getField( Pimpact::W ).initField( Pimpact::ConstField, 2.*dir );

		op->assignField( *x );
		op->apply( *y, *z );

		z2->add( -1, *z2, 1, *z );

		std::cout << "error in "<< dir << " (gradX): " << z2->norm( Belos::InfNorm ) << "\n";
		TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

		// y test
		x->getField( Pimpact::U ).initField( Pimpact::ConstField, 2.*dir );
		x->getField( Pimpact::V ).initField( Pimpact::ConstField, 2.*dir );
		x->getField( Pimpact::W ).initField( Pimpact::ConstField, 2.*dir );
		y->getField( Pimpact::U ).initField( Pimpact::Grad2D_inY );
		y->getField( Pimpact::V ).initField( Pimpact::Grad2D_inY );
		y->getField( Pimpact::W ).initField( Pimpact::Grad2D_inY );
		z->random();
		z2->getField( Pimpact::U ).initField( Pimpact::ConstField, 2.*dir );
		z2->getField( Pimpact::V ).initField( Pimpact::ConstField, 2.*dir );
		z2->getField( Pimpact::W ).initField( Pimpact::ConstField, 2.*dir );

		op->assignField( *x );
		op->apply( *y, *z );

		z2->add( -1, *z2, 1, *z );

		std::cout << "error in "<< dir << " (gradY): " << z2->norm( Belos::InfNorm ) << "\n";
		TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

		// z test
		x->getField( Pimpact::U ).initField( Pimpact::ConstField, 2.*dir );
		x->getField( Pimpact::V ).initField( Pimpact::ConstField, 2.*dir );
		x->getField( Pimpact::W ).initField( Pimpact::ConstField, 2.*dir );
		y->getField( Pimpact::U ).initField( Pimpact::Grad2D_inZ );
		y->getField( Pimpact::V ).initField( Pimpact::Grad2D_inZ );
		y->getField( Pimpact::W ).initField( Pimpact::Grad2D_inZ );
		z->random();
		z2->getField( Pimpact::U ).initField( Pimpact::ConstField, 2.*dir );
		z2->getField( Pimpact::V ).initField( Pimpact::ConstField, 2.*dir );
		z2->getField( Pimpact::W ).initField( Pimpact::ConstField, 2.*dir );

		op->assignField( *x );
		op->apply( *y, *z );

		z2->add( -1, *z2, 1, *z );

		std::cout << "error in "<< dir << " (gradZ): " << z2->norm( Belos::InfNorm ) << "\n";
		TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );
	}

	// x test
	x->initField();
	x->getField( Pimpact::U ).initField( Pimpact::Poiseuille2D_inX );
	y->getField( Pimpact::U ).initField( Pimpact::Grad2D_inX );
	y->getField( Pimpact::V ).initField( Pimpact::Grad2D_inY );
	y->getField( Pimpact::W ).initField( Pimpact::Grad2D_inZ );
	z->random();
	z2->initField();
	z2->getField( Pimpact::U ).initField( Pimpact::Poiseuille2D_inX );

	op->assignField( *x );
	op->apply( *y, *z );

	z2->add( -1, *z2, 1, *z );

	std::cout << "error in (poisX): " << z2->norm( Belos::InfNorm ) << "\n";
	TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

	// y test
	x->initField();
	x->getField( Pimpact::V ).initField( Pimpact::Poiseuille2D_inY );
	y->getField( Pimpact::U ).initField( Pimpact::Grad2D_inX );
	y->getField( Pimpact::V ).initField( Pimpact::Grad2D_inY );
	y->getField( Pimpact::W ).initField( Pimpact::Grad2D_inZ );
	z->random();
	z2->initField();
	z2->getField( Pimpact::V ).initField( Pimpact::Poiseuille2D_inY );

	op->assignField( *x );
	op->apply( *y, *z );

	z2->add( -1, *z2, 1, *z );

	std::cout << "error in (poisY): " << z2->norm( Belos::InfNorm ) << "\n";
	TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

	// z test
	x->initField();
	x->getField( Pimpact::W ).initField( Pimpact::Poiseuille2D_inZ );
	y->getField( Pimpact::U ).initField( Pimpact::Grad2D_inX );
	y->getField( Pimpact::V ).initField( Pimpact::Grad2D_inY );
	y->getField( Pimpact::W ).initField( Pimpact::Grad2D_inZ );
	z->random();
	z2->initField();
	z2->getField( Pimpact::W ).initField( Pimpact::Poiseuille2D_inZ );

	op->assignField( *x );
	op->apply( *y, *z );

	z2->add( -1, *z2, 1, *z );

	std::cout << "error in (poisZ): " << z2->norm( Belos::InfNorm ) << "\n";
	TEST_EQUALITY( z2->norm( Belos::InfNorm, Pimpact::With::noB  )<eps, true );

}



TEUCHOS_UNIT_TEST( BasicOperator, ConvectionDiffusionOp  ) {

	setParameter( SpaceT::sdim );

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

  auto x = Pimpact::create<Pimpact::VectorField>( space );
  auto y = Pimpact::create<Pimpact::VectorField>( space );
  auto z = Pimpact::create<Pimpact::VectorField>( space );
  auto z2 = Pimpact::create<Pimpact::VectorField>( space );
  auto sol = y->clone(Pimpact::ECopy::Shallow);


  auto op = Pimpact::createAdd2Op(
      Pimpact::create< ConvOpT<SpaceT> >(space),
      Pimpact::create<Pimpact::HelmholtzOp>(space) );

	auto op2 = Pimpact::create< ConvDiffOpT<SpaceT> >( space );


  x->getField( Pimpact::U ).initField( Pimpact::ConstField, 2. );
  x->getField( Pimpact::V ).initField( Pimpact::ConstField, 2. );
  x->getField( Pimpact::W ).initField( Pimpact::ConstField, 2. );
  sol->init( Teuchos::tuple(2.,2.,2.) );

  // consistency test in x
  y->getField(Pimpact::U).initField( Pimpact::Grad2D_inX );
  y->getField(Pimpact::V).initField( Pimpact::Grad2D_inX );
  y->getField(Pimpact::W).initField( Pimpact::Grad2D_inX );

  op->assignField( *x );
  op2->assignField( *x );

  // 
  op->apply(   *y, *z  );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

	std::cout << "diff (gradX): " << z2->norm( Belos::InfNorm, Pimpact::With::noB  ) << "\n";
	std::cout << "erro (gradX): " << sol->norm( Belos::InfNorm, Pimpact::With::noB  ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm, Pimpact::With::noB  ) < eps, true );
  TEST_EQUALITY( sol->norm( Belos::InfNorm, Pimpact::With::noB  ) < eps, true );

  // consistency test in y
  y->getField(Pimpact::U).initField( Pimpact::Grad2D_inY );
  y->getField(Pimpact::V).initField( Pimpact::Grad2D_inY );
  y->getField(Pimpact::W).initField( Pimpact::Grad2D_inY );
  sol->init( Teuchos::tuple(2.,2.,2.) );

  op->assignField( *x );
  op2->assignField( *x );

  op->apply(   *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

	std::cout << "diff (gradY): " << z2->norm( Belos::InfNorm, Pimpact::With::noB  ) << "\n";
	std::cout << "erro (gradY): " << sol->norm( Belos::InfNorm, Pimpact::With::noB  ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm, Pimpact::With::noB  ) < eps, true );
  TEST_EQUALITY( sol->norm( Belos::InfNorm, Pimpact::With::noB  ) < eps, true );

  // consistency test in z
  y->getField(Pimpact::U).initField( Pimpact::Grad2D_inZ );
  y->getField(Pimpact::V).initField( Pimpact::Grad2D_inZ );
  y->getField(Pimpact::W).initField( Pimpact::Grad2D_inZ );
  sol->init( Teuchos::tuple(2.,2.,2.) );

  op->assignField( *x );
  op2->assignField( *x );

  op->apply( *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

	std::cout << "diff (gradZ): " << z2->norm( Belos::InfNorm, Pimpact::With::noB  ) << "\n";
	std::cout << "erro (gradZ): " << sol->norm( Belos::InfNorm, Pimpact::With::noB  ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm, Pimpact::With::noB  ) < eps, true );
  TEST_EQUALITY( sol->norm( Belos::InfNorm, Pimpact::With::noB  ) < eps, true );


  // consistency test in pois x
	y->initField();
  y->getField(Pimpact::U).initField( Pimpact::Poiseuille2D_inX );
  sol->init( Teuchos::tuple(0.,0.,0.) );
  z->random();

  op->apply( *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

	std::cout << "diff (poisX): " << z2->norm( Belos::InfNorm, Pimpact::With::noB  ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm, Pimpact::With::noB  )<eps, true );


  // consistency test in pois y
	y->initField();
  y->getField(Pimpact::V).initField( Pimpact::Poiseuille2D_inY );
  sol->init( Teuchos::tuple(0.,0.,0.) );
  z->random();

  op->apply( *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );

	std::cout << "diff (poisY): " << z2->norm( Belos::InfNorm, Pimpact::With::noB  ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm, Pimpact::With::noB  )<eps, true );

  // consistency test in pois Z
	y->initField();
  y->getField(Pimpact::W).initField( Pimpact::Poiseuille2D_inZ );
  z->random();

  op->apply( *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );

	std::cout << "diff (poisZ): " << z2->norm( Belos::InfNorm, Pimpact::With::noB  ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm, Pimpact::With::noB  )<eps, true );

}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionDiffusionSORSmoother ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<Pimpact::Space<ST,OT,sd,d,2> >( pl );

  auto wind = Pimpact::create<Pimpact::VectorField>( space );
  auto x = Pimpact::create<Pimpact::VectorField>( space );
  auto y = Pimpact::create<Pimpact::VectorField>( space );
  auto z = Pimpact::create<Pimpact::VectorField>( space );
  auto z2 = Pimpact::create<Pimpact::VectorField>( space );


  auto op = Pimpact::create<ConvDiffOpT>( space );

  auto pls = Teuchos::parameterList();
  pls->set( "omega", 1. );
  pls->set( "numIter", 1 );
  pls->set( "Ordering", 0 );
  pls->set<short int>( "dir X",  1 );
  pls->set<short int>( "dir Y",  1 );
  pls->set<short int>( "dir Z",  1 );

  auto smoother = Pimpact::create< ConvDiffSORT >( op, pls );

	if( print ) smoother->print();

  wind->getField(Pimpact::U).initField( Pimpact::ConstField, winds );
  wind->getField(Pimpact::V).initField( Pimpact::ConstField, winds );
  wind->getField(Pimpact::W).initField( Pimpact::ConstField, winds );

	// Consistency
	x->random();
	auto xs = x->clone( Pimpact::ECopy::Deep ); 

	op->apply( *x, *y );
	smoother->apply( *y, *x );

	if( write ) x->write(1);
	if( write ) xs->write(2);

	x->add( 1., *xs, -1., *x );

	ST error0 = x->norm( Belos::InfNorm );

	if( space()->rankST()==0 )
		std::cout << "\nConsistency error: " << error0 << "\n";
	if( write ) x->write( 0 );
	if( print ) x->print();

	// X test
	ST n;
  z->initField();

  y->getField( Pimpact::U ).initField( Pimpact::Grad2D_inX );
  y->getField( Pimpact::V ).initField( Pimpact::Grad2D_inX );
  y->getField( Pimpact::W ).initField( Pimpact::Grad2D_inX );

  auto sol = y->clone( Pimpact::ECopy::Deep );

  op->assignField( *wind );

  op->apply( *y, *z );

  {
    y->init(0);
    auto bc = z->clone( Pimpact::ECopy::Shallow );
    op->apply( *y, *bc );
    z->add( 1., *z, -1., *bc );
  }


  y->initField();

  for(int i=0; i<10; ++i) {
    smoother->apply( *z, *y );
    z2->add( -1, *sol, 1, *y );
    n = z2->norm( Belos::InfNorm );
    if( space()->rankST()==0 )
      std::cout << "error: " << n << "\n";
  }

  TEST_EQUALITY_CONST( z2->norm( Belos::InfNorm )/z2->getLength()<eps, true );

	// Y test
  z->initField();

  y->getField( Pimpact::U ).initField( Pimpact::Grad2D_inY );
  y->getField( Pimpact::V ).initField( Pimpact::Grad2D_inY );
  y->getField( Pimpact::W ).initField( Pimpact::Grad2D_inY );

  sol = y->clone( Pimpact::ECopy::Deep );

  op->assignField( *wind );

  op->apply( *y, *z );

  {
    y->init(0);
    auto bc = z->clone( Pimpact::ECopy::Shallow );
    op->apply( *y, *bc );
    z->add( 1., *z, -1., *bc );
  }

  y->initField();

  for(int i=0; i<20; ++i) {
    smoother->apply( *z, *y );

    z2->add( -1, *sol, 1, *y );

    n = z2->norm( Belos::InfNorm, Pimpact::With::noB  );
    if( space()->rankST()==0 )
      std::cout << "error: " << n << "\n";

  }

  TEST_EQUALITY_CONST( z2->norm( Belos::InfNorm, Pimpact::With::noB  )/z2->getLength()<eps, true );

	// Z test
  z->initField();

  y->getField( Pimpact::U ).initField( Pimpact::Grad2D_inZ );
  y->getField( Pimpact::V ).initField( Pimpact::Grad2D_inZ );
  y->getField( Pimpact::W ).initField( Pimpact::Grad2D_inZ );

  sol = y->clone( Pimpact::ECopy::Deep );

  op->assignField( *wind );

  op->apply( *y, *z );

  {
    y->init(0);
    auto bc = z->clone( Pimpact::ECopy::Shallow );
    op->apply( *y, *bc );
    z->add( 1., *z, -1., *bc );
  }

  y->initField();

  for(int i=0; i<20; ++i) {
    smoother->apply( *z, *y );

    z2->add( -1, *sol, 1, *y );

    n = z2->norm( Belos::InfNorm, Pimpact::With::noB  );
    if( space()->rankST()==0 )
      std::cout << "error: " << n << "\n";

  }

  TEST_EQUALITY_CONST( z2->norm( Belos::InfNorm, Pimpact::With::noB  )/z2->getLength()<eps, true );

}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionDiffusionJSmoother ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<Pimpact::Space<ST,OT,sd,d,2> >( pl );

  auto op = Pimpact::create<ConvDiffOpT>( space );

	// init smoother
  auto pls = Teuchos::parameterList();
	pls->set( "omega", 1. );
  pls->set( "numIters", 1 );

  auto smoother =
      Pimpact::create<
        Pimpact::NonlinearSmoother<
          ConvDiffOpT<Pimpact::Space<ST,OT,sd,d,2> > ,
          Pimpact::ConvectionDiffusionJSmoother > > (
              op,
              pls );

	//  init wind
	{
		auto wind = Pimpact::create<Pimpact::VectorField>( space );
		wind->getField(Pimpact::U).initField( Pimpact::ConstField, 1. );
		wind->getField(Pimpact::V).initField( Pimpact::ConstField, 1. );
		wind->getField(Pimpact::W).initField( Pimpact::ConstField, 1. );

		op->assignField( *wind );
	}

	// init initial guess
  auto x = Pimpact::create<Pimpact::VectorField>( space );
	x->random();
	//x->init();

	// init rhs
  auto y = Pimpact::create<Pimpact::VectorField>( space );

	// Consistency
	auto xs = x->clone( Pimpact::ECopy::Deep ); 
	op->apply( *x, *y );

	smoother->apply( *y, *x );
	x->write(1);
	xs->write(2);
	x->add( 1., *xs, -1., *x );
	//x->abs( *x );
	ST error0 = x->norm( Belos::InfNorm );
	if( space()->rankST()==0 )
		std::cout << "\nConsistency error: " << error0 << "\n";
  TEST_EQUALITY_CONST( error0<eps, true );
	if( write ) x->write( 0 );
	if( print ) x->print();

	// Convergence
	y->init( 0. );
	x->random();
	error0 = x->norm( Belos::InfNorm );
	ST error;
	if( space()->rankST()==0 )
		std::cout << "\nerror: " << error0 << "\n";

	//x->write( 0 );
  for(int i=0; i<10; ++i) {
    smoother->apply( *y, *x );

    error = x->norm( Belos::InfNorm )/error0;
    if( space()->rankST()==0 )
      std::cout << "error: " << error << "\n";

		if( write ) x->write( i+1 );
  }

  TEST_EQUALITY_CONST( error<0.1, true );
}




TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicDtConvectionDiffusionOp ) {

	pl->set<bool>( "spectral in time", true );

	ST pi2 = 2.*std::acos(-1.);

	lx = pi2 ;
	ly = pi2 ;
	lz = pi2 ;

	nf = 4;

	setParameter( SpaceT::sdim );
	Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );


	auto x   = Pimpact::createMultiHarmonicVectorField( space );
	auto y   = Pimpact::createMultiHarmonicVectorField( space );
	auto sol = Pimpact::createMultiHarmonicVectorField( space );
	auto err = Pimpact::createMultiHarmonicVectorField( space );

	auto op = Pimpact::createMultiDtConvectionDiffusionOp( space );

	// initializtion

	x->getSField(1).getField(Pimpact::U).initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST {  return(  std::cos(x*pi2)*std::sin(y*pi2) ); } );
	x->getSField(1).getField(Pimpact::V).initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST {  return( -std::sin(x*pi2)*std::cos(y*pi2) ); } );

	if( write ) x->write( 10 );

	// solution init
	
	sol->get0Field().getField(Pimpact::U).initFromFunction(
			[&pi2,&re]( ST x, ST y, ST z ) ->ST {  return(  -std::sin(2.*x*pi2)/4. ); } );
	sol->get0Field().getField(Pimpact::V).initFromFunction(
			[&pi2,&re]( ST x, ST y, ST z ) ->ST {  return(  -std::sin(2.*y*pi2)/4. ); } );

	sol->getCField(1).getField(Pimpact::U).initFromFunction(
	    	[&pi2,&alpha2,&re]( ST x, ST y, ST z ) ->ST {  return(  alpha2*std::cos(x*pi2)*std::sin(y*pi2)/re ); } );
	sol->getCField(1).getField(Pimpact::V).initFromFunction(
	    	[&pi2,&alpha2,&re]( ST x, ST y, ST z ) ->ST {  return( -alpha2*std::sin(x*pi2)*std::cos(y*pi2)/re ); } );

	sol->getSField(1).getField(Pimpact::U).initFromFunction(
				[&pi2,&re]( ST x, ST y, ST z ) ->ST {  return(  2.*std::cos(x*pi2)*std::sin(y*pi2)/re ); } );
	sol->getSField(1).getField(Pimpact::V).initFromFunction(
				[&pi2,&re]( ST x, ST y, ST z ) ->ST {  return( -2.*std::sin(x*pi2)*std::cos(y*pi2)/re ); } );

	sol->getCField(2).getField(Pimpact::U).initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST { return( std::sin(2.*x*pi2)/4. ); } );
	sol->getCField(2).getField(Pimpact::V).initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST { return( std::sin(2.*y*pi2)/4. ); } );

	if( write ) sol->write( 30 );

	op->assignField( *x );
	op->apply( *x, *y );

	if( write ) y->write( 20 );

	err->add( 1., *sol, -1., *y );
	if( write ) err->write( 0 );

	ST error = err->norm()/sol->norm();
	std::cout << "\nerror: " << error << "\n";

}





} // namespace
