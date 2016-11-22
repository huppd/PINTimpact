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
		solv->getFieldPtr(Pimpact::U)->initField( Pimpact::ConstField, 2.*dir );
		solv->getFieldPtr(Pimpact::V)->initField( Pimpact::ConstField, 2.*dir );
		if( 3==SpaceT::sdim )
			solv->getFieldPtr(Pimpact::W)->initField( Pimpact::ConstField, 2.*dir );

		x->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inX );
		x->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inX );
		if( 3==SpaceT::sdim )
			x->getFieldPtr(Pimpact::W)->initField( Pimpact::Grad2D_inX );

		y->random();

		for( int i=0; i<SpaceT::sdim; ++i ) {
			op->apply( u[i], x->getConstField(i), y->getField(i) );
		}

		for( int i=0; i<SpaceT::sdim; ++i ) {
			auto sol = solv->getFieldPtr( i );
			sol->add( 1., *sol, -1., y->getField(i) );
			ST errorInf = sol->norm( Belos::InfNorm );
			std::cout << "error in "<< Pimpact::toString( static_cast<Pimpact::EField>(i) ) << " (gradX): " << errorInf << "\n";
			TEST_EQUALITY( errorInf<eps, true );
			if( write && errorInf>=eps )
				sol->write();
		}

		// dy test
		solv->getFieldPtr(Pimpact::U)->initField( Pimpact::ConstField, 2.*dir );
		solv->getFieldPtr(Pimpact::V)->initField( Pimpact::ConstField, 2.*dir );
		if( 3==SpaceT::sdim )
			solv->getFieldPtr(Pimpact::W)->initField( Pimpact::ConstField, 2.*dir );

		x->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inY );
		x->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inY );
		if( 3==SpaceT::sdim )
			x->getFieldPtr(Pimpact::W)->initField( Pimpact::Grad2D_inY );

		y->random();

		for( int i=0; i<SpaceT::sdim; ++i ) {
			op->apply( u[i], x->getConstField(i), y->getField(i) );
		}

		for( int i=0; i<SpaceT::sdim; ++i ){
			auto sol = solv->getFieldPtr( i );
			sol->add( 1., *sol, -1., y->getField(i) );
			std::cout << "error in "<< Pimpact::toString( static_cast<Pimpact::EField>(i) ) << " (gradY): " << sol->norm( Belos::InfNorm ) << "\n";
			TEST_EQUALITY( sol->norm( Belos::InfNorm )<eps, true );
		}
		x->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inY );
		x->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inY );
		if( 3==SpaceT::sdim )
			x->getFieldPtr(Pimpact::W)->initField( Pimpact::Grad2D_inY );

		// dz test
		if( 3==SpaceT::sdim ) {
			solv->getFieldPtr(Pimpact::U)->initField( Pimpact::ConstField, 2.*dir );
			solv->getFieldPtr(Pimpact::V)->initField( Pimpact::ConstField, 2.*dir );
			solv->getFieldPtr(Pimpact::W)->initField( Pimpact::ConstField, 2.*dir );

			x->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inZ );
			x->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inZ );
			x->getFieldPtr(Pimpact::W)->initField( Pimpact::Grad2D_inZ );

			y->random();

			for( int i=0; i<SpaceT::sdim; ++i ) {
				op->apply( u[i], x->getConstField(i), y->getField(i) );
			}

			for( int i=0; i<SpaceT::sdim; ++i ){
				auto sol = solv->getFieldPtr( i );
				sol->add( 1., *sol, -1., y->getField(i) );
				ST errorInf = sol->norm( Belos::InfNorm );
				std::cout << "error in "<< Pimpact::toString( static_cast<Pimpact::EField>(i) ) << " (gradZ): " << errorInf << "\n";
				TEST_EQUALITY( errorInf<eps, true );
				if( write && errorInf>=eps )
					//sol->write();
					sol->print();
				//y->getFieldPtr(i)->print();
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
		x->getFieldPtr( Pimpact::U )->initField( Pimpact::ConstField, 2.*dir );
		x->getFieldPtr( Pimpact::V )->initField( Pimpact::ConstField, 2.*dir );
		x->getFieldPtr( Pimpact::W )->initField( Pimpact::ConstField, 2.*dir );
		y->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
		y->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX );
		y->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX );
		z->random();
		z2->getFieldPtr( Pimpact::U )->initField( Pimpact::ConstField, 2.*dir );
		z2->getFieldPtr( Pimpact::V )->initField( Pimpact::ConstField, 2.*dir );
		z2->getFieldPtr( Pimpact::W )->initField( Pimpact::ConstField, 2.*dir );

		op->assignField( *x );
		op->apply( *y, *z );

		z2->add( -1, *z2, 1, *z );

		std::cout << "error in "<< dir << " (gradX): " << z2->norm( Belos::InfNorm ) << "\n";
		TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

		// y test
		x->getFieldPtr( Pimpact::U )->initField( Pimpact::ConstField, 2.*dir );
		x->getFieldPtr( Pimpact::V )->initField( Pimpact::ConstField, 2.*dir );
		x->getFieldPtr( Pimpact::W )->initField( Pimpact::ConstField, 2.*dir );
		y->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY );
		y->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
		y->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inY );
		z->random();
		z2->getFieldPtr( Pimpact::U )->initField( Pimpact::ConstField, 2.*dir );
		z2->getFieldPtr( Pimpact::V )->initField( Pimpact::ConstField, 2.*dir );
		z2->getFieldPtr( Pimpact::W )->initField( Pimpact::ConstField, 2.*dir );

		op->assignField( *x );
		op->apply( *y, *z );

		z2->add( -1, *z2, 1, *z );

		std::cout << "error in "<< dir << " (gradY): " << z2->norm( Belos::InfNorm ) << "\n";
		TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

		// z test
		x->getFieldPtr( Pimpact::U )->initField( Pimpact::ConstField, 2.*dir );
		x->getFieldPtr( Pimpact::V )->initField( Pimpact::ConstField, 2.*dir );
		x->getFieldPtr( Pimpact::W )->initField( Pimpact::ConstField, 2.*dir );
		y->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inZ );
		y->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inZ );
		y->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );
		z->random();
		z2->getFieldPtr( Pimpact::U )->initField( Pimpact::ConstField, 2.*dir );
		z2->getFieldPtr( Pimpact::V )->initField( Pimpact::ConstField, 2.*dir );
		z2->getFieldPtr( Pimpact::W )->initField( Pimpact::ConstField, 2.*dir );

		op->assignField( *x );
		op->apply( *y, *z );

		z2->add( -1, *z2, 1, *z );

		std::cout << "error in "<< dir << " (gradZ): " << z2->norm( Belos::InfNorm ) << "\n";
		TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );
	}

	// x test
	x->initField();
	x->getFieldPtr( Pimpact::U )->initField( Pimpact::Poiseuille2D_inX );
	y->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
	y->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
	y->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );
	z->random();
	z2->initField();
	z2->getFieldPtr( Pimpact::U )->initField( Pimpact::Poiseuille2D_inX );

	op->assignField( *x );
	op->apply( *y, *z );

	z2->add( -1, *z2, 1, *z );

	std::cout << "error in (poisX): " << z2->norm( Belos::InfNorm ) << "\n";
	TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

	// y test
	x->initField();
	x->getFieldPtr( Pimpact::V )->initField( Pimpact::Poiseuille2D_inY );
	y->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
	y->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
	y->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );
	z->random();
	z2->initField();
	z2->getFieldPtr( Pimpact::V )->initField( Pimpact::Poiseuille2D_inY );

	op->assignField( *x );
	op->apply( *y, *z );

	z2->add( -1, *z2, 1, *z );

	std::cout << "error in (poisY): " << z2->norm( Belos::InfNorm ) << "\n";
	TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

	// z test
	x->initField();
	x->getFieldPtr( Pimpact::W )->initField( Pimpact::Poiseuille2D_inZ );
	y->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
	y->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
	y->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );
	z->random();
	z2->initField();
	z2->getFieldPtr( Pimpact::W )->initField( Pimpact::Poiseuille2D_inZ );

	op->assignField( *x );
	op->apply( *y, *z );

	z2->add( -1, *z2, 1, *z );

	std::cout << "error in (poisZ): " << z2->norm( Belos::InfNorm ) << "\n";
	TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

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


  x->getFieldPtr( Pimpact::U )->initField( Pimpact::ConstField, 2. );
  x->getFieldPtr( Pimpact::V )->initField( Pimpact::ConstField, 2. );
  x->getFieldPtr( Pimpact::W )->initField( Pimpact::ConstField, 2. );
  sol->init( Teuchos::tuple(2.,2.,2.) );

  // consistency test in x
  y->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inX );
  y->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inX );
  y->getFieldPtr(Pimpact::W)->initField( Pimpact::Grad2D_inX );

  op->assignField( *x );
  op2->assignField( *x );

  // 
  op->apply(   *y, *z  );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

	std::cout << "diff (gradX): " << z2->norm( Belos::InfNorm ) << "\n";
	std::cout << "erro (gradX): " << sol->norm( Belos::InfNorm ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm ) < eps, true );
  TEST_EQUALITY( sol->norm( Belos::InfNorm ) < eps, true );

  // consistency test in y
  y->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inY );
  y->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inY );
  y->getFieldPtr(Pimpact::W)->initField( Pimpact::Grad2D_inY );
  sol->init( Teuchos::tuple(2.,2.,2.) );

  op->assignField( *x );
  op2->assignField( *x );

  op->apply(   *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

	std::cout << "diff (gradY): " << z2->norm( Belos::InfNorm ) << "\n";
	std::cout << "erro (gradY): " << sol->norm( Belos::InfNorm ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm ) < eps, true );
  TEST_EQUALITY( sol->norm( Belos::InfNorm ) < eps, true );

  // consistency test in z
  y->getFieldPtr(Pimpact::U)->initField( Pimpact::Grad2D_inZ );
  y->getFieldPtr(Pimpact::V)->initField( Pimpact::Grad2D_inZ );
  y->getFieldPtr(Pimpact::W)->initField( Pimpact::Grad2D_inZ );
  sol->init( Teuchos::tuple(2.,2.,2.) );

  op->assignField( *x );
  op2->assignField( *x );

  op->apply( *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

	std::cout << "diff (gradZ): " << z2->norm( Belos::InfNorm ) << "\n";
	std::cout << "erro (gradZ): " << sol->norm( Belos::InfNorm ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm ) < eps, true );
  TEST_EQUALITY( sol->norm( Belos::InfNorm ) < eps, true );


  // consistency test in pois x
	y->initField();
  y->getFieldPtr(Pimpact::U)->initField( Pimpact::Poiseuille2D_inX );
  sol->init( Teuchos::tuple(0.,0.,0.) );
  z->random();

  op->apply( *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );
  sol->add(-1, *sol,1, *z );

	std::cout << "diff (poisX): " << z2->norm( Belos::InfNorm ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );


  // consistency test in pois y
	y->initField();
  y->getFieldPtr(Pimpact::V)->initField( Pimpact::Poiseuille2D_inY );
  sol->init( Teuchos::tuple(0.,0.,0.) );
  z->random();

  op->apply( *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );

	std::cout << "diff (poisY): " << z2->norm( Belos::InfNorm ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

  // consistency test in pois Z
	y->initField();
  y->getFieldPtr(Pimpact::W)->initField( Pimpact::Poiseuille2D_inZ );
  z->random();

  op->apply( *y, *z );
  op2->apply(  *y, *z2 );

  z2->add( -1, *z2, 1, *z );

	std::cout << "diff (poisZ): " << z2->norm( Belos::InfNorm ) << "\n";
  TEST_EQUALITY( z2->norm( Belos::InfNorm )<eps, true );

}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionDiffusionSORSmoother ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<Pimpact::Space<ST,OT,sd,d,2> >( pl );

  auto wind = Pimpact::create<Pimpact::VectorField>( space );
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

  wind->getFieldPtr(Pimpact::U)->initField( Pimpact::ConstField, winds );
  wind->getFieldPtr(Pimpact::V)->initField( Pimpact::ConstField, winds );
  wind->getFieldPtr(Pimpact::W)->initField( Pimpact::ConstField, winds );

	// X test
  z->initField();

  y->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
  y->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inX );
  y->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inX );

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

  double n;

  for(int i=0; i<20; ++i) {
    smoother->apply( *z, *y );
    z2->add( -1, *sol, 1, *y );
    n = z2->norm( Belos::InfNorm );
    if( space()->rankST()==0 )
      std::cout << "error: " << n << "\n";


  }

  TEST_EQUALITY_CONST( z2->norm( Belos::InfNorm )/z2->getLength()<eps, true );

	// Y test
  z->initField();

  y->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inY );
  y->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
  y->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inY );

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

    n = z2->norm( Belos::InfNorm );
    if( space()->rankST()==0 )
      std::cout << "error: " << n << "\n";

  }

  TEST_EQUALITY_CONST( z2->norm( Belos::InfNorm )/z2->getLength()<eps, true );

	// Z test
  z->initField();

  y->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inZ );
  y->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inZ );
  y->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );

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

    n = z2->norm( Belos::InfNorm );
    if( space()->rankST()==0 )
      std::cout << "error: " << n << "\n";

  }

  TEST_EQUALITY_CONST( z2->norm( Belos::InfNorm )/z2->getLength()<eps, true );

}




TEUCHOS_UNIT_TEST( BasicOperator, ConvectionDiffusionJSmoother ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<Pimpact::Space<ST,OT,sd,d,2> >( pl );

  auto op = Pimpact::create<ConvDiffOpT>( space );

	// init smoother
  auto pls = Teuchos::parameterList();
//  pls->set( "omega", 0.7 );
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
		wind->getFieldPtr(Pimpact::U)->initField( Pimpact::ConstField, 1. );
		wind->getFieldPtr(Pimpact::V)->initField( Pimpact::ConstField, 1. );
		wind->getFieldPtr(Pimpact::W)->initField( Pimpact::ConstField, 1. );

		op->assignField( *wind );
	}

	// init initial guess
  auto x = Pimpact::create<Pimpact::VectorField>( space );
	x->random();

	// init rhs
  auto y = Pimpact::create<Pimpact::VectorField>( space );


	ST error0 = x->norm( Belos::InfNorm );
	ST error;
	if( space()->rankST()==0 )
		std::cout << "\nerror: " << error0 << "\n";

	//x->write( 0 );
  for(int i=0; i<10; ++i) {
    smoother->apply( *y, *x );

    error = x->norm( Belos::InfNorm )/error0;
    if( space()->rankST()==0 )
      std::cout << "error: " << error << "\n";

		//x->write( i+1 );
  }

  TEST_EQUALITY_CONST( error<0.1, true );

}



//TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicConvectionOp ) {

	//const int d = 4;

	//setParameter( SpaceT::sdim );

	//pl->set<bool>( "spectral in time", true );

	//Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

	//auto vel = Pimpact::create<Pimpact::VectorField>( space );

	//auto mv1 = Pimpact::createMultiHarmonicVectorField( space );
	//auto mv2 = Pimpact::createMultiHarmonicVectorField( space );

	//auto op = Pimpact::createMultiHarmonicConvectionOp( space );

	//op->assignField( *mv1 );
	//op->apply( *mv1, *mv2 );

//}



//TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicDtConvectionDiffusionOp ) {

	//const int d = 4;

	//setParameter( SpaceT::sdim );
	//pl->set<bool>( "spectral in time", true );

	//OT nf = 4;
  //pl->set<ST>( "Re", 1.e0 );
  //pl->set<ST>( "alpha2", 12.0 );
  //pl->set<OT>( "nf", nf );

  //Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

	//ST re = space->getDomainSize()->getRe();
	//ST alpha2 = space->getDomainSize()->getAlpha2();

  //auto vel = Pimpact::create<Pimpact::VectorField>( space );

  //auto wind = Pimpact::createMultiHarmonicVectorField( space );
  //auto x    = Pimpact::createMultiHarmonicVectorField( space );
  //auto y1   = Pimpact::createMultiHarmonicVectorField( space );
  //auto y2   = Pimpact::createMultiHarmonicVectorField( space );
  //auto diff = Pimpact::createMultiHarmonicVectorField( space );


	//auto op1 =
		//Pimpact::createAdd2Op(
				//Pimpact::createMultiDtHelmholtz( space, alpha2/re, 1./re ),
				//Pimpact::createMultiHarmonicConvectionOp( space ) );

  //auto op2 =
		//Pimpact::createMultiDtConvectionDiffusionOp( space );

	//// init const wind 
	//wind->get0FieldPtr()->initField( Pimpact::ConstFlow, 1., 1., 1. );
	//for( int i=0; i<nf; ++i ) {
		//wind->getCFieldPtr(i)->initField( Pimpact::ConstFlow, 1., 1., 1. );
		//wind->getSFieldPtr(i)->initField( Pimpact::ConstFlow, 1., 1., 1. );
	//}

	//// init const dt
	//x->get0FieldPtr()->initField( Pimpact::ConstFlow, 1., 1., 1. );
	//for( int i=0; i<nf; ++i ) {
		//x->getCFieldPtr(i)->initField( Pimpact::ConstFlow, 1., 1., 1. );
		//x->getSFieldPtr(i)->initField( Pimpact::ConstFlow, 1., 1., 1. );
	//}

	//op1->assignField( *wind );
	//op2->assignField( *wind );

  //op1->apply( *x, *y1 );
  //op2->apply( *x, *y2 );

	//diff->add( 1., *y1, -1.,  *y2 );

	//std::cout << "\nconst error(dt): " << diff->norm( Belos::InfNorm ) << "\n";
  //TEST_EQUALITY_CONST( diff->norm( Belos::InfNorm ) <eps, true );

	//// init const udx 
	//x->get0FieldPtr()->initField( Pimpact::ConstFlow,0 );
////	x->get0FieldPtr()->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
////	x->get0FieldPtr()->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
////	x->get0FieldPtr()->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );
	//for( int i=0; i<nf; ++i ) {
////		x->getCFieldPtr(i)->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX );
////		x->getCFieldPtr(i)->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY );
////		x->getCFieldPtr(i)->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ );
		//x->getCFieldPtr(i)->initField( Pimpact::ConstFlow, 0. );
		//x->getSFieldPtr(i)->initField( Pimpact::ConstFlow, 0. );
		//x->getSFieldPtr(i)->getFieldPtr( Pimpact::U )->initField( Pimpact::Grad2D_inX, 1. );
////		x->getSFieldPtr(i)->getFieldPtr( Pimpact::V )->initField( Pimpact::Grad2D_inY, -1. );
////		x->getSFieldPtr(i)->getFieldPtr( Pimpact::W )->initField( Pimpact::Grad2D_inZ, -1. );
	//}

  //op1->apply( *x, *y1 );
  //op2->apply( *x, *y2 );

	//diff->add( 1., *y1, -1.,  *y2 );

	//ST diffs = diff->norm( Belos::InfNorm );
	//std::cout << "const error(du): " << diffs << "\n";
  //TEST_EQUALITY_CONST( diffs<eps, true );
	//if( diffs>=eps ) {
		////y1->write(0);
		////y2->write(100);
		////diff->write(1000);
	//}


	//// init const ddx 
	//x->get0FieldPtr()->initField( Pimpact::PoiseuilleFlow2D_inX );
	//for( int i=0; i<nf; ++i ) {
		//x->getCFieldPtr(i)->initField( Pimpact::PoiseuilleFlow2D_inX );
		//x->getSFieldPtr(i)->initField( Pimpact::PoiseuilleFlow2D_inX );
	//}

  //op1->apply( *x, *y1 );
  //op2->apply( *x, *y2 );

	//diff->add( 1., *y1, -1.,  *y2 );

	//std::cout << "error(ddx): " << diff->norm( Belos::InfNorm ) << "\n";
  //TEST_EQUALITY_CONST( diff->norm( Belos::InfNorm ) <eps, true );

//}





} // namespace
