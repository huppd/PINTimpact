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

const ST pi2 = std::atan(1.)*8.;

using SpaceT = Pimpact::Space<ST,OT,sd,d,dNC>;

using SF = typename Pimpact::ScalarField<SpaceT>;
using VF = typename Pimpact::VectorField<SpaceT>;
using MSF = typename Pimpact::ModeField<SF>;
using MVF = typename Pimpact::ModeField<VF>;

template<class T>
using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;
template<class T>
using ConvDiffJT = Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionJSmoother >;
template<class T>
using ConvDiffSORT = Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionSORSmoother >;

template<class T> using ConvOpT = Pimpact::NonlinearOp<Pimpact::ConvectionSOp<T> >;

using ConvDiffOp2D = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<D2> >;
using ConvDiffOp3D = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<D3> >;

using ConvDiffJ2D = Pimpact::NonlinearSmoother<ConvDiffOp2D,Pimpact::ConvectionDiffusionJSmoother >;
using ConvDiffSOR2D = Pimpact::NonlinearSmoother<ConvDiffOp2D,Pimpact::ConvectionDiffusionSORSmoother >;

using ConvDiffJ3D = Pimpact::NonlinearSmoother<ConvDiffOp3D,Pimpact::ConvectionDiffusionJSmoother >;
using ConvDiffSOR3D = Pimpact::NonlinearSmoother<ConvDiffOp3D,Pimpact::ConvectionDiffusionSORSmoother >;



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, ConvectionSOp, SpaceT ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

	Pimpact::ScalarField<SpaceT> u[3][3] = { 
		{ { space, true, Pimpact::F::U },
			{ space, true, Pimpact::F::U },
			{ space, true, Pimpact::F::U }},{
			{ space, true, Pimpact::F::V },
			{ space, true, Pimpact::F::V },
			{ space, true, Pimpact::F::V }},{
			{ space, true, Pimpact::F::W },
			{ space, true, Pimpact::F::W },
			{ space, true, Pimpact::F::W }}
	};

	Pimpact::VectorField<SpaceT> x   ( space );
	Pimpact::VectorField<SpaceT> y   ( space );
	Pimpact::VectorField<SpaceT> solv( space );

	auto op = Pimpact::create<Pimpact::ConvectionSOp>( space ) ;

	if( print ) op->print();

	std::cout << "\n";
	for( int dir=-1; dir<2; dir+= 2 ) {

		for( int i=0; i<SpaceT::sdim; ++i )
			for( int j=0; j<SpaceT::sdim; ++j )
				u[i][j].initField( Pimpact::ConstField, 2.*dir );

		// dx test
		solv(Pimpact::F::U).initField( Pimpact::ConstField, 2.*dir );
		solv(Pimpact::F::V).initField( Pimpact::ConstField, 2.*dir );
		if( 3==SpaceT::sdim )
			solv(Pimpact::F::W).initField( Pimpact::ConstField, 2.*dir );

		x(Pimpact::F::U).initField( Pimpact::Grad2D_inX );
		x(Pimpact::F::V).initField( Pimpact::Grad2D_inX );
		if( 3==SpaceT::sdim )
			x(Pimpact::F::W).initField( Pimpact::Grad2D_inX );

		y.random();

		for( Pimpact::F i=Pimpact::F::U; i<SpaceT::sdim; ++i )
			op->apply( u[static_cast<int>(i)], x(i), y(i) );

		for( Pimpact::F i=Pimpact::F::U; i<SpaceT::sdim; ++i ) {
			Pimpact::ScalarField<SpaceT>& sol = solv( i );
			sol.add( 1., sol, -1., y(i) );
			ST errorInf = sol.norm( Belos::InfNorm, Pimpact::B::N );
			std::cout << "error in "<< Pimpact::toString( static_cast<Pimpact::F>(i) ) << " (gradX): " << errorInf << "\n";
			TEST_EQUALITY( errorInf<eps, true );
			if( write && errorInf>=eps )
				sol.write();
		}

		// dy test
		solv(Pimpact::F::U).initField( Pimpact::ConstField, 2.*dir );
		solv(Pimpact::F::V).initField( Pimpact::ConstField, 2.*dir );
		if( 3==SpaceT::sdim )
			solv(Pimpact::F::W).initField( Pimpact::ConstField, 2.*dir );

		x(Pimpact::F::U).initField( Pimpact::Grad2D_inY );
		x(Pimpact::F::V).initField( Pimpact::Grad2D_inY );
		if( 3==SpaceT::sdim )
			x(Pimpact::F::W).initField( Pimpact::Grad2D_inY );

		y.random();

		for( Pimpact::F i=Pimpact::F::U; i<SpaceT::sdim; ++i )
			op->apply( u[static_cast<int>(i)], x(i), y(i) );

		for( Pimpact::F i=Pimpact::F::U; i<SpaceT::sdim; ++i ){
			Pimpact::ScalarField<SpaceT>& sol = solv( i );
			sol.add( 1., sol, -1., y(i) );
			std::cout << "error in "<< Pimpact::toString( static_cast<Pimpact::F>(i) ) << " (gradY): " << sol.norm( Belos::InfNorm, Pimpact::B::N ) << "\n";
			TEST_EQUALITY( sol.norm( Belos::InfNorm, Pimpact::B::N )<eps, true );
		}
		x(Pimpact::F::U).initField( Pimpact::Grad2D_inY );
		x(Pimpact::F::V).initField( Pimpact::Grad2D_inY );
		if( 3==SpaceT::sdim )
			x(Pimpact::F::W).initField( Pimpact::Grad2D_inY );

		// dz test
		if( 3==SpaceT::sdim ) {
			solv(Pimpact::F::U).initField( Pimpact::ConstField, 2.*dir );
			solv(Pimpact::F::V).initField( Pimpact::ConstField, 2.*dir );
			solv(Pimpact::F::W).initField( Pimpact::ConstField, 2.*dir );

			x(Pimpact::F::U).initField( Pimpact::Grad2D_inZ );
			x(Pimpact::F::V).initField( Pimpact::Grad2D_inZ );
			x(Pimpact::F::W).initField( Pimpact::Grad2D_inZ );

			y.random();

			for( Pimpact::F i=Pimpact::F::U; i<SpaceT::sdim; ++i )
				op->apply( u[static_cast<int>(i)], x(i), y(i) );

			for( Pimpact::F i=Pimpact::F::U; i<SpaceT::sdim; ++i ) {
				Pimpact::ScalarField<SpaceT>& sol = solv( i );
				sol.add( 1., sol, -1., y(i) );
				ST errorInf = sol.norm( Belos::InfNorm, Pimpact::B::N );
				std::cout << "error in "<< Pimpact::toString( static_cast<Pimpact::F>(i) ) << " (gradZ): " << errorInf << "\n";
				TEST_EQUALITY( errorInf<eps, true );
				if( write && errorInf>=eps )
					//sol.write();
					sol.print();
				//y(i).print();
			}
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionSOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionSOp, D3 )





TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, ConvectionVOp, SpaceT ) {

	setParameter( SpaceT::sdim );

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

	Pimpact::VectorField<SpaceT> x( space );
	Pimpact::VectorField<SpaceT> y( space );
	Pimpact::VectorField<SpaceT> z( space );
	Pimpact::VectorField<SpaceT> z2(space );

  auto op = Pimpact::create<ConvOpT>( space  ) ;

	for(int dir=-1; dir<2; dir+=2 ) {
		// x test
		for( Pimpact::F f = Pimpact::F::U; f<SpaceT::sdim; ++ f ) {
			x( f ).init( 2.*dir );
			y( f ).initField( Pimpact::Grad2D_inX );
			z2( f ).init( 2.*dir, Pimpact::B::N );
		}
		z.random();

		op->assignField( x );
		op->apply( y, z );

		z2.add( -1, z2, 1, z );

		if( print ) z2.print();
		ST error = z2.norm(Belos::InfNorm,Pimpact::B::N);
		if( 0==space->rankST() )
			std::cout << "\nerror in "<< dir << " (gradX): " << error << "\n";
		TEST_EQUALITY( error<eps, true );

		// y test
		x( Pimpact::F::U ).initField( Pimpact::ConstField, 2.*dir );
		x( Pimpact::F::V ).initField( Pimpact::ConstField, 2.*dir );
		x( Pimpact::F::W ).initField( Pimpact::ConstField, 2.*dir );
		y( Pimpact::F::U ).initField( Pimpact::Grad2D_inY );
		y( Pimpact::F::V ).initField( Pimpact::Grad2D_inY );
		y( Pimpact::F::W ).initField( Pimpact::Grad2D_inY );
		z.random();
		z2( Pimpact::F::U ).initField( Pimpact::ConstField, 2.*dir );
		z2( Pimpact::F::V ).initField( Pimpact::ConstField, 2.*dir );
		z2( Pimpact::F::W ).initField( Pimpact::ConstField, 2.*dir );

		op->assignField( x );
		op->apply( y, z );

		z2.add( -1, z2, 1, z );

		error = z2.norm(Belos::InfNorm,Pimpact::B::N);
		if( 0==space->rankST() )
			std::cout << "error in "<< dir << " (gradY): " << error << "\n";
		TEST_EQUALITY( error<eps, true );

		// z test
		if( 3==SpaceT::sdim ) {
			x( Pimpact::F::U ).initField( Pimpact::ConstField, 2.*dir );
			x( Pimpact::F::V ).initField( Pimpact::ConstField, 2.*dir );
			x( Pimpact::F::W ).initField( Pimpact::ConstField, 2.*dir );
			y( Pimpact::F::U ).initField( Pimpact::Grad2D_inZ );
			y( Pimpact::F::V ).initField( Pimpact::Grad2D_inZ );
			y( Pimpact::F::W ).initField( Pimpact::Grad2D_inZ );
			z.random();
			z2( Pimpact::F::U ).initField( Pimpact::ConstField, 2.*dir );
			z2( Pimpact::F::V ).initField( Pimpact::ConstField, 2.*dir );
			z2( Pimpact::F::W ).initField( Pimpact::ConstField, 2.*dir );

			op->assignField( x );
			op->apply( y, z );

			z2.add( -1, z2, 1, z );

			error = z2.norm(Belos::InfNorm,Pimpact::B::N);
			if( 0==space->rankST() )
				std::cout << "error in "<< dir << " (gradZ): " << error << "\n";
			TEST_EQUALITY( error<eps, true );
		}
	}

	// x test
	x.init();
	x( Pimpact::F::U ).initField( Pimpact::Poiseuille2D_inX );
	y( Pimpact::F::U ).initField( Pimpact::Grad2D_inX );
	y( Pimpact::F::V ).initField( Pimpact::Grad2D_inY );
	y( Pimpact::F::W ).initField( Pimpact::Grad2D_inZ );
	z.random();
	z2.init();
	z2( Pimpact::F::U ).initField( Pimpact::Poiseuille2D_inX );

	op->assignField( x );
	op->apply( y, z );

	z2.add( -1, z2, 1, z );

	ST error = z2.norm(Belos::InfNorm,Pimpact::B::N);
	if( 0==space->rankST() )
		std::cout << "error in (poisX): " << error << "\n";
	TEST_EQUALITY( error<eps, true );

	// y test
	x.init();
	x( Pimpact::F::V ).initField( Pimpact::Poiseuille2D_inY );
	y( Pimpact::F::U ).initField( Pimpact::Grad2D_inX );
	y( Pimpact::F::V ).initField( Pimpact::Grad2D_inY );
	y( Pimpact::F::W ).initField( Pimpact::Grad2D_inZ );
	z.random();
	z2.init();
	z2( Pimpact::F::V ).initField( Pimpact::Poiseuille2D_inY );

	op->assignField( x );
	op->apply( y, z );

	z2.add( -1, z2, 1, z );

	error = z2.norm(Belos::InfNorm,Pimpact::B::N);
	if( 0==space->rankST() )
		std::cout << "error in (poisY): " << error << "\n";
	TEST_EQUALITY( error<eps, true );

	// z test
	if( 3==SpaceT::sdim ) {
		x.init();
		x( Pimpact::F::W ).initField( Pimpact::Poiseuille2D_inZ );
		y( Pimpact::F::U ).initField( Pimpact::Grad2D_inX );
		y( Pimpact::F::V ).initField( Pimpact::Grad2D_inY );
		y( Pimpact::F::W ).initField( Pimpact::Grad2D_inZ );
		z.random();
		z2.init();
		z2( Pimpact::F::W ).initField( Pimpact::Poiseuille2D_inZ );

		op->assignField( x );
		op->apply( y, z );

		z2.add( -1, z2, 1, z );

		error = z2.norm(Belos::InfNorm,Pimpact::B::N);
		if( 0==space->rankST() )
			std::cout << "error in (poisZ): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionVOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionVOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, ConvectionDiffusionOp, SpaceT ) {

	setParameter( SpaceT::sdim );

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

	Pimpact::VectorField<SpaceT> x( space );
	Pimpact::VectorField<SpaceT> y( space );
	Pimpact::VectorField<SpaceT> z( space );
	Pimpact::VectorField<SpaceT> z2( space );
	Pimpact::VectorField<SpaceT> sol( space );


  auto op = Pimpact::createAdd2Op(
      Pimpact::create< ConvOpT<SpaceT> >(space),
      Pimpact::create<Pimpact::HelmholtzOp>(space) );

	auto op2 = Pimpact::create< ConvDiffOpT<SpaceT> >( space );


  x( Pimpact::F::U ).initField( Pimpact::ConstField, 2. );
  x( Pimpact::F::V ).initField( Pimpact::ConstField, 2. );
  x( Pimpact::F::W ).initField( Pimpact::ConstField, 2. );
  sol.init( Teuchos::tuple(2.,2.,2.) );

  // consistency test in x
  y(Pimpact::F::U).initField( Pimpact::Grad2D_inX );
  y(Pimpact::F::V).initField( Pimpact::Grad2D_inX );
  y(Pimpact::F::W).initField( Pimpact::Grad2D_inX );

  op->assignField( x );
  op2->assignField( x );

  // 
  op->apply(  y, z  );
  op2->apply( y, z2 );

  z2.add( -1, z2, 1, z );
  sol.add(-1, sol,1, z );

	std::cout << "diff (gradX): " << z2.norm( Belos::InfNorm, Pimpact::B::N  ) << "\n";
	std::cout << "erro (gradX): " << sol.norm( Belos::InfNorm, Pimpact::B::N  ) << "\n";
  TEST_EQUALITY( z2.norm( Belos::InfNorm, Pimpact::B::N  ) < eps, true );
  TEST_EQUALITY( sol.norm( Belos::InfNorm, Pimpact::B::N  ) < eps, true );

  // consistency test in y
  y(Pimpact::F::U).initField( Pimpact::Grad2D_inY );
  y(Pimpact::F::V).initField( Pimpact::Grad2D_inY );
  y(Pimpact::F::W).initField( Pimpact::Grad2D_inY );
  sol.init( Teuchos::tuple(2.,2.,2.) );

  op->assignField( x );
  op2->assignField( x );

  op->apply(   y, z );
  op2->apply(  y, z2 );

  z2.add( -1, z2, 1, z );
  sol.add(-1, sol,1, z );

	std::cout << "diff (gradY): " << z2.norm( Belos::InfNorm, Pimpact::B::N  ) << "\n";
	std::cout << "erro (gradY): " << sol.norm( Belos::InfNorm, Pimpact::B::N  ) << "\n";
  TEST_EQUALITY( z2.norm( Belos::InfNorm, Pimpact::B::N  ) < eps, true );
  TEST_EQUALITY( sol.norm( Belos::InfNorm, Pimpact::B::N  ) < eps, true );

  // consistency test in z
	if( 3==SpaceT::sdim ) {
		y(Pimpact::F::U).initField( Pimpact::Grad2D_inZ );
		y(Pimpact::F::V).initField( Pimpact::Grad2D_inZ );
		y(Pimpact::F::W).initField( Pimpact::Grad2D_inZ );
		sol.init( Teuchos::tuple(2.,2.,2.) );

		op->assignField( x );
		op2->assignField( x );

		op->apply(  y, z );
		op2->apply( y, z2 );

		z2.add( -1, z2, 1, z );
		sol.add(-1, sol,1, z );

		std::cout << "diff (gradZ): " << z2.norm( Belos::InfNorm, Pimpact::B::N  ) << "\n";
		std::cout << "erro (gradZ): " << sol.norm( Belos::InfNorm, Pimpact::B::N  ) << "\n";
		TEST_EQUALITY( z2.norm( Belos::InfNorm, Pimpact::B::N  ) < eps, true );
		TEST_EQUALITY( sol.norm( Belos::InfNorm, Pimpact::B::N  ) < eps, true );
	}


  // consistency test in pois x
	y.init();
  y(Pimpact::F::U).initField( Pimpact::Poiseuille2D_inX );
  sol.init( Teuchos::tuple(0.,0.,0.) );
  z.random();

  op->apply( y, z );
  op2->apply(  y, z2 );

  z2.add( -1, z2, 1, z );
  sol.add(-1, sol,1, z );

	std::cout << "diff (poisX): " << z2.norm( Belos::InfNorm, Pimpact::B::N  ) << "\n";
  TEST_EQUALITY( z2.norm( Belos::InfNorm, Pimpact::B::N  )<eps, true );


  // consistency test in pois y
	y.init();
  y(Pimpact::F::V).initField( Pimpact::Poiseuille2D_inY );
  sol.init( Teuchos::tuple(0.,0.,0.) );
  z.random();

  op->apply( y, z );
  op2->apply(  y, z2 );

  z2.add( -1, z2, 1, z );

	std::cout << "diff (poisY): " << z2.norm( Belos::InfNorm, Pimpact::B::N  ) << "\n";
  TEST_EQUALITY( z2.norm( Belos::InfNorm, Pimpact::B::N  )<eps, true );

  // consistency test in pois Z
	if( 3==SpaceT::sdim ) {
		y.init();
		y(Pimpact::F::W).initField( Pimpact::Poiseuille2D_inZ );
		z.random();

		op->apply( y, z );
		op2->apply(  y, z2 );

		z2.add( -1, z2, 1, z );

		std::cout << "diff (poisZ): " << z2.norm( Belos::InfNorm, Pimpact::B::N  ) << "\n";
		TEST_EQUALITY( z2.norm( Belos::InfNorm, Pimpact::B::N  )<eps, true );
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionDiffusionOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionDiffusionOp, D3 )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, ConvectionDiffusionSmoother, SType ) {

	setParameter( SType::SpaceT::sdim );

  auto space = Pimpact::create<typename SType::SpaceT>( pl );

  auto op = Pimpact::create<ConvDiffOpT>( space );

	// init smoother
  auto pls = Teuchos::parameterList();
	//pls->set( "omega", 0.5 );
	pls->set( "numIters", 1 );

	auto smoother = Pimpact::create<SType>( op, pls );

	//  init wind
	{

		Pimpact::VectorField<typename SType::SpaceT> wind( space );

		auto windfunc = [&pi2]( ST x ) -> ST {
			return( -std::cos(pi2*x/2) );
		};
		wind(Pimpact::F::U).initFromFunction(
				[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(x) ); } );
		wind(Pimpact::F::V).initFromFunction(
				[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(y) ); } );
		if( 3==SpaceT::sdim )
			wind(Pimpact::F::W).initFromFunction(
					[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(z) ); } );

		op->assignField( wind );
	}

	// init initial guess
	Pimpact::VectorField<typename SType::SpaceT> x( space );
	x.random();
	//x.init();

	// init rhs
	Pimpact::VectorField<typename SType::SpaceT> y( space );

	// Consistency
	Pimpact::VectorField<typename SType::SpaceT> xs( space );
	xs = x;
	op->apply( x, y );

	smoother->apply( y, x );
	//x.write(1);
	//xs.write(2);
	x.add( 1., xs, -1., x );
	//x.abs( x );
	ST error0 = x.norm( Belos::InfNorm );
	if( space()->rankST()==0 )
		std::cout << "\nConsistency error: " << error0 << "\n";
  TEST_EQUALITY_CONST( error0<eps, true );
	if( write ) x.write( 0 );
	if( print ) x.print();

	// Convergence
	y.init( 0. );
	x.random();
	error0 = x.norm( Belos::InfNorm );
	ST error;
	if( space()->rankST()==0 )
		std::cout << "\ninit error: " << error0 << "\n\n";

	if( space()->rankST()==0 )
		std::cout << "error:\n" << 1. << "\n";
	//x.write( 0 );
  for(int i=0; i<ns; ++i) {
    smoother->apply( y, x );

    error = x.norm( Belos::InfNorm )/error0;
    if( space()->rankST()==0 )
      std::cout << error << "\n";

		if( write ) x.write( i+1 );
  }

  TEST_EQUALITY_CONST( error<0.1, true );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionDiffusionSmoother, ConvDiffJ2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionDiffusionSmoother, ConvDiffSOR2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionDiffusionSmoother, ConvDiffJ3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ConvectionDiffusionSmoother, ConvDiffSOR3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ModeOperator, ModeNonlinearOp, SpaceT ) {

	pl->set<bool>( "spectral in time", true );

	ST pi2 = 2.*std::acos(-1.);

	lx = pi2 ;
	ly = pi2 ;
	lz = pi2 ;

	nf = 4;

	setParameter( SpaceT::sdim );
	Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );


	Pimpact::ModeField< Pimpact::VectorField<SpaceT> > x( space );
	Pimpact::ModeField< Pimpact::VectorField<SpaceT> > y( space );
	Pimpact::ModeField< Pimpact::VectorField<SpaceT> > sol( space );
	Pimpact::ModeField< Pimpact::VectorField<SpaceT> > err( space );

	auto zeroOp = Pimpact::create<ConvDiffOpT>( space );
	
	//Teuchos::RCP<const Pimpact::ModeNonlinearOp<ConvDiffOpT<SpaceT> > op =
	auto op =
		Teuchos::rcp( new Pimpact::ModeNonlinearOp< ConvDiffOpT<SpaceT> >(zeroOp));

	ST iRe = 1./space->getDomainSize()->getRe();
	ST a2 = space->getDomainSize()->getAlpha2()*iRe;

	// computing zero mode of z
	// set paramteters
	auto para = Teuchos::parameterList();
	para->set<ST>( "mulI", a2  );
	para->set<ST>( "mulC", 1.  );
	para->set<ST>( "mulL", iRe );
	op->setParameter( para );

	// initializtion
	{
		Pimpact::VectorField<SpaceT> wind( space );
		wind(Pimpact::F::U).init( 1. );
		wind(Pimpact::F::V).init( 1. );
		////wind(Pimpact::F::W).init( 1. );

		zeroOp->assignField( wind );
	}

	x.getSField()(Pimpact::F::U).initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::cos(x*pi2)*std::sin(y*pi2) ); } );
	x.getSField()(Pimpact::F::V).initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST { return( -std::sin(x*pi2)*std::cos(y*pi2) ); } );

	if( write ) x.write( 10 );

	// solution init
	
	sol.getCField()(Pimpact::F::U).initFromFunction(
	    	[&pi2,&alpha2,&re,&space]( ST x, ST y, ST z ) ->ST {
					if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
						(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
						(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
						(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
						(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
						(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
						return( 0. );
					else
						return( alpha2*std::cos(x*pi2)*std::sin(y*pi2)/re ); } );

	sol.getCField()(Pimpact::F::V).initFromFunction(
	    	[&pi2,&alpha2,&re,&space]( ST x, ST y, ST z ) ->ST {
					if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
						(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
						(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
						(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
						(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
						(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
						return( 0. );
					else
						return( -alpha2*std::sin(x*pi2)*std::cos(y*pi2)/re ); } );

	sol.getSField()(Pimpact::F::U).initFromFunction(
				[&pi2,&re,&space]( ST x, ST y, ST z ) ->ST {
					if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
						(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
						(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
						(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
						(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
						(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
						return(  std::cos( std::min(std::max(x,0.),1.)*pi2)*std::sin(std::min(std::max(y,0.),1.)*pi2) );
					else
						return(
							-    std::sin(x*pi2)*std::sin(y*pi2)
							+    std::cos(x*pi2)*std::cos(y*pi2)
							+ 2.*std::cos(x*pi2)*std::sin(y*pi2)/re ); } );
	sol.getSField()(Pimpact::F::V).initFromFunction(
				[&pi2,&re,&space]( ST x, ST y, ST z ) ->ST {
					if( (( x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
						(  ( x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
						(  ( y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
						(  ( y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
						(  ( z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
						(  ( z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
						return(  -std::sin(std::min(std::max(x,0.),1.)*pi2)*std::cos(std::min(std::max(y,0.),1.)*pi2) );
					else
						return(
							-    std::cos(x*pi2)*std::cos(y*pi2)
							+    std::sin(x*pi2)*std::sin(y*pi2)
							- 2.*std::sin(x*pi2)*std::cos(y*pi2)/re ); } );

	if( write ) sol.write( 30 );

	op->apply( x, y );

	if( write ) y.write( 20 );

	err.add( 1., sol, -1., y );
	if( write ) err.write( 0 );
	if( print ) err.print(   );

	ST error = err.norm(Belos::InfNorm)/sol.norm(Belos::InfNorm);
	std::cout << "\nerror: " << error << "\n";
	if( 1==domain )
		TEST_EQUALITY( error<(1./nx/ny), true );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ModeOperator, ModeNonlinearOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ModeOperator, ModeNonlinearOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiHarmonicOperator, MultiHarmonicDtConvectionDiffusionOp, SpaceT ) {

	pl->set<bool>( "spectral in time", true );

	ST pi2 = 2.*std::acos(-1.);

	lx = pi2 ;
	ly = pi2 ;
	lz = pi2 ;

	nf = 4;

	setParameter( SpaceT::sdim );
	Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );


	Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> > x( space );
	Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> > y( space );
	Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> > sol( space );
	Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> > err( space );

	auto op = Pimpact::createMultiDtConvectionDiffusionOp( space );

	// initializtion

	x.getSField(1)(Pimpact::F::U).initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::cos(x*pi2)*std::sin(y*pi2) ); } );
	x.getSField(1)(Pimpact::F::V).initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST { return( -std::sin(x*pi2)*std::cos(y*pi2) ); } );

	if( write ) x.write( 10 );

	// solution init
	
	sol.get0Field()(Pimpact::F::U).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST {
				if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
					(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
					(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
					(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
					(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
					(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
					return( 0. );
				else
					return(  -std::sin(2.*x*pi2)/4. ); } );
	sol.get0Field()(Pimpact::F::V).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST {
				if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
					(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
					(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
					(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
					(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
					(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
					return( 0. );
				else
					return(  -std::sin(2.*y*pi2)/4. ); } );

	sol.getCField(1)(Pimpact::F::U).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST {
				if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
					(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
					(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
					(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
					(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
					(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
					return( 0. );
				else
					return(  alpha2*std::cos(x*pi2)*std::sin(y*pi2)/re ); } );
	sol.getCField(1)(Pimpact::F::V).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST {
				if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
					(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
					(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
					(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
					(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
					(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
					return( 0. );
				else
					return( -alpha2*std::sin(x*pi2)*std::cos(y*pi2)/re ); } );

	sol.getSField(1)(Pimpact::F::U).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST {
				if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
					(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
					(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
					(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
					(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
					(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
				return(  std::cos(std::min(std::max(x,0.),1.)*pi2)*std::sin(y*pi2) );
				else
					return(  2.*std::cos(x*pi2)*std::sin(y*pi2)/re ); } );
	sol.getSField(1)(Pimpact::F::V).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST {
				if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
					(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
					(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
					(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
					(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
					(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
					return( -std::sin(x*pi2)*std::cos(std::max(std::min(y,1.),0.)*pi2) );
				else
					return( -2.*std::sin(x*pi2)*std::cos(y*pi2)/re ); } );

	sol.getCField(2)(Pimpact::F::U).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST {
				if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
					(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
					(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
					(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
					(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
					(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
					return( 0. );
				else
					return( std::sin(2.*x*pi2)/4. ); } );
	sol.getCField(2)(Pimpact::F::V).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST {
				if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
					(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
					(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
					(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
					(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
					(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
					return( 0. );
				else
					return( std::sin(2.*y*pi2)/4. ); } );

	if( write ) sol.write( 30 );

	op->assignField( x );
	op->apply( x, y );

	if( write ) y.write( 20 );

	err.add( 1., sol, -1., y );
	if( write ) err.write( 0 );

	ST error = err.norm()/sol.norm();
	std::cout << "\nerror: " << error << "\n";
	TEST_EQUALITY( error<(1./nx/ny), true );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicOperator, MultiHarmonicDtConvectionDiffusionOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiHarmonicOperator, MultiHarmonicDtConvectionDiffusionOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, ConvectionSOp, SpaceT ) { 

	int rank = 0;
	std::string label;
	ST pi2 = std::atan(1.)*8.;

	setParameter( SpaceT::sdim );

	for( Pimpact::F field=Pimpact::F::U; field<SpaceT::sdim; ++field ) {
		for( int dir=0; dir<SpaceT::sdim; ++dir ) {

			//  grid size
			pl->set<OT>( "nx", 9 );
			pl->set<OT>( "ny", 9 );
			pl->set<OT>( "nz", 9 );
			pl->set<OT>( "nf", 1 );

			std::vector<ST> error2( ns );
			std::vector<ST> errorInf( ns );
			std::vector<ST> dofs( ns );

			if( 0==rank )	
				std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

			for( OT n=0; n<ns; ++n ) {

				if( 0==dir ) 
					pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
				else if( 1==dir )
					pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
				else if( 2==dir )
					pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

				// grid stretching
				setStretching();

				auto space = Pimpact::create<SpaceT>( pl );
				rank = space->rankST();

				Pimpact::ScalarField<SpaceT> wind[3][3] = { 
					{ { space, true, Pimpact::F::U },
						{ space, true, Pimpact::F::U },
						{ space, true, Pimpact::F::U }},{
						{ space, true, Pimpact::F::V },
						{ space, true, Pimpact::F::V },
						{ space, true, Pimpact::F::V }},{
						{ space, true, Pimpact::F::W },
						{ space, true, Pimpact::F::W },
						{ space, true, Pimpact::F::W }}
				};
				auto windfunc = [&pi2]( ST x ) -> ST {
					return( -std::cos(pi2*x/2) );
				};
				wind[static_cast<int>(field)][0].initFromFunction(
						[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(x) ); } );
				wind[static_cast<int>(field)][1].initFromFunction(
						[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(y) ); } );
				wind[static_cast<int>(field)][2].initFromFunction(
						[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(z) ); } );

				Pimpact::VectorField<SpaceT> x( space );
				Pimpact::VectorField<SpaceT> y( space );
				Pimpact::VectorField<SpaceT> sol( space );

				// init 

				x.init();
				sol.init();
				auto inifunc = [&pi2]( ST x ) -> ST {  
					return(  std::cos(pi2*x) );
				};
				auto solfunc = [&pi2]( ST x ) -> ST {
					return( -pi2*std::sin(pi2*x) );
				};
 
				if( 0==dir ) {
					x(field).initFromFunction(
							[&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(x) ); } );
					sol(field).initFromFunction(
							[&field,&dir,&space,&solfunc,&windfunc]( ST x, ST y, ST z ) ->ST {
								if( field==dir && (
									( space->bcl(Pimpact::X)==Pimpact::BC::Dirichlet && x<Teuchos::ScalarTraits<ST>::eps()   ) ||
									( space->bcu(Pimpact::X)==Pimpact::BC::Dirichlet && x>1-Teuchos::ScalarTraits<ST>::eps() ) )
									) 
									return( 0. );
								else
									return( windfunc(x)*solfunc(x) ); } );
				}
				else if( 1==dir ) {
					x(field).initFromFunction( [&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(y) ); } );
					sol(field).initFromFunction(
							[&field,&dir,&space,&windfunc,&solfunc]( ST x, ST y, ST z ) ->ST {
								if( field==dir && (
									( space->bcl(Pimpact::Y)==Pimpact::BC::Dirichlet && y<Teuchos::ScalarTraits<ST>::eps()   ) ||
									( space->bcu(Pimpact::Y)==Pimpact::BC::Dirichlet && y>1-Teuchos::ScalarTraits<ST>::eps() ) )
									) 
									return( 0. );
								else
									return( windfunc(y)*solfunc(y) ); } );
				}
				else if( 2==dir ) {
					x(field).initFromFunction(
							[&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(z) ); } );
					sol(field).initFromFunction(
							[&field,&dir,&space,&windfunc,&solfunc]( ST x, ST y, ST z ) ->ST {
								if( field==dir && (
									( space->bcl(Pimpact::Z)==Pimpact::BC::Dirichlet && z<Teuchos::ScalarTraits<ST>::eps()   ) ||
									( space->bcu(Pimpact::Z)==Pimpact::BC::Dirichlet && z>1-Teuchos::ScalarTraits<ST>::eps() ) )
									) 
									return( 0. );
								else
									return( windfunc(z)*solfunc(z) ); } );
				}
				if( write ) x.write( 0 );
				if( write ) sol.write( 1 );
				if( print ) wind[static_cast<int>(field)][0].print();
				if( print ) wind[static_cast<int>(field)][1].print();
				if( print ) wind[static_cast<int>(field)][2].print();

				auto op = Pimpact::create<Pimpact::ConvectionSOp>( space ) ;
				label = op->getLabel();

				op->apply( wind[static_cast<int>(field)], x(field), y(field) );

				// compute error
				//if( print ) y.print();
				if( write ) y.write( n+10 );
				y(field).add( 1., sol(field), -1., y(field) );
				//if( print ) sol->print();
				if( write ) y.write( n+20 );
				if( write ) sol.write( n+30 );
				error2[n]   = std::log10( y(field).norm(Belos::TwoNorm) / sol(field).norm(Belos::TwoNorm) );
				errorInf[n] = std::log10( y(field).norm(Belos::InfNorm) / sol(field).norm(Belos::InfNorm) );
				dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
				if( 0==rank )	
					std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";
			}
			// compute order
			ST order2 = order<ST>( dofs, error2 );
			if( 0==rank )	
				std::cout << label << "(" << field << "): order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

			ST orderInf = order<ST>( dofs, errorInf );
			if( 0==rank )	
				std::cout << label << "(" << field << "): order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";
			// test
			TEST_EQUALITY( -order2  >2., true );
			TEST_EQUALITY( -orderInf>2., true );
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, ConvectionSOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, ConvectionSOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, ConvectionVOp, SpaceT ) { 

	int rank = 0;
	std::string label;
	ST pi2 = std::atan(1.)*8.;

	setParameter( SpaceT::sdim );

	for( Pimpact::F field=Pimpact::F::U; field<SpaceT::sdim; ++field ) {
		for( int dir=0; dir<SpaceT::sdim; ++dir ) {

			//  grid size
			pl->set<OT>( "nx", 9 );
			pl->set<OT>( "ny", 9 );
			pl->set<OT>( "nz", 9 );
			pl->set<OT>( "nf", 1 );

			std::vector<ST> error2( ns );
			std::vector<ST> errorInf( ns );
			std::vector<ST> dofs( ns );

			if( 0==rank )	
				std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

			for( OT n=0; n<ns; ++n ) {

				if( 0==dir ) 
					pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
				else if( 1==dir )
					pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
				else if( 2==dir )
					pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

				// grid stretching
				setStretching();

				auto space = Pimpact::create<SpaceT>( pl );
				rank = space->rankST();

				Pimpact::VectorField<SpaceT> wind(space);

				auto windfunc = [&pi2]( ST x ) -> ST {
					return( -std::cos(pi2*x/2) );
				};
				wind(Pimpact::F::U).initFromFunction(
						[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(x) ); } );
				wind(Pimpact::F::V).initFromFunction(
						[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(y) ); } );
				wind(Pimpact::F::W).initFromFunction(
						[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(z) ); } );

				Pimpact::VectorField<SpaceT> x( space );
				Pimpact::VectorField<SpaceT> y( space );
				Pimpact::VectorField<SpaceT> sol( space );

				// init 

				x.init();
				sol.init();
				auto inifunc = [&pi2]( ST x ) -> ST {  
					return(  std::cos(pi2*x) );
				};
				auto solfunc = [&pi2]( ST x ) -> ST {
					return( -pi2*std::sin(pi2*x) );
				};
 
				if( 0==dir ) {
					x(field).initFromFunction(
							[&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(x) ); } );
					sol(field).initFromFunction(
							[&field,&dir,&space,&solfunc,&windfunc]( ST x, ST y, ST z ) ->ST {
								if( field==dir && (
									( space->bcl(Pimpact::X)==Pimpact::BC::Dirichlet && x<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									( space->bcu(Pimpact::X)==Pimpact::BC::Dirichlet && x>=1-Teuchos::ScalarTraits<ST>::eps() ) )
									) 
									return( 0. );
								else
									return( windfunc(x)*solfunc(x) ); } );
				}
				else if( 1==dir ) {
					x(field).initFromFunction( [&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(y) ); } );
					sol(field).initFromFunction(
							[&field,&dir,&space,&windfunc,&solfunc]( ST x, ST y, ST z ) ->ST {
								if( field==dir && (
									( space->bcl(Pimpact::Y)==Pimpact::BC::Dirichlet && y<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									( space->bcu(Pimpact::Y)==Pimpact::BC::Dirichlet && y>=1-Teuchos::ScalarTraits<ST>::eps() ) )
									) 
									return( 0. );
								else
									return( windfunc(y)*solfunc(y) ); } );
				}
				else if( 2==dir ) {
					x(field).initFromFunction(
							[&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(z) ); } );
					sol(field).initFromFunction(
							[&field,&dir,&space,&windfunc,&solfunc]( ST x, ST y, ST z ) ->ST {
								if( field==dir && (
									( space->bcl(Pimpact::Z)==Pimpact::BC::Dirichlet && z<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									( space->bcu(Pimpact::Z)==Pimpact::BC::Dirichlet && z>=1-Teuchos::ScalarTraits<ST>::eps() ) )
									) 
									return( 0. );
								else
									return( windfunc(z)*solfunc(z) ); } );
				}
				if( write ) x.write( 0 );
				if( write ) sol.write( 1 );
				if( write ) wind.write( 2 );

				auto op = Pimpact::create<ConvOpT>( space  ) ;
				label = op->getLabel();

				op->assignField( wind );
				op->apply( x, y );

				// compute error
				//if( print ) y.print();
				if( write ) y.write( n+10 );
				y(field).add( 1., sol(field), -1., y(field) );
				//if( print ) sol->print();
				if( write ) y.write( n+20 );
				if( write ) sol.write( n+30 );
				error2[n]   = std::log10( y(field).norm(Belos::TwoNorm) / sol(field).norm(Belos::TwoNorm) );
				errorInf[n] = std::log10( y(field).norm(Belos::InfNorm) / sol(field).norm(Belos::InfNorm) );
				dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
				if( 0==rank )	
					std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";
			}
			// compute order
			ST order2 = order<ST>( dofs, error2 );
			if( 0==rank )	
				std::cout << label << "(" << field << "): order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

			ST orderInf = order<ST>( dofs, errorInf );
			if( 0==rank )	
				std::cout << label << "(" << field << "): order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";
			// test
			TEST_EQUALITY( -order2  >2., true );
			TEST_EQUALITY( -orderInf>2., true );
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, ConvectionVOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, ConvectionVOp, D3 )




/// \todo fix BC for outside points in solution-function
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, ConvectionDiffusionOp, SpaceT ) { 

	int rank = 0;
	std::string label;
	ST pi2 = std::atan(1.)*8.;

	setParameter( SpaceT::sdim );

	for( Pimpact::F field=Pimpact::F::U; field<SpaceT::sdim; ++field ) {
		for( int dir=0; dir<SpaceT::sdim; ++dir ) {

			//  grid size
			pl->set<OT>( "nx", 9 );
			pl->set<OT>( "ny", 9 );
			pl->set<OT>( "nz", 9 );
			pl->set<OT>( "nf", 1 );

			std::vector<ST> error2( ns );
			std::vector<ST> errorInf( ns );
			std::vector<ST> dofs( ns );

			if( 0==rank )	
				std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

			for( OT n=0; n<ns; ++n ) {

				if( 0==dir ) 
					pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
				else if( 1==dir )
					pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
				else if( 2==dir )
					pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

				// grid stretching
				setStretching();

				auto space = Pimpact::create<SpaceT>( pl );
				rank = space->rankST();
				ST ire = 1./space->getDomainSize()->getRe();

				Pimpact::VectorField<SpaceT> wind(space);

				auto windfunc = [&pi2]( ST x ) -> ST {
					return( -std::cos(pi2*x/2) );
				};
				wind(Pimpact::F::U).initFromFunction(
						[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(x) ); } );
				wind(Pimpact::F::V).initFromFunction(
						[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(y) ); } );
				wind(Pimpact::F::W).initFromFunction(
						[&windfunc]( ST x, ST y, ST z) ->ST { return( windfunc(z) ); } );

				Pimpact::VectorField<SpaceT> x( space );
				Pimpact::VectorField<SpaceT> y( space );
				Pimpact::VectorField<SpaceT> sol( space );

				// init 

				x.init();
				sol.init();
				auto inifunc = [&pi2]( ST x ) -> ST {  
					return(  std::cos(pi2*x) );
				};
				auto solfunc = [&pi2]( ST x ) -> ST {
					return( -pi2*std::sin(pi2*x) );
				};
				auto solfunc2 = [&pi2,&ire,&space]( ST x ) -> ST {  return(
						std::cos(x*pi2)*pi2*pi2*ire );
				};
 
				if( 0==dir ) {
					x(field).initFromFunction(
							[&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(x) ); } );
					sol(field).initFromFunction(
							[&field,&dir,&space,&inifunc,&solfunc,&solfunc2,&windfunc]( ST x, ST y, ST z ) ->ST {
								if( (space->bcl(Pimpact::X)==Pimpact::BC::Dirichlet && x<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									(  space->bcu(Pimpact::X)==Pimpact::BC::Dirichlet && x>=1-Teuchos::ScalarTraits<ST>::eps() ) ||
									(  space->bcl(Pimpact::Z)==Pimpact::BC::Dirichlet && z<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									(  space->bcu(Pimpact::Z)==Pimpact::BC::Dirichlet && z>=1-Teuchos::ScalarTraits<ST>::eps() ) ||
								  (  space->bcl(Pimpact::Y)==Pimpact::BC::Dirichlet && y<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									(  space->bcu(Pimpact::Y)==Pimpact::BC::Dirichlet && y>=1-Teuchos::ScalarTraits<ST>::eps() ) ) 
									return( inifunc(x) );
								else
									return( windfunc(x)*solfunc(x)+solfunc2(x) ); } );
				}
				else if( 1==dir ) {
					x(field).initFromFunction( [&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(y) ); } );
					sol(field).initFromFunction(
							[&field,&dir,&space,&inifunc,&solfunc,&solfunc2,&windfunc]( ST x, ST y, ST z ) ->ST {
								if( (space->bcl(Pimpact::X)==Pimpact::BC::Dirichlet && x<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									(  space->bcu(Pimpact::X)==Pimpact::BC::Dirichlet && x>=1-Teuchos::ScalarTraits<ST>::eps() ) ||
									(  space->bcl(Pimpact::Z)==Pimpact::BC::Dirichlet && z<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									(  space->bcu(Pimpact::Z)==Pimpact::BC::Dirichlet && z>=1-Teuchos::ScalarTraits<ST>::eps() ) ||
								  (  space->bcl(Pimpact::Y)==Pimpact::BC::Dirichlet && y<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									(  space->bcu(Pimpact::Y)==Pimpact::BC::Dirichlet && y>=1-Teuchos::ScalarTraits<ST>::eps() ) ) 
									return( inifunc(y) );
								else
									return( windfunc(y)*solfunc(y)+solfunc2(y) ); } );
				}
				else if( 2==dir ) {
					x(field).initFromFunction(
							[&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(z) ); } );
					sol(field).initFromFunction(
							[&field,&dir,&space,&inifunc,&solfunc,&solfunc2,&windfunc]( ST x, ST y, ST z ) ->ST {
								if( (space->bcl(Pimpact::X)==Pimpact::BC::Dirichlet && x<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									(  space->bcu(Pimpact::X)==Pimpact::BC::Dirichlet && x>=1-Teuchos::ScalarTraits<ST>::eps() ) ||
									(  space->bcl(Pimpact::Z)==Pimpact::BC::Dirichlet && z<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									(  space->bcu(Pimpact::Z)==Pimpact::BC::Dirichlet && z>=1-Teuchos::ScalarTraits<ST>::eps() ) ||
								  (  space->bcl(Pimpact::Y)==Pimpact::BC::Dirichlet && y<=Teuchos::ScalarTraits<ST>::eps()   ) ||
									(  space->bcu(Pimpact::Y)==Pimpact::BC::Dirichlet && y>=1-Teuchos::ScalarTraits<ST>::eps() ) ) 
									return( inifunc(z) );
								else
									return( windfunc(z)*solfunc(z)+solfunc2(z) ); } );
				}
				if( write ) x.write( 0 );
				if( write ) sol.write( 1 );
				if( write ) wind.write( 2 );

				auto op = Pimpact::create< ConvDiffOpT<SpaceT> >( space );
				label = op->getLabel();

				op->assignField( wind );
				op->apply( x, y );

				// compute error
				//if( print ) y.print();
				if( write ) y.write( n+10 );
				y(field).add( 1., sol(field), -1., y(field) );
				if( print ) y(field).print();
				if( write ) y.write( n+20 );
				if( write ) sol.write( n+30 );
				error2[n]   = std::log10( y(field).norm(Belos::TwoNorm) / sol(field).norm(Belos::TwoNorm) );
				errorInf[n] = std::log10( y(field).norm(Belos::InfNorm) / sol(field).norm(Belos::InfNorm) );
				dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
				if( 0==rank )	
					std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";
			}
			// compute order
			ST order2 = order<ST>( dofs, error2 );
			if( 0==rank )	
				std::cout << label << "(" << field << "): order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

			ST orderInf = order<ST>( dofs, errorInf );
			if( 0==rank )	
				std::cout << label << "(" << field << "): order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";
			// test
			TEST_EQUALITY( -order2  >2., true );
			TEST_EQUALITY( -orderInf>2., true );
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, ConvectionDiffusionOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, ConvectionDiffusionOp, D3 )



} // namespace
