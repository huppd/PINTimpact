#include <functional>
#include <iostream>
#include <cmath>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Fields.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_TransferOp.hpp"
#include "Pimpact_VectorFieldOpWrap.hpp"

#include "Pimpact_DivGradNullSpace.hpp"

#include "Pimpact_Test.hpp"





namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, DivOp, SpaceT ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

	Pimpact::ScalarField<SpaceT> p( space );
	Pimpact::VectorField<SpaceT> vel( space );
	Pimpact::ScalarField<SpaceT> sol( space );

	// zero test
	vel.init();

	auto op = Pimpact::create<Pimpact::DivOp>( space );

	if( print ) op->print();

	op->apply( vel, p );

	TEST_EQUALITY( p.norm( Belos::InfNorm )<eps, true );

	// random test
	vel.random();

	op->apply( vel, p );

	TEST_EQUALITY( p.norm( Belos::InfNorm )>eps, true );

	// circle XY test
	vel( Pimpact::F::U ).initField( Pimpact::Grad2D_inY, -1. );
	vel( Pimpact::F::V ).init();
	vel( Pimpact::F::W ).initField( Pimpact::Grad2D_inX, 1. );

	op->apply( vel, p );

	TEST_EQUALITY( p.norm( Belos::InfNorm )<eps, true );

	// circle XZ test
	vel( Pimpact::F::U ).initField( Pimpact::Grad2D_inY, -1. );
	vel( Pimpact::F::V ).initField( Pimpact::Grad2D_inX, 1. );
	vel( Pimpact::F::W ).init();

	op->apply( vel, p );

	TEST_EQUALITY( p.norm( Belos::InfNorm )<eps, true );

	// Grad test
	vel( Pimpact::F::U ).initField( Pimpact::Grad2D_inX );
	vel( Pimpact::F::V ).initField( Pimpact::Grad2D_inY, -2. );
	vel( Pimpact::F::W ).initField( Pimpact::Grad2D_inZ );

	p.random();

	op->apply( vel, p );

	sol.initField( Pimpact::ConstField, SpaceT::sdim );
	sol.add( 1., sol, -1., p );

	ST error = sol.norm( Belos::InfNorm );
	if( 0==rank )	std::cout << "GradXYZ error: " << error << "\n";
	TEST_EQUALITY( error>eps, true );

	// Grad test
	vel( Pimpact::F::U ).initField( Pimpact::Grad2D_inZ );
	vel( Pimpact::F::V ).initField( Pimpact::Grad2D_inX );
	vel( Pimpact::F::W ).initField( Pimpact::Grad2D_inY );

	p.random();

	op->apply( vel, p );

	error = p.norm( Belos::InfNorm );
	if( 0==rank )	std::cout << "GradZXY error: " << error << "\n";
	TEST_EQUALITY( error<eps, true );

	// Grad test
	vel( Pimpact::F::U ).initField( Pimpact::Grad2D_inY );
	vel( Pimpact::F::V ).initField( Pimpact::Grad2D_inZ );
	vel( Pimpact::F::W ).initField( Pimpact::Grad2D_inX );

	p.random();

	op->apply( vel, p );

	error = p.norm( Belos::InfNorm );
	if( 0==rank ) std::cout << "erro: " << error << "\n";
	TEST_EQUALITY( error<eps, true );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, InterpolateV2SOp, SpaceT ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

	if( print ) space->print();

	Pimpact::ScalarField<SpaceT> p( space );
	Pimpact::ScalarField<SpaceT> sol( space );
	Pimpact::VectorField<SpaceT> vel( space );

  auto op = Pimpact::createInterpolateV2S( space );
	
	if( print )
		op->print();

	if( 0==space->rankST() ) std::cout << "\n";
	for( Pimpact::F i=Pimpact::F::U; i<SpaceT::sdim; ++i ) {

		for( int bla=0; bla<=((domain!=1)?3:0); ++bla ) {

			Pimpact::EScalarField type =	static_cast<Pimpact::EScalarField>( bla );

			vel( i ).initField( type );
			sol.initField( type );

			p.random();
			op->apply( vel(i), p );
			sol.add( 1., sol, -1., p );

			ST error = sol.norm( Belos::InfNorm );
			if( 0==space->rankST() ) std::cout << "field " << i <<
				", error("<< type <<"): " << error << "\n";
			TEST_EQUALITY( error<eps, true );
			if( error>= eps ){
				std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
				if( print )
				sol.print(
						*Pimpact::createOstream( "error_int"+Pimpact::toString(i)+"2S_"+Pimpact::toString(type)+"_r"+r+".txt" ));
				if( write ) sol.write();
			}
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, InterpolateV2SOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, InterpolateV2SOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, InterpolateS2VOp, SpaceT ) {

	setParameter( SpaceT::sdim );

	auto space = Pimpact::create<SpaceT>( pl );

	Pimpact::ScalarField<SpaceT> p   ( space );
	Pimpact::VectorField<SpaceT> vel ( space );
	Pimpact::VectorField<SpaceT> solv( space );

	auto op = Pimpact::create<Pimpact::InterpolateS2V>( space );

	if( print )
		op->print();

	if( 0==rank ) std::cout << "\n";

	for( Pimpact::F i=Pimpact::F::U; i<SpaceT::sdim; ++i ) {

		Pimpact::ScalarField<SpaceT>& sol = solv( i );

		for( int bla=0; bla<=(domain!=1?6:0); ++bla ) {

			Pimpact::EScalarField type =	static_cast<Pimpact::EScalarField>(bla);

			vel( i ).random();
			sol.initField( type );
			p.initField( type );

			op->apply( p, vel( i ) );
			sol.add( 1., sol, -1., vel(i) );

			ST error = sol.norm( Belos::InfNorm );
			if( 0==rank ) std::cout << "field " << i
					<< ", error("<< type <<"): "
					<< error << "\n";

			TEST_EQUALITY( error<eps, true );
			if( error>=eps ) {
				std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
				sol.print(
						*Pimpact::createOstream( "error_int_S2"+Pimpact::toString( i )+"_"+Pimpact::toString(type)+"_r"+r+".txt" ));
			}

		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, InterpolateS2VOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, InterpolateS2VOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, TransferOp, SpaceT ) {

	setParameter( SpaceT::sdim );

  using FSpaceT = Pimpact::Space<ST,OT,SpaceT::sdim,d,4>;
  using CSpaceT = Pimpact::Space<ST,OT,SpaceT::sdim,d,2>;

  auto fSpace = Pimpact::create<FSpaceT>( pl );
  auto cSpace = Pimpact::create<CSpaceT>( pl );

	Pimpact::ScalarField<FSpaceT> fx( fSpace );
  Pimpact::ScalarField<CSpaceT> cx( cSpace );

	Pimpact::ScalarField<FSpaceT> fxs( fSpace );
  Pimpact::ScalarField<CSpaceT> cxs( cSpace );

  auto op = Pimpact::create< Pimpact::TransferOp<FSpaceT,CSpaceT> >( fSpace, cSpace );

  // test
  fx.initField( Pimpact::Poiseuille2D_inX );
  cx.random();
  cxs.initField( Pimpact::Poiseuille2D_inX );

  op->apply( fx, cx );

	cxs.add( 1., cxs, -1., cx );
  TEST_EQUALITY( cxs.norm(Belos::InfNorm)<eps, true );


  cx.initField( Pimpact::Poiseuille2D_inX );
  fx.random();
  fxs.initField( Pimpact::Poiseuille2D_inX );

  op->apply( cx, fx );

	fxs.add( 1., fxs, -1., fx );
	
  TEST_EQUALITY( fxs.norm(Belos::InfNorm)<eps, true );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, TransferOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, TransferOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, VectorFieldOpWrap, SpaceT ) {

	setParameter( SpaceT::sdim );

  using FSpaceT = Pimpact::Space<ST,OT,3,d,4>;
  using CSpaceT = Pimpact::Space<ST,OT,3,d,2>;

  auto fSpace = Pimpact::create<FSpaceT>( pl );
  auto cSpace = Pimpact::create<CSpaceT>( pl );

	Pimpact::VectorField<FSpaceT>fx( fSpace );
  Pimpact::VectorField<CSpaceT>cx( cSpace );

	auto op = Pimpact::create< Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<FSpaceT,CSpaceT> > >( fSpace, cSpace );

  // test
  fx.init();
  fx( Pimpact::F::U ).initField( Pimpact::Poiseuille2D_inX );
  cx.random();

  op->apply( fx, cx );

  TEST_FLOATING_EQUALITY( fx.norm(Belos::OneNorm), cx.norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx.norm(Belos::TwoNorm), cx.norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx.norm(Belos::InfNorm), cx.norm(Belos::InfNorm), eps );

  fx.random();

  op->apply( cx, fx );

  TEST_FLOATING_EQUALITY( fx.norm(Belos::OneNorm), cx.norm(Belos::OneNorm), eps );
  TEST_FLOATING_EQUALITY( fx.norm(Belos::TwoNorm), cx.norm(Belos::TwoNorm), eps );
  TEST_FLOATING_EQUALITY( fx.norm(Belos::InfNorm), cx.norm(Belos::InfNorm), eps );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, VectorFieldOpWrap, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, VectorFieldOpWrap, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, GradOp, SpaceT ) {

	if( domain!=1 ) {

		setParameter( SpaceT::sdim );

		auto space = Pimpact::create<SpaceT>( pl );

		Pimpact::ScalarField<SpaceT> p( space );
		Pimpact::VectorField<SpaceT> v( space );
		Pimpact::VectorField<SpaceT> sol( space );

		auto op = Pimpact::create<Pimpact::GradOp>( space );

		if( print )
			op->print();

		// op in x test
		p.initField( Pimpact::Grad2D_inX );
		sol.init();
		sol( Pimpact::F::U ).initField( Pimpact::ConstField, 1. );
		v.random();

		op->apply( p, v );

		sol.add( 1., sol, -1., v );

		ST error = sol.norm( Belos::InfNorm, Pimpact::B::N );
		if( 0==rank ) std::cout << "error: " << error << "\n";
		TEST_EQUALITY( error<eps, true );

		// op in y test
		p.initField( Pimpact::Grad2D_inY );
		sol.init();
		sol( Pimpact::F::V ).initField( Pimpact::ConstField, 1. );
		v.random();

		op->apply( p, v );

		sol.add( 1., sol, -1., v );

		error = sol.norm( Belos::InfNorm, Pimpact::B::N  );
		if( 0==rank ) std::cout << "error: " << error << "\n";
		TEST_EQUALITY( error<eps, true );

		// op in z test
		p.initField( Pimpact::Grad2D_inZ );
		sol.init();
		sol( Pimpact::F::W ).initField( Pimpact::ConstField, 1. );
		v.random();

		op->apply( p, v );

		sol.add( 1., sol, -1., v );

		error = sol.norm( Belos::InfNorm, Pimpact::B::N  );
		if( 0==rank ) std::cout << "error: " << error << "\n";
		TEST_EQUALITY( error<eps, true );
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, GradOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, GradOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, HelmholtzOp, SpaceT ) {

	if( domain!= 1 ) {

		setParameter( SpaceT::sdim );

		ST mulI = 5.;
		ST mulL = 3.;

		pl->set<ST>( "alpha2", mulI );
		pl->set<ST>( "Re", 1./mulL );

		auto space = Pimpact::create<SpaceT>( pl );

		Pimpact::VectorField<SpaceT> x(space);
		Pimpact::VectorField<SpaceT> b(space);
		Pimpact::VectorField<SpaceT> bs(space);

		auto op = Pimpact::create<Pimpact::HelmholtzOp>( space );

		if( print )
			op->print();

		// test in x direction
		x.init();
		x( Pimpact::F::U ).initField( Pimpact::Poiseuille2D_inX );
		bs( Pimpact::F::U ).init( 8./std::pow(space->getDomainSize()->getSize(Pimpact::Y),2) );
		bs( Pimpact::F::V ).init( 0. );
		bs( Pimpact::F::W ).init( 0. );

		bs.add( mulI, x, mulL, bs );

		auto para = Teuchos::parameterList();
		para->set<ST>( "mulI", mulI );
		op->setParameter( para );
		op->apply( x, b );

		bs.add( 1., bs, -1, b);

		ST error = bs.norm( Belos::InfNorm, Pimpact::B::N );
		if( 0==rank ) std::cout << "error(poiseuilleX): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			if( print ) bs.print( *Pimpact::createOstream(
						"error_helm_poiX_r"+r+".txt" ));
		}

		// test in y direction
		x.init();
		x( Pimpact::F::V ).initField( Pimpact::Poiseuille2D_inY );
		bs( Pimpact::F::U ).init( 0. );
		bs( Pimpact::F::V ).init( 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2) );
		bs( Pimpact::F::W ).init( 0. );
		bs.add( mulI, x, mulL, bs );

		op->apply( x, b );

		bs.add( 1., bs, -1, b );

		error = bs.norm( Belos::InfNorm, Pimpact::B::N );
		if( 0==rank ) std::cout << "error(poiseuilleY): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			if( print ) bs.print(
					*Pimpact::createOstream( "error_helm_poiY_r"+r+".txt" ));
		}

		// test in z direction
		x.init();
		x( Pimpact::F::W ).initField( Pimpact::Poiseuille2D_inZ );

		bs( Pimpact::F::U ).init( 0. );
		bs( Pimpact::F::V ).init( 0. );
		bs( Pimpact::F::W ).init( 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2) );

		bs.add( mulI, x, mulL, bs );

		op->apply( x, b );

		bs.add( 1., bs, -1, b);

		error = bs.norm( Belos::InfNorm, Pimpact::B::N );
		if( 0==rank ) std::cout << "error(poiseuilleZ): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			if( print ) bs.print(
					*Pimpact::createOstream( "error_helm_poiZ_r"+r+".txt" ));
		}

		// the circle XY test
		x.init();
		x( Pimpact::F::U ).initField( Pimpact::Grad2D_inY, -1. );
		x( Pimpact::F::W ).initField( Pimpact::Grad2D_inX,  1. );
		bs.init();
		bs( Pimpact::F::U ).initField( Pimpact::Grad2D_inY, -1. );
		bs( Pimpact::F::W ).initField( Pimpact::Grad2D_inX,  1. );
		bs.scale( mulI );

		op->apply( x, b );

		bs.add( 1., bs, -1, b );

		error = bs.norm( Belos::InfNorm, Pimpact::B::N );
		if( 0==rank ) std::cout << "error(circleXY): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			if( print ) bs.print(
					*Pimpact::createOstream( "error_helm_circXY_r"+r+".txt" ));
		}

		// the circle XZ test
		x.init();
		x( Pimpact::F::U ).initField( Pimpact::Grad2D_inY, -1. );
		x( Pimpact::F::V ).initField( Pimpact::Grad2D_inX,  1. );
		bs.init();
		bs( Pimpact::F::U ).initField( Pimpact::Grad2D_inY, -1. );
		bs( Pimpact::F::V ).initField( Pimpact::Grad2D_inX,  1. );

		bs.scale( mulI );

		op->apply( x, b );

		bs.add( 1., bs, -1, b );

		error = bs.norm( Belos::InfNorm, Pimpact::B::N );
		if( 0==rank ) std::cout << "error(circleXZ): " << error << "\n";
		TEST_EQUALITY( error<eps, true );
		if( error>= eps ){
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			if( print ) bs.print(
					*Pimpact::createOstream( "error_helm_circXZ_r"+r+".txt" ));
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, HelmholtzOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, HelmholtzOp, D3 )



//#ifdef NDEBUG
TEUCHOS_UNIT_TEST( MaatrixTest, DivGradOp2M ) {

	OT nx = 7;
	OT ny = 7;
	OT nz = 7;

	//OT nx = 9;
	//OT ny = 9;
	//OT nz = 9;
	
	//OT nx = 33;
	//OT ny = 33;
	//OT nz = 33;

	Pimpact::setBoundaryConditions( pl, domain );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	setStretching();

	// processor grid size
	pl->set( "npx", npx );
	pl->set( "npy", npy );
	pl->set( "npz", npz );
	pl->set( "npf", npf );

	auto space = Pimpact::create<D3>( pl );

	Pimpact::ScalarField<D3> x   ( space );
	Pimpact::ScalarField<D3> x2  ( space );
	Pimpact::ScalarField<D3> b   ( space );
	Pimpact::ScalarField<D3> ones( space );
	Pimpact::ScalarField<D3> diag( space );

	auto op = Pimpact::create<Pimpact::DivGradOp>( space );
	if( print ) op->print();

	ones.init( 1. );
	op->applyInvDiag( ones, diag );

	if( write ) diag.write(999);

	OT nxS = space->ei(Pimpact::F::S,Pimpact::X) - space->si(Pimpact::F::S,Pimpact::X) + 1;
	OT nyS = space->ei(Pimpact::F::S,Pimpact::Y) - space->si(Pimpact::F::S,Pimpact::Y) + 1;
	OT nzS = space->ei(Pimpact::F::S,Pimpact::Z) - space->si(Pimpact::F::S,Pimpact::Z) + 1;
	OT nS = nxS*nyS*nzS;


	Teuchos::SerialDenseMatrix<OT,ST> DJG( nS, nS );

	// test DJG
	if( npx*npy*npz==1 ) {
		Teuchos::RCP<std::ostream> output = Pimpact::createOstream( "DivGradOp.txt" );
		*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

		for( OT k=space->si(Pimpact::F::S,Pimpact::Z); k<=space->ei(Pimpact::F::S,Pimpact::Z); ++k )
			for( OT j=space->si(Pimpact::F::S,Pimpact::Y); j<=space->ei(Pimpact::F::S,Pimpact::Y); ++j )
				for( OT i=space->si(Pimpact::F::S,Pimpact::X); i<=space->ei(Pimpact::F::S,Pimpact::X); ++i ) {
					OT II = i-space->si(Pimpact::F::S,Pimpact::X)
						+ nxS*( j-space->si(Pimpact::F::S,Pimpact::Y) )
						+ nxS*nyS*( k-space->si(Pimpact::F::S,Pimpact::Z) );
					x.init( 0. );
					x(i,j,k) = 1.;
					op->apply( x, b );

					for( OT kk=space->si(Pimpact::F::S,Pimpact::Z); kk<=space->ei(Pimpact::F::S,Pimpact::Z); ++kk )
						for( OT jj=space->si(Pimpact::F::S,Pimpact::Y); jj<=space->ei(Pimpact::F::S,Pimpact::Y); ++jj )
							for( OT ii=space->si(Pimpact::F::S,Pimpact::X); ii<=space->ei(Pimpact::F::S,Pimpact::X); ++ii ) {
								OT JJ = ii-space->si(Pimpact::F::S,Pimpact::X)
									+ nxS*( jj-space->si(Pimpact::F::S,Pimpact::Y) )
									+ nxS*nyS*( kk-space->si(Pimpact::F::S,Pimpact::Z) );
								x2.init( 0. );
								x2(ii,jj,kk) = 1.;

								DJG(II,JJ) = x2.dot( b );
							}
					ST errorDiag = std::abs( 1./std::abs(DJG(II,II))-diag(i,j,k)) /
						std::abs(1./std::abs(DJG(II,II)) );
					if( errorDiag>=eps ) {
						std::cout << "diag("<<i<<", "<<j<< ", " <<k<<")\t" << errorDiag <<
							"\n";
					}
					TEST_EQUALITY( errorDiag<eps, true );
				}
		DJG.print( *output );
	}

	// test DJG^T
	Teuchos::SerialDenseMatrix<OT,ST> DJGT( nS, nS );
	if( npx*npy*npz==1 ) {
		Teuchos::RCP<std::ostream> output = Pimpact::createOstream( "DivGradOpT.txt" );
		*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

		for( OT k=space->si(Pimpact::F::S,Pimpact::Z); k<=space->ei(Pimpact::F::S,Pimpact::Z); ++k )
			for( OT j=space->si(Pimpact::F::S,Pimpact::Y); j<=space->ei(Pimpact::F::S,Pimpact::Y); ++j )
				for( OT i=space->si(Pimpact::F::S,Pimpact::X); i<=space->ei(Pimpact::F::S,Pimpact::X); ++i ) {
					OT II = i-space->si(Pimpact::F::S,Pimpact::X)
						+ nxS*( j-space->si(Pimpact::F::S,Pimpact::Y) )
						+ nxS*nyS*( k-space->si(Pimpact::F::S,Pimpact::Z) );
					x.init( 0. );
					x(i,j,k) = 1.;
					op->apply( x, b, Belos::TRANS );

					for( OT kk=space->si(Pimpact::F::S,Pimpact::Z); kk<=space->ei(Pimpact::F::S,Pimpact::Z); ++kk )
						for( OT jj=space->si(Pimpact::F::S,Pimpact::Y); jj<=space->ei(Pimpact::F::S,Pimpact::Y); ++jj )
							for( OT ii=space->si(Pimpact::F::S,Pimpact::X); ii<=space->ei(Pimpact::F::S,Pimpact::X); ++ii ) {
								OT JJ = ii-space->si(Pimpact::F::S,Pimpact::X)
									+ nxS*( jj-space->si(Pimpact::F::S,Pimpact::Y) )
									+ nxS*nyS*( kk-space->si(Pimpact::F::S,Pimpact::Z) );
								x2.init( 0. );
								x2(ii,jj,kk) = 1.;
								DJGT(JJ,II) = x2.dot( b );
							}
				}
		DJGT.print( *output );
	}
	DJG -= DJGT;

	ST errOne = DJG.normOne();
	std::cout << "\n\n||DJG-DJG^TT||_1 = " << errOne << "\n";
	TEST_EQUALITY( errOne<eps, true );

	ST errInf = DJG.normInf();
	std::cout << "||DJG-DJG^TT||_infty = " << errInf << "\n";
	TEST_EQUALITY( errInf<eps, true );
}
//#endif


#ifdef NDEBUG
TEUCHOS_UNIT_TEST( MatrixTest, DivOp2M ) {

	OT nx = 7;
	OT ny = 7;
	OT nz = 7;

	//OT nx = 9;
	//OT ny = 9;
	//OT nz = 9;
	
	//OT nx = 13;
	//OT ny = 13;
	//OT nz = 13;

	//OT nx = 17;
	//OT ny = 17;
	//OT nz = 17;

  Pimpact::setBoundaryConditions( pl, domain );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	setStretching();

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

  auto space = Pimpact::create<D3>( pl );

	Pimpact::VectorField<D3> x   ( space );
	Pimpact::ScalarField<D3> x2  ( space );
	Pimpact::ScalarField<D3> b   ( space );
	Pimpact::ScalarField<D3> ones( space );

  auto op = Pimpact::create<Pimpact::DivOp>( space );


	OT nxS = space->ei(Pimpact::F::S,Pimpact::X,Pimpact::B::Y) - space->si(Pimpact::F::S,Pimpact::X,Pimpact::B::Y) + 1;
	OT nyS = space->ei(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y) - space->si(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y) + 1;
	OT nzS = space->ei(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y) - space->si(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y) + 1;
	OT nS = nxS*nyS*nzS;

	OT nxU = space->ei(Pimpact::F::U,Pimpact::X,Pimpact::B::Y) - space->si(Pimpact::F::U,Pimpact::X,Pimpact::B::Y) + 1;
	OT nyU = space->ei(Pimpact::F::U,Pimpact::Y,Pimpact::B::Y) - space->si(Pimpact::F::U,Pimpact::Y,Pimpact::B::Y) + 1;
	OT nzU = space->ei(Pimpact::F::U,Pimpact::Z,Pimpact::B::Y) - space->si(Pimpact::F::U,Pimpact::Z,Pimpact::B::Y) + 1;
	OT nU = nxU*nyU*nzU;

	OT nxV = space->ei(Pimpact::F::V,Pimpact::X,Pimpact::B::Y) - space->si(Pimpact::F::V,Pimpact::X,Pimpact::B::Y) + 1;
	OT nyV = space->ei(Pimpact::F::V,Pimpact::Y,Pimpact::B::Y) - space->si(Pimpact::F::V,Pimpact::Y,Pimpact::B::Y) + 1;
	OT nzV = space->ei(Pimpact::F::V,Pimpact::Z,Pimpact::B::Y) - space->si(Pimpact::F::V,Pimpact::Z,Pimpact::B::Y) + 1;
	OT nV = nxV*nyV*nzV;

	OT nxW = space->ei(Pimpact::F::W,Pimpact::X,Pimpact::B::Y) - space->si(Pimpact::F::W,Pimpact::X,Pimpact::B::Y) + 1;
	OT nyW = space->ei(Pimpact::F::W,Pimpact::Y,Pimpact::B::Y) - space->si(Pimpact::F::W,Pimpact::Y,Pimpact::B::Y) + 1;
	OT nzW = space->ei(Pimpact::F::W,Pimpact::Z,Pimpact::B::Y) - space->si(Pimpact::F::W,Pimpact::Z,Pimpact::B::Y) + 1;
	OT nW = nxW*nyW*nzW;

	Teuchos::Tuple<OT,3> NX = Teuchos::tuple( nxU,nxV,nxW );
	Teuchos::Tuple<OT,3> NY = Teuchos::tuple( nyU,nyV,nyW );
	Teuchos::Tuple<OT,3> NZ = Teuchos::tuple( nzU,nzV,nzW );
	Teuchos::Tuple<OT,3> offset = Teuchos::tuple( 0, nU, nU+nV );

	Teuchos::SerialDenseMatrix<OT,ST> D  ( nS      , nU+nV+nW );

	// test D

	Teuchos::RCP<std::ostream> output = Pimpact::createOstream( "DivOp.txt" );
	*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

	for( int dir=0; dir<D3::sdim; ++dir ) {
		Pimpact::F fdir = static_cast<Pimpact::F>(dir);
		for( OT k=space->si(fdir,Pimpact::Z,Pimpact::B::Y); k<=space->ei(fdir,Pimpact::Z,Pimpact::B::Y); ++k )
			for( OT j=space->si(fdir,Pimpact::Y,Pimpact::B::Y); j<=space->ei(fdir,Pimpact::Y,Pimpact::B::Y); ++j )
				for( OT i=space->si(fdir,Pimpact::X,Pimpact::B::Y); i<=space->ei(fdir,Pimpact::X,Pimpact::B::Y); ++i ) {
					OT II =               i-space->si(fdir,Pimpact::X,Pimpact::B::Y)
						+ NX[dir]*(         j-space->si(fdir,Pimpact::Y,Pimpact::B::Y) )
						+ NX[dir]*NY[dir]*( k-space->si(fdir,Pimpact::Z,Pimpact::B::Y) ) + offset[dir];
					x.init();
					x(fdir)(i,j,k) = 1.;
					op->apply( x, b );

					for( OT kk=space->si(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y); kk<=space->ei(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y); ++kk )
						for( OT jj=space->si(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y); jj<=space->ei(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y); ++jj )
							for( OT ii=space->si(Pimpact::F::S,Pimpact::X,Pimpact::B::Y); ii<=space->ei(Pimpact::F::S,Pimpact::X,Pimpact::B::Y); ++ii ) {
								OT JJ =       ii-space->si(Pimpact::F::S,Pimpact::X,Pimpact::B::Y)
									+ nxS*(     jj-space->si(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y) )
									+ nxS*nyS*( kk-space->si(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y) );
								x2.init();
								x2(ii,jj,kk) = 1.;

								D(JJ,II) = x2.dot( b );
								//*output << D(JJ,II) << "\t";
							}
					//*output << "\n";
				}
	}
	D.print( *output );
}
#endif



#ifdef NDEBUG
TEUCHOS_UNIT_TEST( MatrixTest, GradOp2M ) {

	OT nx = 7;
	OT ny = 7;
	OT nz = 7;

	//OT nx = 9;
	//OT ny = 9;
	//OT nz = 9;
	
	//OT nx = 13;
	//OT ny = 13;
	//OT nz = 13;
	
	//OT nx = 17;
	//OT ny = 17;
	//OT nz = 17;

  Pimpact::setBoundaryConditions( pl, domain );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	setStretching();

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

  auto space = Pimpact::create<D3>( pl );

	Pimpact::VectorField<D3> x ( space );
	Pimpact::ScalarField<D3> x2( space );
	Pimpact::VectorField<D3> b ( space );

  auto op = Pimpact::create<Pimpact::GradOp>( space );


	OT nxS = space->ei(Pimpact::F::S,Pimpact::X,Pimpact::B::Y) - space->si(Pimpact::F::S,Pimpact::X,Pimpact::B::Y) + 1;
	OT nyS = space->ei(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y) - space->si(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y) + 1;
	OT nzS = space->ei(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y) - space->si(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y) + 1;
	OT nS = nxS*nyS*nzS;

	OT nxU = space->ei(Pimpact::F::U,Pimpact::X,Pimpact::B::Y) - space->si(Pimpact::F::U,Pimpact::X,Pimpact::B::Y) + 1;
	OT nyU = space->ei(Pimpact::F::U,Pimpact::Y,Pimpact::B::Y) - space->si(Pimpact::F::U,Pimpact::Y,Pimpact::B::Y) + 1;
	OT nzU = space->ei(Pimpact::F::U,Pimpact::Z,Pimpact::B::Y) - space->si(Pimpact::F::U,Pimpact::Z,Pimpact::B::Y) + 1;
	OT nU = nxU*nyU*nzU;

	OT nxV = space->ei(Pimpact::F::V,Pimpact::X,Pimpact::B::Y) - space->si(Pimpact::F::V,Pimpact::X,Pimpact::B::Y) + 1;
	OT nyV = space->ei(Pimpact::F::V,Pimpact::Y,Pimpact::B::Y) - space->si(Pimpact::F::V,Pimpact::Y,Pimpact::B::Y) + 1;
	OT nzV = space->ei(Pimpact::F::V,Pimpact::Z,Pimpact::B::Y) - space->si(Pimpact::F::V,Pimpact::Z,Pimpact::B::Y) + 1;
	OT nV = nxV*nyV*nzV;

	OT nxW = space->ei(Pimpact::F::W,Pimpact::X,Pimpact::B::Y) - space->si(Pimpact::F::W,Pimpact::X,Pimpact::B::Y) + 1;
	OT nyW = space->ei(Pimpact::F::W,Pimpact::Y,Pimpact::B::Y) - space->si(Pimpact::F::W,Pimpact::Y,Pimpact::B::Y) + 1;
	OT nzW = space->ei(Pimpact::F::W,Pimpact::Z,Pimpact::B::Y) - space->si(Pimpact::F::W,Pimpact::Z,Pimpact::B::Y) + 1;
	OT nW = nxW*nyW*nzW;

	Teuchos::Tuple<OT,3> NX = Teuchos::tuple( nxU,nxV,nxW );
	Teuchos::Tuple<OT,3> NY = Teuchos::tuple( nyU,nyV,nyW );
	Teuchos::Tuple<OT,3> NZ = Teuchos::tuple( nzU,nzV,nzW );
	Teuchos::Tuple<OT,3> offset = Teuchos::tuple( 0, nU, nU+nV );

	Teuchos::SerialDenseMatrix<OT,ST> G  ( nU+nV+nW, nS);

	// test G

	Teuchos::RCP<std::ostream> output = Pimpact::createOstream( "GradOp.txt" );
	*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

	for( OT kk=space->si(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y); kk<=space->ei(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y); ++kk )
		for( OT jj=space->si(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y); jj<=space->ei(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y); ++jj )
			for( OT ii=space->si(Pimpact::F::S,Pimpact::X,Pimpact::B::Y); ii<=space->ei(Pimpact::F::S,Pimpact::X,Pimpact::B::Y); ++ii ) {
				OT JJ =       ii-space->si(Pimpact::F::S,Pimpact::X,Pimpact::B::Y)
					+ nxS*(     jj-space->si(Pimpact::F::S,Pimpact::Y,Pimpact::B::Y) )
					+ nxS*nyS*( kk-space->si(Pimpact::F::S,Pimpact::Z,Pimpact::B::Y) );
				x2.init();
				x2(ii,jj,kk) = 1.;
				op->applyG( x2, b );

				for( Pimpact::F fdir=Pimpact::F::U; fdir<D3::sdim; ++fdir ) {
					int dir = static_cast<int>( fdir );
					for( OT k=space->si(fdir,Pimpact::Z,Pimpact::B::Y); k<=space->ei(fdir,Pimpact::Z,Pimpact::B::Y); ++k )
						for( OT j=space->si(fdir,Pimpact::Y,Pimpact::B::Y); j<=space->ei(fdir,Pimpact::Y,Pimpact::B::Y); ++j )
							for( OT i=space->si(fdir,Pimpact::X,Pimpact::B::Y); i<=space->ei(fdir,Pimpact::X,Pimpact::B::Y); ++i ) {
								OT II =               i-space->si(fdir,Pimpact::X,Pimpact::B::Y)
									+ NX[dir]*(         j-space->si(fdir,Pimpact::Y,Pimpact::B::Y) )
									+ NX[dir]*NY[dir]*( k-space->si(fdir,Pimpact::Z,Pimpact::B::Y) ) + offset[dir];
								x.init();
								x(fdir)(i,j,k) = 1.;

								G(II,JJ) = x.dot( b );
							}
				}
			}
	*output << G;

	// test J
	Teuchos::RCP<std::ostream> outputJ = Pimpact::createOstream( "JOp.txt" );
	*output << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) ;

	Teuchos::SerialDenseMatrix<OT,ST> J  ( nU+nV+nW, nU+nV+nW );
	for( Pimpact::F fdir=Pimpact::F::U; fdir<D3::sdim; ++fdir ) {
		int dir = static_cast<int>(fdir);
		for( OT kk=space->si(fdir,Pimpact::Z,Pimpact::B::Y); kk<=space->ei(fdir,Pimpact::Z,Pimpact::B::Y); ++kk )
			for( OT jj=space->si(fdir,Pimpact::Y,Pimpact::B::Y); jj<=space->ei(fdir,Pimpact::Y,Pimpact::B::Y); ++jj )
				for( OT ii=space->si(fdir,Pimpact::X,Pimpact::B::Y); ii<=space->ei(fdir,Pimpact::X,Pimpact::B::Y); ++ii ) {
					OT JJ =           ii-space->si(fdir,Pimpact::X,Pimpact::B::Y)
						+ NX[dir]*(     jj-space->si(fdir,Pimpact::Y,Pimpact::B::Y) )
						+ NX[dir]*NY[dir]*( kk-space->si(fdir,Pimpact::Z,Pimpact::B::Y) )+offset[dir];
					b.init();
					b(fdir)(ii,jj,kk) = 1.;
					op->applyJ( b );

					for( OT k=space->si(fdir,Pimpact::Z,Pimpact::B::Y); k<=space->ei(fdir,Pimpact::Z,Pimpact::B::Y); ++k )
						for( OT j=space->si(fdir,Pimpact::Y,Pimpact::B::Y); j<=space->ei(fdir,Pimpact::Y,Pimpact::B::Y); ++j )
							for( OT i=space->si(fdir,Pimpact::X,Pimpact::B::Y); i<=space->ei(fdir,Pimpact::X,Pimpact::B::Y); ++i ) {
								OT II =               i-space->si(fdir,Pimpact::X,Pimpact::B::Y)
									+ NX[dir]*(         j-space->si(fdir,Pimpact::Y,Pimpact::B::Y) )
									+ NX[dir]*NY[dir]*( k-space->si(fdir,Pimpact::Z,Pimpact::B::Y) ) + offset[dir];
								x.init();
								x(fdir)(i,j,k) = 1.;

								J(II,JJ) = x.dot( b );
							}
				}
	}
	*outputJ << J ;
}
#endif



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, DivGradO2Op, SpaceT ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  Pimpact::ScalarField<SpaceT> x ( space );
  Pimpact::ScalarField<SpaceT> x2( space );
  Pimpact::ScalarField<SpaceT> b ( space );
  Pimpact::ScalarField<SpaceT> diff( space );


  auto op = Pimpact::create<Pimpact::DivGradOp>( space );
  auto op2 = Pimpact::create<Pimpact::DivGradO2Op>( space );

	if( print ) {
		op->print();
		space->getInterpolateV2S()->print();
		op2->print();
	}

  // zero test
  b.init();
  x.random();

  op->apply( b, x );

	ST error = x.norm( Belos::InfNorm );
	if( 0==rank )
		std::cout << "error(zero): " << error << "\n";
  TEST_EQUALITY( error<eps, true );

  b.init();
  x.random();

  op2->apply( b, x );

	error = x.norm( Belos::InfNorm );
	if( 0==rank )
		std::cout << "error(zero) O2: " << error << "\n";
  TEST_EQUALITY( error<eps, true );

	// one test
  b.init( 1. );
  x.random();

  op->apply( b, x );

	error = x.norm( Belos::InfNorm );
	if( 0==rank )
		std::cout << "error(const): " << error << "\n";
	TEST_EQUALITY( error<eps, true );

  b.init( 1. );
  x.random();

  op2->apply( b, x );

	error = x.norm( Belos::InfNorm );
	if( 0==rank )
		std::cout << "error(const) O2: " << error << "\n";
	if( error>= eps ) {
		//x->print();
	}
	TEST_EQUALITY( error<eps, true );

	// consistency test
	b.random();
	op->apply( b, x );
	op2->apply( b, x2 );

	diff.add( 1., x2, -1., x );
	ST errInf = diff.norm( Belos::InfNorm );
	ST err2 = diff.norm( Belos::TwoNorm )/std::sqrt( diff.getLength() );
	if( 0==rank )
		std::cout << "consistency error random: inf: " << errInf << ", two: " << err2 << "\n";
	TEST_EQUALITY( errInf<eps, true );
	TEST_EQUALITY( err2<eps, true );
	if( errInf>=eps && write ) {
		diff.write();
	}

	for( int dir=0; dir<=6; ++dir ) {

		Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>( dir );

		b.initField( type );
		op->apply( b, x );
		op2->apply( b, x2 );

		diff.add( 1., x2, -1., x );
		errInf = diff.norm( Belos::InfNorm );
		err2 = diff.norm( Belos::TwoNorm )/std::sqrt( diff.getLength() );
		if( 0==rank )
			std::cout << "consistency error(" << type << "): inf: " << errInf << ", two: " << err2 << "\n";
		TEST_EQUALITY( errInf<eps, true );
		TEST_EQUALITY( err2<eps, true );
		if( errInf>=eps && err2>=eps ) {
			std::string r = std::to_string( static_cast<long long>( rank ) ); // long long needed on brutus(intel)
			x2.print( *Pimpact::createOstream(
						"error_dgo2_"+Pimpact::toString(type)+"_r"+r+".txt" ) );
			if( write ) {
				x.write(dir+1);
				x2.write( (dir+1)*10 );
				diff.write( (dir+1)*100 );
			}
		}
	}

	// InvDiag consistency test
	x.init( 1. );
	op->applyInvDiag( x, b );
	op2->applyInvDiag( x, x2 );
	if( write ) {
		b.write();
		x2.write(1);
	}

	diff.add( 1., x2, -1., b );
	if( write ) diff.write(2);
	error = diff.norm( Belos::InfNorm );
	if( 0==space->rankST() ) 
		std::cout << "diff InvDiag: " << error << "\n";
	TEST_EQUALITY( error<eps, true );
}

using D22 = Pimpact::Space<ST,OT,2,d,2>;
using D32 = Pimpact::Space<ST,OT,3,d,2>;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Op, D22 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Op, D32 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, DivGradTransposeOp, SpaceT ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  Pimpact::ScalarField<SpaceT> xp( space );
  auto xv = Pimpact::create<Pimpact::VectorField>( space );
  auto bp = Pimpact::create<Pimpact::ScalarField>( space );
  auto bp2 = Pimpact::create<Pimpact::ScalarField>( space );
  auto bv = Pimpact::create<Pimpact::VectorField>( space );
  auto bv2 = Pimpact::create<Pimpact::VectorField>( space );


  auto div = Pimpact::create<Pimpact::DivOp>( space );
  auto grad = Pimpact::create<Pimpact::GradOp>( space );

	auto divGrad = Pimpact::createDivGradOp( div, grad );

	if( print ) {
		div->print();
		grad->print();
	}

	xp.init(1.);
	xv->init(1.);
	////xp->random();

	div->apply(  xp, *bv  );
	std::cout << "||div^T(ones)||" << bv->norm() << "\n";
	grad->apply( xp, *bv2 );
	std::cout << "||grad(ones)||" << bv2->norm() << "\n";

	if( write ) {
		bv->write();
		bv2->write(1);
	}

	grad->apply( *xv, *bp);
	std::cout << "||grad^T(ones)||" << bp->norm() << "\n";
	div->apply(  *xv, *bp2  );
	std::cout << "||div(ones)||" << bp2->norm() << "\n";

	if( write ) {
		bp->write();
		bp2->write(1);
	}

	divGrad->apply( xp, *bp );
	divGrad->apply( xp, *bp2, Belos::TRANS );
	std::cout << "||divGrad(ones): " << bp->norm() << "\n";
	std::cout << "||divGrad^T(ones): " << bp2->norm() << "\n";

	if( write ) {
		bp->write(2);
		bp2->write(3);
	}

	ST pi2 = std::atan(1.)*8.;
	xp.initFromFunction(
			[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
	//xp.random();
	divGrad->apply( xp, *bp );
	divGrad->apply( xp, *bp2, Belos::TRANS );

	if( write ) {
		bp->write(1);
		bp2->write(2);
	}
	bp->add( 1., *bp, -1., *bp2 );
	std::cout << "difference(divgrad, divgrad^T): " << bp->norm() << "\n";

	if( write ) bp->write(3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradTransposeOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradTransposeOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, DivGradO2Smoother, SType ) {

	setParameter( SType::SpaceT::sdim );

	Teuchos::RCP<const typename SType::SpaceT> space = Pimpact::create<typename SType::SpaceT>( pl );

	auto coord = space->getCoordinatesGlobal();
	if( 0==rank ) {
		Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "coord.txt" );
		coord->print( *fstream );
	}

	auto x = Pimpact::create<Pimpact::ScalarField>( space );

  auto b = Pimpact::create<Pimpact::ScalarField>( space );

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );
	//op->print2Mat();

	Teuchos::RCP<Teuchos::ParameterList> ppl = Teuchos::parameterList();

	// compute EV
	ST evMax;
	ST evMin;

	auto ev = 
		Teuchos::rcp( new Pimpact::TeuchosEigenvalues<Pimpact::DivGradO2Op<typename SType::SpaceT> >( op ) );

	if( 0==rank )	std::cout << "\n";
	for( int i=0; i<SType::SpaceT::sdim; ++i ) {
		ev->computeEV( static_cast<Pimpact::ECoord>(i), evMax, evMin );
		if( 0==rank )	
			std::cout << static_cast<Pimpact::ECoord>(i) << ": "
				<< evMax << "\t" <<evMin << "\n";
	}

	ev->computeEV( evMax, evMin );
	if( 0==rank )	
		std::cout << "glob : " << evMax << "\t" <<evMin << "\n";

	if( nx<=20 && ny<=20 && nz<=20 ) {
		ev->computeFullEV( evMax, evMin );
		if( 0==rank )	
			std::cout << "glob2: " << evMax << "\t" <<evMin << "\n";
	}

	ppl->set<int>( "numIters", sweeps );

	// good result with no stretch + own
	ppl->set<ST>( "max EV", evMax );
	ppl->set<ST>( "min EV", evMin*1.1 );
	//ppl->set<ST>( "max EV", evMax*1.  );
	//ppl->set<ST>( "min EV", evMin*2.1 );

	ppl->set<bool>( "with output", true );
	
	// JSmoother
  auto smoother = Pimpact::create<SType>( op, ppl );

	// --- zero rhs test --- 
	x->random();
	b->init( 0. );

	//Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "conv_"+smoother->getLabel()+".txt", 0 );
	if( 0==rank ) {
		if( 0==rank ) {
			std::cout << "\nstep\terror\t\trate\n";
			std::cout << 0 << "\t" << 1. << "\n";
		}
	}
	ST err0 = x->norm( Belos::InfNorm );
	ST errP = err0;

	if( write ) x->write( 0 );
	for( int i=1; i<=nIter; ++i ) {
		smoother->apply( *b, *x );
		if( write ) x->write(i);
		ST err = x->norm( Belos::InfNorm );
		if( 0==rank ) {
				std::cout << i << "\t" << err/err0 << "\t" << err/errP << "\n";
				//fstream << i << "\t" << err/err0 << "\t" << err/errP << "\n";
		}
		errP = err;
	}

	// --- grad test ---
	auto sol = x->clone();
	auto e = x->clone();
	auto res = x->clone();
	for( int dir=1; dir<=3; ++dir ) {

		Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

		x->initField( type );
		sol->initField( type );
		sol->level();

		op->apply( *x, *b );

		//x->init( 0 );
		x->random();

		// residual
		op->apply( *x, *res );
		res->add( -1, *b, 1., *res );
		ST residual = res->norm();
		ST res0 = res->norm();
		ST resP = res0;
		//error
		res->add( 1., *sol, -1., *x );
		ST err0 = res->norm();
		ST errP = err0;
		//ST err0 = 1.;
		//ST errP =1.;

		if( 0==rank ) {
			std::cout << "\n\n\t\t\t--- " << type << " test ---\n";
			std::cout << "\tresidual:\trate:\t\t\terror:\t\trate:\n";
			std::cout << "\t"  << 1.<< "\t\t\t\t" << 1.  << "\n";
		}
		for( int i=0; i<nIter; ++i ) {
			smoother->apply( *b, *x );
			x->level();

			op->apply( *x, *res );
			res->add( -1, *b, 1., *res );
			residual = res->norm();
			res->add( 1., *sol, -1., *x );
			ST err = res->norm();

			if( 0==rank )
				std::cout << "\t" << residual/res0 << "\t" <<  residual/resP << "\t\t" << err/err0 << "\t" <<  err/errP  << "\n";
			resP = residual;
			errP = err;
		}
		TEST_EQUALITY( residual/res0<1.e-1, true );
	}

	// --- consistency test ---
	for( int dir=1; dir<=6; ++dir ) {

		x->initField( static_cast<Pimpact::EScalarField>(dir) );
		x->level();
		Pimpact::ScalarField<typename SType::SpaceT> xp(space);
		xp = *x;

		op->apply( *x, *b );

		smoother->apply( *b, *x );
		x->level();

		xp.add( 1., xp, -1., *x );
		ST err2 = xp.norm( Belos::InfNorm )/std::sqrt( static_cast<ST>(xp.getLength()) );
		ST errInf = xp.norm( Belos::InfNorm );

		if( 0==rank ) {
			std::cout << "consistency for " << dir << ": ||" << err2 << "||_2, ||" << errInf << "||_inf\n";
		}

		TEST_EQUALITY( err2<eps, true );
		TEST_EQUALITY( errInf<eps, true );
		if( err2>eps && errInf>eps ) 
			if( write ) x->write( nIter + dir + 1 );

	}

}

using JT2D   = Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op<D2> >;
using CheT2D = Pimpact::Chebyshev< Pimpact::DivGradO2Op<D2> >;
using LT2D   = Pimpact::DivGradO2LSmoother< Pimpact::DivGradO2Op<D2> >;

using JT3D   = Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op<D3> >;
using CheT3D = Pimpact::Chebyshev< Pimpact::DivGradO2Op<D3> >;
using LT3D   = Pimpact::DivGradO2LSmoother< Pimpact::DivGradO2Op<D3> >;

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, JT2D   )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, CheT2D )
//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, LT2D   ) // debug this

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, JT3D   )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, CheT3D )
//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Smoother, LT3D   ) // debug this




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, DivGradO2Inv, SpaceT ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

	Pimpact::ScalarField<SpaceT> x ( space );
	Pimpact::ScalarField<SpaceT> xp( space );
	Pimpact::ScalarField<SpaceT> b ( space );

  auto op = Pimpact::create<Pimpact::DivGradO2Op>( space );

	if( print )
		op->print2Mat();

  auto ppl = Teuchos::parameterList();

  auto solver = Pimpact::create<Pimpact::DivGradO2Inv>( op, ppl );

  // --- consistency test ---
	for( int field=0; field<=(domain!=1?6:0); ++field ) {
		x.initField( static_cast<Pimpact::EScalarField>(field) );
		x.level();
		xp.random();

		op->apply( x, b );
		solver->apply( b, xp );

		xp.add( 1., xp, -1., x );
		xp.level();

		ST err2 = xp.norm( Belos::InfNorm )/std::sqrt( static_cast<ST>(xp.getLength()) );
		ST errInf = xp.norm( Belos::InfNorm );
		if( errInf>=eps ) {
			if( write ) xp.write();
			if( print ) xp.print();
			//std::cout << "rank: " << rank << "\te: " << xp.norm( Belos::InfNorm, false ) << "\n";
			//if( 0==rank )
				//xp.print();
		}


		if( 0==rank )
			std::cout << "consistency for " << static_cast<Pimpact::EScalarField>(field) << ": ||" << err2 << "||_2, ||" << errInf << "||_inf\n";

		if( npx*npy*npz==1 ) {
			TEST_EQUALITY( err2<eps, true );
			TEST_EQUALITY( errInf<eps, true );
		}
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Inv, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradO2Inv, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, ForcingOp, SpaceT ) {

	setParameter( SpaceT::sdim );

	auto space = Pimpact::create<SpaceT>( pl );

	auto vel = Pimpact::create<Pimpact::VectorField>( space );
	auto force = Pimpact::create<Pimpact::VectorField>( space );
	auto res = Pimpact::create<Pimpact::VectorField>( space );

	vel->init(1.);
	//force->initField( Pimpact::EVectorField(11) );
	force->random();
	res->random();

	auto op = Pimpact::createForcingOp( force, 1. );

	op->apply( *vel, *res );
	vel->add( 1., *res, -1., *force );

	TEST_EQUALITY( vel->norm( Belos::InfNorm )<eps, true );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ForcingOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, ForcingOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiOperator, Add2Op, SpaceT ) {

	setParameter( SpaceT::sdim );

	auto space = Pimpact::create<SpaceT>( pl );

	auto velx = Pimpact::create<Pimpact::VectorField>(space);

	auto x = Teuchos::rcp( new Pimpact::MultiField<Pimpact::VectorField<SpaceT> >(space,10) );
	auto y = Teuchos::rcp( new Pimpact::MultiField<Pimpact::VectorField<SpaceT> >(space,10) );

	auto op = Pimpact::createOperatorBase(
			Pimpact::createMultiOpWrap(
				Pimpact::createAdd2Op(
					Pimpact::create<Pimpact::HelmholtzOp>( space ),
					Pimpact::create<Pimpact::HelmholtzOp>( space )
					)
				)
			);

	x->random();
	op->apply( *x, *y);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiOperator, Add2Op, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiOperator, Add2Op, D3 )


template<class T> using WUP = Pimpact::MultiOpUnWrap<Pimpact::MultiOpWrap<T> >;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiOperator, MulitOpUnWrap, SpaceT ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  auto x = Pimpact::create<Pimpact::VectorField>(space);
  auto y = Pimpact::create<Pimpact::VectorField>(space);

  auto op = Pimpact::createOperatorBase(
			Pimpact::createMultiOpWrap(
				Pimpact::create<Pimpact::HelmholtzOp>( space)  ) );

  auto opUW = Pimpact::create<Pimpact::MultiOpUnWrap>( op );

  auto opUW2 =
      Pimpact::create<WUP>(
          Pimpact::create<Pimpact::HelmholtzOp>(space)
      );


  //x->initField(Pimpact::RankineVortexD2 );
	x->random();

  opUW2->apply( *x, *y );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiOperator, MulitOpUnWrap, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiOperator, MulitOpUnWrap, D3 )



TEUCHOS_UNIT_TEST( MultiHarmonicOperator, MultiHarmonicOpWrap ) {

	using SpaceT = Pimpact::Space<ST,OT,3,4,dNC>;

	Pimpact::setBoundaryConditions( pl, domain );

	pl->set<bool>( "spectral in time", true );

	pl->set( "lx", lx );
	pl->set( "ly", ly );
	pl->set( "lz", lz );

	//  grid size
	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	// grid stretching
	setStretching();

  // processor grid size
  pl->set( "npx", npx );
  pl->set( "npy", npy );
  pl->set( "npz", npz );
  pl->set( "npf", npf );

  auto space = Pimpact::create<Pimpact::Space<ST,OT,3,4,dNC> >( pl );


  auto op = Pimpact::createMultiHarmonicOpWrap< Pimpact::HelmholtzOp<SpaceT> >(
      Pimpact::create<Pimpact::HelmholtzOp>(space) );

  Pimpact::MultiHarmonicField< Pimpact::VectorField<Pimpact::Space<ST,OT,3,4,dNC> > > x( space );
  Pimpact::MultiHarmonicField< Pimpact::VectorField<Pimpact::Space<ST,OT,3,4,dNC> > > y( space );

  op->apply( x, y );
}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, DivOp, SpaceT ) {

	setParameter( SpaceT::sdim );

	//  grid size
	for( int dir=0; dir<SpaceT::sdim; ++dir ) {

		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );

		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> dofs( ns );

		ST pi2 = 4.*std::acos(0.);

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

			Pimpact::VectorField<SpaceT> vel( space );
			Pimpact::ScalarField<SpaceT> p  ( space );
			Pimpact::ScalarField<SpaceT> sol( space );

			// init 
			if( 0==dir ) {
				vel(Pimpact::F::U).initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
				sol.initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST { return( -std::sin(x*pi2)*pi2 ); } );
			}
			else if( 1==dir ) {
				vel(Pimpact::F::V).initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol.initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST { return( -std::sin(y*pi2)*pi2 ); } );
			}
			else if( 2==dir ) {
				vel(Pimpact::F::W).initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
				sol.initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST { return( -std::sin(z*pi2)*pi2 ); } );
			}


			if( print ) vel(Pimpact::F::W).print();
			auto op = Pimpact::create<Pimpact::DivOp>( space );

			op->apply( vel, p );

			// compute error
			p.add( 1., sol, -1., p );
			if( write ) p.write(n);
			error2[n] = std::log10( p.norm( Belos::TwoNorm ) / sol.norm( Belos::TwoNorm ) );
			errorInf[n] = std::log10( p.norm( Belos::InfNorm ) / sol.norm( Belos::InfNorm ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
			if( 0==rank )	
				std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";

		}
		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << "DivOp: order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << "DivOp: order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";
		// test
		TEST_EQUALITY( -order2>3., true );
		TEST_EQUALITY( -orderInf>3., true );
	}
	
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, DivOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, DivOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, InterpolateV2SOp, SpaceT ) { 

	int rank = 0;
	setParameter( SpaceT::sdim );

	for( Pimpact::F dir=Pimpact::F::U; dir<SpaceT::sdim; ++dir ) {

		//  grid size
		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );

		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> dofs( ns );

		ST pi2 = std::atan(1)*8;

		if( 0==rank )	std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

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

			Pimpact::VectorField<SpaceT> vel( space );
			Pimpact::ScalarField<SpaceT> p  ( space );
			Pimpact::ScalarField<SpaceT> sol( space );

			// init 
			if( 0==dir ) {
				vel(Pimpact::F::U).initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
				sol.initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
			}
			else if( 1==dir ) {
				vel(Pimpact::F::V).initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol.initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
			}
			else if( 2==dir ) {
				vel(Pimpact::F::W).initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
				sol.initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
			}

			auto op = Pimpact::createInterpolateV2S( space );

			op->apply( vel(dir), p );

			// compute error
			p.add( 1., sol, -1., p );

			if( write ) p.write( n );
			if( print ) p.print();

			error2[n]   = std::log10( p.norm( Belos::TwoNorm ) / sol.norm( Belos::TwoNorm ) );
			errorInf[n] = std::log10( p.norm( Belos::InfNorm ) / sol.norm( Belos::InfNorm ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );

			if( 0==rank )	
				std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";

		}
		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << "InterpolateV2SOp: order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << "InterpolateV2SOp: order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";
		// test
		TEST_EQUALITY( -order2  >4., true );
		TEST_EQUALITY( -orderInf>4., true );
	}
	
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, InterpolateV2SOp, D2 ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, InterpolateV2SOp, D3 ) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, InterpolateS2VOp, SpaceT ) { 

	setParameter( SpaceT::sdim );

	for( Pimpact::F dir=Pimpact::F::U; dir<SpaceT::sdim; ++dir ) {

		//  grid size
		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );


		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> dofs( ns );

		ST pi2 = std::atan(1)*8;

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

			Pimpact::VectorField<SpaceT> vel( space );
			Pimpact::ScalarField<SpaceT> p  ( space );
			Pimpact::VectorField<SpaceT> sol( space );

			// init 
			if( 0==dir ) {
				p.initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
				sol(Pimpact::F::U).initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
			}
			else if( 1==dir ) {
				p.initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol(Pimpact::F::V).initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
			}
			else if( 2==dir ) {
				p.initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
				sol(Pimpact::F::W).initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
			}

			auto op = Pimpact::create<Pimpact::InterpolateS2V>( space );

			op->apply( p, vel(dir) );

			// compute error
			vel(dir).add( 1., sol(dir), -1., vel(dir) );
			if( write ) vel.write(n);

			error2[n]   = std::log10( vel(dir).norm( Belos::TwoNorm ) / sol(dir).norm( Belos::TwoNorm ) );
			errorInf[n] = std::log10( vel(dir).norm( Belos::InfNorm ) / sol(dir).norm( Belos::InfNorm ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
			if( 0==rank )	
				std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";

		}
		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << "InterpolateS2V: order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << "InterpolateS2V: order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";
		// test
		TEST_EQUALITY( -order2  >4., true );
		TEST_EQUALITY( -orderInf>3.9, true );
	}
	
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, InterpolateS2VOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, InterpolateS2VOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, extrapolateBC, SpaceT ) { 

	if( domain!=1 ) {

		setParameter( SpaceT::sdim );

		for( Pimpact::F dir=Pimpact::F::U; dir<SpaceT::sdim; ++dir ) {

			//  grid size
			pl->set<OT>( "nx", 9 );
			pl->set<OT>( "ny", 9 );
			pl->set<OT>( "nz", 9 );
			pl->set<OT>( "nf", 1 );

			std::vector<ST> error2( ns );
			std::vector<ST> errorInf( ns );
			std::vector<ST> dofs( ns );

			ST pi2 = std::atan(1)*8;

			if( 0==rank )	
				std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

			for( OT n=0; n<ns; ++n ) {

				if( Pimpact::F::U==dir ) 
					pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
				else if( Pimpact::F::V==dir )
					pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
				else if( Pimpact::F::W==dir )
					pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

				// grid stretching
				setStretching();

				auto space = Pimpact::create<SpaceT>( pl );

				Pimpact::VectorField<SpaceT> vel( space );
				Pimpact::VectorField<SpaceT> sol( space );

				// init 
				if( Pimpact::F::U==dir ) {
					vel(Pimpact::F::U).initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::sin(x*pi2) ); } );
					sol(Pimpact::F::U).initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::sin(x*pi2) ); } );
				}
				else if( Pimpact::F::V==dir ) {
					vel(Pimpact::F::V).initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return( std::sin(y*pi2) ); } );
					sol(Pimpact::F::V).initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return( std::sin(y*pi2) ); } );
				}
				else if( Pimpact::F::W==dir ) {
					vel(Pimpact::F::W).initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::sin(z*pi2) ); } );
					sol(Pimpact::F::W).initFromFunction(
							[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::sin(z*pi2) ); } );
				}
				vel.extrapolateBC();

				// compute error
				vel(dir).add( 1., sol(dir), -1., vel(dir), Pimpact::B::Y );
				if( write ) vel.write( n );

				error2[n]   = std::log10( vel(dir).norm( Belos::TwoNorm ) / sol(dir).norm( Belos::TwoNorm ) );
				errorInf[n] = std::log10( vel(dir).norm( Belos::InfNorm ) / sol(dir).norm( Belos::InfNorm ) );

				dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
				if( 0==rank )	
					std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";
			}
			// compute order
			ST order2 = order<ST>( dofs, error2 );
			if( 0==rank )	
				std::cout << "extrapolateBC: order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

			ST orderInf = order<ST>( dofs, errorInf );
			if( 0==rank )	
				std::cout << "extrapolateBC: order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";
			// test
			TEST_EQUALITY( -order2  >3., true );
			TEST_EQUALITY( -orderInf>3., true );
		}
	}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, extrapolateBC, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, extrapolateBC, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, GradOp, SpaceT ) { 

	setParameter( SpaceT::sdim );

	for( Pimpact::F dir=Pimpact::F::U; dir<SpaceT::sdim; ++dir ) {

		//  grid size
		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );

		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> dofs( ns );

		ST pi2 = std::atan(1)*8;

		if( 0==rank )	
			std::cout << "\n\nN\t||e||_2\t\t||e||_inf\n";

		for( OT n=0; n<ns; ++n ) {

			if( Pimpact::F::U==dir ) 
				pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
			else if( Pimpact::F::V==dir )
				pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
			else if( Pimpact::F::W==dir )
				pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

			setStretching();

			auto space = Pimpact::create<SpaceT>( pl );

			Pimpact::VectorField<SpaceT> vel( space );
			auto p   = Pimpact::create<Pimpact::ScalarField>( space );
			Pimpact::VectorField<SpaceT> sol( space );

			// init 
			if( Pimpact::F::U==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST { return(  std::cos(x*pi2) ); } );
				sol(Pimpact::F::U).initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {
							if( (std::abs( y )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
									(std::abs( y-ly )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
							    (std::abs( z )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
									(std::abs( z-lz )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 )  )
								return( -std::sin(x*pi2)*pi2*0.1 ); 
							else
								return( -std::sin(x*pi2)*pi2 );
						} );
			}
			else if( Pimpact::F::V==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol(Pimpact::F::V).initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {
							if( (std::abs( x )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
									(std::abs( x-lx )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
							    (std::abs( z )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
									(std::abs( z-lz )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 )  )
								return( -std::sin(y*pi2)*pi2*0.1 ); 
							else
								return( -std::sin(y*pi2)*pi2 );
						} );
			}
			else if( Pimpact::F::W==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return(  std::cos(z*pi2) ); } );
				sol(Pimpact::F::W).initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST { 
							if( (std::abs( x )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
									(std::abs( x-lx )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
							    (std::abs( y )   <Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
									(std::abs( y-ly )<Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 )  )
								return( -std::sin(z*pi2)*pi2*0.1 ); 
							else
								return( -std::sin(z*pi2)*pi2 );
						}  );
			}

			auto op = Pimpact::create<Pimpact::GradOp>( space );

			op->apply(  *p, vel );

			// compute error
			if( write ) sol.write( n );
			vel(dir).add( 1., sol(dir), -1., vel(dir), Pimpact::B::Y );
			
			error2[n]   = std::log10( vel(dir).norm( Belos::TwoNorm, Pimpact::B::Y ) / sol(dir).norm( Belos::TwoNorm, Pimpact::B::Y ) );
			errorInf[n] = std::log10( vel(dir).norm( Belos::InfNorm, Pimpact::B::Y ) / sol(dir).norm( Belos::InfNorm, Pimpact::B::Y ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
			if( 0==rank )	
				std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";

		}
		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << "GradOp: order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << "GradOp: order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";
		// test
		TEST_EQUALITY( -order2  >3., true );
		TEST_EQUALITY( -orderInf>3., true );
	}
	
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, GradOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, GradOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, DivGradOp, OperatorT ) { 

	ST pi2 = std::atan(1)*8;
	std::string label;

	setParameter( OperatorT::SpaceT::sdim );

	for( int dir=0; dir<OperatorT::SpaceT::sdim; ++dir ) {

		//  grid size
		pl->set<OT>( "nx", 9 );
		pl->set<OT>( "ny", 9 );
		pl->set<OT>( "nz", 9 );
		pl->set<OT>( "nf", 1 );

		// grid stretching
		setStretching();

		std::vector<ST> error2( ns );
		std::vector<ST> errorInf( ns );
		std::vector<ST> error2Ref( ns );
		std::vector<ST> errorInfRef( ns );
		std::vector<ST> dofs( ns );

		// reference
		Teuchos::RCP<Teuchos::ParameterList> plRef = Teuchos::parameterList();
		*plRef = *pl;

		if( 0==dir ) 
			plRef->set<OT>( "nx", 8*std::pow(2,ns)+1 );
		else if( 1==dir )
			plRef->set<OT>( "ny", 8*std::pow(2,ns)+1 );
		else if( 2==dir )
			plRef->set<OT>( "nz", 8*std::pow(2,ns)+1 );

		// grid stretching
		setStretching();

		auto spaceRef = Pimpact::create<typename OperatorT::SpaceT>( plRef );

		auto yRef = Pimpact::create<Pimpact::ScalarField>( spaceRef );
		auto xRef   = Pimpact::create<Pimpact::ScalarField>( spaceRef );

		auto opRef = Pimpact::create<OperatorT>( spaceRef );
		//opRef->print();
		//auto opRef = Pimpact::create<Pimpact::DivGradOp>( spaceRef );

		// init 
		if( 0==dir ) {
			xRef->initFromFunction(
					[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
		}
		else if( 1==dir ) {
			xRef->initFromFunction(
					[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
		}
		else if( 2==dir ) {
			xRef->initFromFunction(
					[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
		}


		yRef->random();
		opRef->apply( *xRef, *yRef );
		if( 0==rank )	
			std::cout << "\n\n"<<std::setw(5)<< "N\t||e||_2\t\t||e||_inf\n";

		// loop of discretization
		for( OT n=0; n<ns; ++n ) {

			pl->set<OT>( "nx", 9 );
			pl->set<OT>( "ny", 9 );
			pl->set<OT>( "nz", 9 );

			if( 0==dir ) 
				pl->set<OT>( "nx", 8*std::pow(2,n)+1 );
			else if( 1==dir )
				pl->set<OT>( "ny", 8*std::pow(2,n)+1 );
			else if( 2==dir )
				pl->set<OT>( "nz", 8*std::pow(2,n)+1 );

			// grid stretching
			setStretching();

			auto space = Pimpact::create<typename OperatorT::SpaceT>( pl );

			auto y = Pimpact::create<Pimpact::ScalarField>( space );
			auto p   = Pimpact::create<Pimpact::ScalarField>( space );
			auto e   = Pimpact::create<Pimpact::ScalarField>( space );
			auto eRef= Pimpact::create<Pimpact::ScalarField>( space );
			auto sol = y->clone( Pimpact::ECopy::Shallow );

			// init 
			if( 0==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(x*pi2) ); } );
				sol->initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {  
							if( (std::abs( y )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
									(std::abs( y-ly )<Teuchos::ScalarTraits<ST>::eps()  && space->bcu(1)>0 ) ||
							    (std::abs( z )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
									(std::abs( z-lz )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 )  )
								return( - pi2*pi2*0.1*std::cos(x*pi2) );
							else
								return( - pi2*pi2*std::cos(x*pi2) );
						}
					);
			}
			else if( 1==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(y*pi2) ); } );
				sol->initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {
							if( (std::abs( x )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
									(std::abs( x-lx )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
							    (std::abs( z )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
									(std::abs( z-lz )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 )  )
								return( - pi2*pi2*0.1*std::cos(y*pi2) );
							else
								return( - pi2*pi2*std::cos(y*pi2) );
						}
					);
			}
			else if( 2==dir ) {
				p->initFromFunction(
						[&pi2]( ST x, ST y, ST z ) ->ST {  return( std::cos(z*pi2) ); } );
				sol->initFromFunction(
						[&pi2,&space]( ST x, ST y, ST z ) ->ST {  
							if( (std::abs( x )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
									(std::abs( x-lx )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
							    (std::abs( y )   < Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
									(std::abs( y-ly )< Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 )  )
								return( - pi2*pi2*0.1*std::cos(z*pi2) );
							else
								return( - pi2*pi2*std::cos(z*pi2) );
						}
					);
			}

			auto op = Pimpact::create<OperatorT>( space );
			label = op->getLabel();

			y->random();
			op->apply( *p, *y );

			// compute error
			e->add( 1., *sol, -1., *y );
			error2[n] = std::log10( e->norm( Belos::TwoNorm ) / sol->norm( Belos::TwoNorm ) );
			errorInf[n] = std::log10( e->norm( Belos::InfNorm ) / sol->norm( Belos::InfNorm ) );
			dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
			eRef->add( 1., *yRef, -1., *y );
			error2Ref[n] = std::log10( eRef->norm( Belos::TwoNorm ) / std::sqrt( dofs[n] ) );
			errorInfRef[n] = std::log10( eRef->norm( Belos::InfNorm ) );
			if( 0==rank )	 {
				std::cout  << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\t";
				std::cout                                  << std::pow(10.,error2Ref[n]) << "\t" << std::pow(10.,errorInfRef[n]) << "\n";
			}
			if( print ) e->print();
			if( write ) {
				y->write(0);
				sol->write(1);
				e->write(2);
				eRef->write(3);
				yRef->write(4);
			}
		}

		// compute order
		ST order2 = order<ST>( dofs, error2 );
		if( 0==rank )	
			std::cout << label << ": order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

		ST orderInf = order<ST>( dofs, errorInf );
		if( 0==rank )	
			std::cout << label << ": order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";

		ST order2Ref = order<ST>( dofs, error2Ref );
		if( 0==rank )	
			std::cout << label << ": ref order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2Ref << "\n";

		ST orderInfRef = order<ST>( dofs, errorInfRef );
		if( 0==rank )	
			std::cout << label << ": ref order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInfRef << "\n";
		// test
		TEST_EQUALITY( -order2  >2., true );
		TEST_EQUALITY( -orderInf>2., true );
	}
	
}

using DivGradOpT2D = Pimpact::DivGradOp<D2>;
using DivGradO2OpT2D = Pimpact::DivGradO2Op<D2>;

using DivGradOpT3D = Pimpact::DivGradOp<D3>;
using DivGradO2OpT3D = Pimpact::DivGradO2Op<D3>;

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, DivGradOp, DivGradOpT2D   )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, DivGradOp, DivGradO2OpT2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, DivGradOp, DivGradOpT3D   )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, DivGradOp, DivGradO2OpT3D )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Convergence, HelmholtzOp, SpaceT ) { 

	int rank = 0;

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

				Pimpact::VectorField<SpaceT> x( space );
				Pimpact::VectorField<SpaceT> y( space );
				Pimpact::VectorField<SpaceT> sol( space );

				// init 
				ST pi2 = std::atan(1)*8;
				ST ire = 1./space->getDomainSize()->getRe();

				x.init();
				sol.init();
				auto inifunc = [&pi2]( ST x ) -> ST {  
					return(  std::cos(x*pi2) );
				};
				auto solfunc = [&pi2,&ire,&space]( ST x ) -> ST {
					return( std::cos(x*pi2)*pi2*pi2*ire );
				};
 
				if( 0==dir ) {
					x(field).initFromFunction(
							[&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(x) ); } );
					sol(field).initFromFunction(
							[&space,&inifunc,&solfunc]( ST x, ST y, ST z ) ->ST {

								if( space->bcl(Pimpact::X)>0 && x<=Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::X) )
										return( inifunc(std::min(std::max(x,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcl(Pimpact::X) )
										return( 0. );
								}
								else if( space->bcu(Pimpact::X)>0 && x>=1.-Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::X) )
										return( inifunc(std::min(std::max(x,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcu(Pimpact::X) )
										return( 0. );
								}
								else if( space->bcl(Pimpact::Y)>0 && y<=Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::Y) )
										return( inifunc(std::min(std::max(x,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcl(Pimpact::Y) )
										return( 0. );
								}
								else if( space->bcu(Pimpact::Y)>0 && y>=1.-Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::Y) )
										return( inifunc(std::min(std::max(x,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcu(Pimpact::Y) )
										return( 0. );
								}
								else if( space->bcl(Pimpact::Z)>0 && z<=Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::Z) )
										return( inifunc(std::min(std::max(x,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcl(Pimpact::Z) )
										return( 0. );
								}
								else if( space->bcu(Pimpact::Z)>0 && z>=1.-Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::Z) )
										return( inifunc(std::min(std::max(x,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcu(Pimpact::Z) )
										return( 0. );
								}
								return( solfunc(x) ); } );
				}
				else if( 1==dir ) {
					x(field).initFromFunction( [&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(y) ); } );
					sol(field).initFromFunction(
							[&space,&inifunc,&solfunc]( ST x, ST y, ST z ) ->ST {
								if( space->bcl(Pimpact::X)>0 && x<=Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::X) )
										return( inifunc(std::min(std::max(y,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcl(Pimpact::X) )
										return( 0. );
								}
								else if( space->bcu(Pimpact::X)>0 && x>=1.-Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::X) )
										return( inifunc(std::min(std::max(y,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcu(Pimpact::X) )
										return( 0. );
								}
								else if( space->bcl(Pimpact::Y)>0 && y<=Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::Y) )
										return( inifunc(std::min(std::max(y,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcl(Pimpact::Y) )
										return( 0. );
								}
								else if( space->bcu(Pimpact::Y)>0 && y>=1.-Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::Y) )
										return( inifunc(std::min(std::max(y,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcu(Pimpact::Y) )
										return( 0. );
								}
								else if( space->bcl(Pimpact::Z)>0 && z<=Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::Z) )
										return( inifunc(std::min(std::max(y,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcl(Pimpact::Z) )
										return( 0. );
								}
								else if( space->bcu(Pimpact::Z)>0 && z>=1.-Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::Z) )
										return( inifunc(std::min(std::max(y,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcu(Pimpact::Z) )
										return( 0. );
								}
								return( solfunc(y) ); } );
				}
				else if( 2==dir ) {
					x(field).initFromFunction(
							[&inifunc]( ST x, ST y, ST z ) ->ST {  
								return( inifunc(z) ); } );
					sol(field).initFromFunction(
							[&space,&inifunc,&solfunc]( ST x, ST y, ST z ) ->ST {
								if( space->bcl(Pimpact::X)>0 && x<=Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::X) )
										return( inifunc(std::min(std::max(z,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcl(Pimpact::X) )
										return( 0. );
								}
								else if( space->bcu(Pimpact::X)>0 && x>=1.-Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::X) )
										return( inifunc(std::min(std::max(z,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcu(Pimpact::X) )
										return( 0. );
								}
								else if( space->bcl(Pimpact::Y)>0 && y<=Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::Y) )
										return( inifunc(std::min(std::max(z,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcl(Pimpact::Y) )
										return( 0. );
								}
								else if( space->bcu(Pimpact::Y)>0 && y>=1.-Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::Y) )
										return( inifunc(std::min(std::max(z,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcu(Pimpact::Y) )
										return( 0. );
								}
								else if( space->bcl(Pimpact::Z)>0 && z<=Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcl(Pimpact::Z) )
										return( inifunc(std::min(std::max(z,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcl(Pimpact::Z) )
										return( 0. );
								}
								else if( space->bcu(Pimpact::Z)>0 && z>=1.-Teuchos::ScalarTraits<ST>::eps() ) {
									if( Pimpact::BC::Dirichlet==space->bcu(Pimpact::Z) )
										return( inifunc(std::min(std::max(z,0.),1.)) );
									else if( Pimpact::BC::Neumann==space->bcu(Pimpact::Z) )
										return( 0. );
								}
								return( solfunc(z) ); } );
				}

				auto op = Pimpact::create<Pimpact::HelmholtzOp>( space );

				op->apply( x, y );

				// compute error
				if( print ) y(field).print();
				y(field).add( 1., sol(field), -1., y(field) );
				//if( print ) sol->print();
				if( write ) y.write( n );
				error2[n]   = std::log10( y(field).norm( Belos::TwoNorm ) / sol(field).norm( Belos::TwoNorm ) );
				errorInf[n] = std::log10( y(field).norm( Belos::InfNorm ) / sol(field).norm( Belos::InfNorm ) );
				dofs[n] = std::log10( 8.*std::pow(2.,n)+1. );
				if( 0==rank )	
					std::cout << std::pow(10.,dofs[n]) << "\t" << std::pow(10.,error2[n]) << "\t" << std::pow(10.,errorInf[n]) << "\n";
			}
			// compute order
			ST order2 = order<ST>( dofs, error2 );
			if( 0==rank )	
				std::cout << "Helmholtz(" << field << "): order two norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << order2 << "\n";

			ST orderInf = order<ST>( dofs, errorInf );
			if( 0==rank )	
				std::cout << "Helmholtz(" << field << "): order inf norm in "<< static_cast<Pimpact::ECoord>(dir) << "-dir: " << orderInf << "\n";
			// test
			TEST_EQUALITY( -order2  >3., true );
			TEST_EQUALITY( -orderInf>3., true );
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, HelmholtzOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Convergence, HelmholtzOp, D3 )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BasicOperator, DivGradNullSpace, SpaceT ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

	auto null   = Pimpact::create<Pimpact::ScalarField>( space );
	auto res   = Pimpact::create<Pimpact::ScalarField>( space );
	null->random();

	auto op = Pimpact::create<Pimpact::DivOp>( space );

	if( print )
		op->print();

	// DJG
	Pimpact::DivGradNullSpace<Pimpact::DivOp<SpaceT> > compNullspace;// = Pimpact::create<Pimpact::DivGradNullSpace>( );

	compNullspace.computeNullSpace( op, *null );

	auto divGrad = Pimpact::create<Pimpact::DivGradOp>( space );
	if( write )
		null->write( 0 );
	divGrad->apply( *null, *res, Belos::TRANS );
	ST error = res->norm();
	std::cout << "\nerror: " << error << "\n";

	TEST_EQUALITY( error<eps, true );

	// full sys

	compNullspace.computeNullSpace( op, *null, false );

	auto resV   = Pimpact::create<Pimpact::VectorField>( space );
	if( write )
		null->write( 1 );
	op->apply( *null, *resV );
	if( write ) {
		Pimpact::VectorField<SpaceT> resV2( space );
		resV2.abs( *resV );
		resV2.write();
	}
	if( print )
		resV->print();
	ST errorV = resV->norm();
	std::cout << "\nerrorV: " << errorV << "\n";
	//TEST_EQUALITY( errorV<eps, true ); // not working yet as D^T is including BC
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradNullSpace, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( BasicOperator, DivGradNullSpace, D3 )

} // namespace
