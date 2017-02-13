#include <iostream>
#include <cmath>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Fields.hpp"

#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_MultiGrid.hpp"
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

template<class T1,class T2>
using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T>
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T>
using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

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

template<class T> using POP2  = Pimpact::PrecInverseOp< T, ConvDiffJT >;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ModeSolver, ModeNonlinearOp, SpaceT ) {

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
	Pimpact::ModeField< Pimpact::VectorField<SpaceT> > rhs( space );
	Pimpact::ModeField< Pimpact::VectorField<SpaceT> > sol( space );
	Pimpact::ModeField< Pimpact::VectorField<SpaceT> > err( space );

	Teuchos::RCP< Teuchos::ParameterList > param = Teuchos::parameterList();// = Pimpact::createLinSolverParameter( solvName, 1.e-6 );

	std::string solvName = "TFQMR";
	param->get<std::string>( "Solver name", solvName );

	param->sublist("Solver").set( "Output Style", Belos::Brief );
	param->sublist("Solver").set( "Output Frequency", -1 );
	param->sublist("Solver").set( "Verbosity",			        
			Belos::Errors +
			Belos::Warnings +
			Belos::IterationDetails +
			Belos::OrthoDetails +
			Belos::FinalSummary +
			//Belos::TimingDetails +
			Belos::StatusTestDetails +
			Belos::Debug );
	param->sublist("Solver").set( "Maximum Iterations", nMax );
	param->sublist("Solver").set( "Convergence Tolerance", 1.e-6 );
	//param->set( "Flexible Gmres", false );
	
	auto zeroOp = Pimpact::create<ConvDiffOpT>( space );
	
	//Teuchos::RCP<const Pimpact::ModeNonlinearOp<ConvDiffOpT<SpaceT> > op =
	auto op =
		Teuchos::rcp( new Pimpact::ModeNonlinearOp< ConvDiffOpT<SpaceT> >(zeroOp));

	auto modeInv =
		Pimpact::createMultiOpUnWrap( Pimpact::createInverseOp( op, param ) );


	using FSpaceT = SpaceT;
	using CSpaceT = Pimpact::Space<ST,OT,SpaceT::sdim,SpaceT::dimension,2>;

	using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;
	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, 10 );
	auto zeroInv =
		Pimpact::createMultiOpUnWrap(
				Pimpact::createInverseOp( zeroOp, param ) );

	auto mgConvDiff =
		Pimpact::createMultiGrid<
		Pimpact::VectorField,
		TransVF,
		RestrVF,
		InterVF,
		ConvDiffOpT,
		ConvDiffOpT,
		//ConvDiffSORT,
		ConvDiffJT,
		//ConvDiffJT
		POP2
			//POP3
			> ( mgSpaces, Teuchos::parameterList() ) ;

	//if( 0==space->rankST() )
	//mgConvDiff->print();

	zeroInv->getOperatorPtr()->setRightPrec( Pimpact::createMultiOperatorBase(mgConvDiff) );
	modeInv->getOperatorPtr()->setRightPrec( Pimpact::createMultiOperatorBase( Pimpact::create<Pimpact::EddyPrec>(zeroInv) ) );
	//modeInv->getOperatorPtr()->setLeftPrec( Pimpact::createMultiOperatorBase( Pimpact::create<Pimpact::EddyPrec>(zeroInv) ) );

	ST iRe = 1./space->getDomainSize()->getRe();
	ST a2 = space->getDomainSize()->getAlpha2()*iRe;

	//std::cout << "a2: " << a2 << "\n";
	//std::cout << "iRe: " << iRe << "\n";
	// computing zero mode of z
	// set paramteters
	auto para = Teuchos::parameterList();
	para->set<ST>( "mulI", a2  );
	para->set<ST>( "mulC", 1.  );
	para->set<ST>( "mulL", iRe );
	op->setParameter( para );
	modeInv->setParameter( para );

	// initializtion
	{
		Pimpact::VectorField<SpaceT> wind( space );
		wind(Pimpact::F::U).init( 1. );
		wind(Pimpact::F::V).init( 1. );
		////wind(Pimpact::F::W).init( 1. );

		zeroOp->assignField( wind );
		zeroInv->assignField( wind );
	}

	auto initFunC = []( ST x, ST y ) ->ST { return( std::pow((y-0.5),2) ); };
	auto initFunS = []( ST x, ST y ) ->ST { return( std::pow((x-0.5),1) ); };
	auto deriFunC = [=]( ST y ) ->ST { return( 2.*(y-0.5)/ly - iRe*2./ly/ly ); };
	auto deriFunS = [=]( ST x ) ->ST { return( 1./lx ); };

	x.getCField()(Pimpact::F::U).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST { return( initFunC(x,y) ); } );
	x.getCField()(Pimpact::F::V).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST { return( initFunC(x,y) ); } );
	x.getSField()(Pimpact::F::U).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST { return( initFunS(x,y) ); } );
	x.getSField()(Pimpact::F::V).initFromFunction(
			[=]( ST x, ST y, ST z ) ->ST { return( initFunS(x,y) ); } );

	sol = x;
	if( write ) x.write( 10 );

	// solution init
	rhs.getCField()(Pimpact::F::U).initFromFunction(
	    	[=]( ST x, ST y, ST z ) ->ST {
					if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
						(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
						(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
						(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
						(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
						(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
						return( initFunC( std::min(std::max(x,0.),1.),std::min(std::max(y,0.),1.)) );
					else
						return( a2*initFunS(x,y) + deriFunC(y) ); } );

	rhs.getCField()(Pimpact::F::V).initFromFunction(
	    	[=]( ST x, ST y, ST z ) ->ST {
					if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
						(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
						(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
						(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
						(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
						(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
						return( initFunC( std::min(std::max(x,0.),1.),std::min(std::max(y,0.),1.)) );
					else
						return( a2*initFunS(x,y) + deriFunC(y) ); } );

	rhs.getSField()(Pimpact::F::U).initFromFunction(
				[=]( ST x, ST y, ST z ) ->ST {
					if( ((x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
						(  (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
						(  (y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
						(  (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
						(  (z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
						(  (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
						return( initFunS( std::min(std::max(x,0.),1.),std::min(std::max(y,0.),1.)) );
					else
						return( -a2*initFunC(x,y) +deriFunS(x) ); } );

	rhs.getSField()(Pimpact::F::V).initFromFunction(
				[=]( ST x, ST y, ST z ) ->ST {
					if( (( x   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0 ) ||
						(  ( x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0 ) ||
						(  ( y   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0 ) ||
						(  ( y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0 ) ||
						(  ( z   )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0 ) ||
						(  ( z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0 ) )
						return( initFunS( std::min(std::max(x,0.),1.),std::min(std::max(y,0.),1.)) );
					else
						return( -a2*initFunC(x,y) + deriFunS(x) ); } );

	//if( write ) rhs.write( 30 );

	op->apply( sol, y );
	//if( write ) y.write( 20 );
	if( 2==print ) y.print();
	if( 3==print ) rhs.print();

	err.add( 1., y, -1., rhs );
	//if( write ) err.write( 0 );
	if( 1==print ) err.print(   );

	ST error = err.norm(Belos::InfNorm)/rhs.norm(Belos::InfNorm);
	std::cout << "\nresidual: " << error << "\n";
	if( 1==domain )
		TEST_EQUALITY( error<(1./nx/ny), true );

	x.init();
	modeInv->apply( rhs, x );

	//if( write ) x.write( 20 );

	err.add( 1., sol, -1., x );
	if( write ) err.write( 0 );
	if( print ) err.print(   );

	error = err.norm(Belos::InfNorm);
	std::cout << "\nerror: " << error << "\n";
	if( 1==domain )
		TEST_EQUALITY( error<1.e-6, true );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ModeSolver, ModeNonlinearOp, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ModeSolver, ModeNonlinearOp, D3 )



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




} // namespace
