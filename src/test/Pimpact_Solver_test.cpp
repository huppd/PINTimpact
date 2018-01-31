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

using SpaceT = Pimpact::Space<ST, OT, sd, d, dNC>;


using SF = typename Pimpact::ScalarField<SpaceT>;
using VF = typename Pimpact::VectorField<SpaceT>;
using MSF = typename Pimpact::ModeField<SF>;
using MVF = typename Pimpact::ModeField<VF>;

template<class T1, class T2>
using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1, T2> >;
template<class T>
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T>
using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

template<class T>
using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;
template<class T>
using ConvDiffJT = Pimpact::NonlinearSmoother<T, Pimpact::ConvectionDiffusionJSmoother >;
template<class T>
using ConvDiffSORT = Pimpact::NonlinearSmoother<T, Pimpact::ConvectionDiffusionSORSmoother >;

template<class T> using ConvOpT = Pimpact::NonlinearOp<Pimpact::ConvectionSOp<T> >;

using ConvDiffOp2D = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<D2> >;
using ConvDiffOp3D = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<D3> >;

using ConvDiffJ2D = Pimpact::NonlinearSmoother<ConvDiffOp2D, Pimpact::ConvectionDiffusionJSmoother >;
using ConvDiffSOR2D = Pimpact::NonlinearSmoother<ConvDiffOp2D, Pimpact::ConvectionDiffusionSORSmoother >;

using ConvDiffJ3D = Pimpact::NonlinearSmoother<ConvDiffOp3D, Pimpact::ConvectionDiffusionJSmoother >;
using ConvDiffSOR3D = Pimpact::NonlinearSmoother<ConvDiffOp3D, Pimpact::ConvectionDiffusionSORSmoother >;

template<class T> using POP2  = Pimpact::PrecInverseOp<T, ConvDiffJT >;
template<class T> using MOP = Pimpact::InverseOp<T >;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(ModeSolver, ModeNonlinearOp, SpaceT) {

	pl->set<bool>("spectral in time", true);

	ST pi2 = 2.*std::acos(-1.);

	lx = pi2 ;
	ly = pi2 ;
	lz = pi2 ;

	nf = 4;

	setParameter(SpaceT::sdim);
	Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(pl);


	Pimpact::ModeField<Pimpact::VectorField<SpaceT> > x(space);
	Pimpact::ModeField<Pimpact::VectorField<SpaceT> > y(space);
	Pimpact::ModeField<Pimpact::VectorField<SpaceT> > rhs(space);
	Pimpact::ModeField<Pimpact::VectorField<SpaceT> > sol(space);
	Pimpact::ModeField<Pimpact::VectorField<SpaceT> > err(space);

	Teuchos::RCP<Teuchos::ParameterList > param = Teuchos::parameterList();// = Pimpact::createLinSolverParameter(solvName, 1.e-6);

	//std::string solvName = "TFQMR";
	std::string solvName = "GMRES";
	param->get<std::string>("Solver name", solvName);

	param->sublist("Solver").set("Output Style", Belos::Brief);
	param->sublist("Solver").set("Output Frequency", -1);
	param->sublist("Solver").set("Verbosity",			        
			Belos::Errors +
			Belos::Warnings +
			Belos::IterationDetails +
			Belos::OrthoDetails +
			Belos::FinalSummary +
			//Belos::TimingDetails +
			Belos::StatusTestDetails +
			Belos::Debug);
	param->sublist("Solver").set("Maximum Iterations", nMax);
	param->sublist("Solver").set("Convergence Tolerance", 1.e-6);
	//param->set("Flexible Gmres", false);
	
	auto zeroOp = Pimpact::create<ConvDiffOpT>(space);
	
	//Teuchos::RCP<const Pimpact::ModeNonlinearOp<ConvDiffOpT<SpaceT> > modeOp =
	auto modeOp =
		Teuchos::rcp(new Pimpact::ModeNonlinearOp<ConvDiffOpT<SpaceT> >(zeroOp));

	auto modeInv = Pimpact::createInverseOp(modeOp, param);


	using FSpaceT = SpaceT;
	using CSpaceT = Pimpact::Space<ST, OT, SpaceT::sdim, SpaceT::dimension, 2>;

	using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT, CSpaceT>;
	auto mgSpaces = Pimpact::createMGSpaces<CS>(space, maxGrids);

	auto zeroInv = Pimpact::createInverseOp(zeroOp, param);

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
		//POP2
		//POP3
		MOP
			> (mgSpaces, zeroOp, Teuchos::parameterList()) ;

	//if(0==space->rankST())
	//mgConvDiff->print();

	//zeroInv->setRightPrec(Pimpact::createMultiOperatorBase(mgConvDiff));
	//modeInv->setRightPrec(Pimpact::createMultiOperatorBase(Pimpact::create<Pimpact::ModePrec>(zeroInv)));
	//modeInv->setLeftPrec(Pimpact::createMultiOperatorBase(Pimpact::create<Pimpact::ModePrec>(zeroInv)));

	ST iRe = 1./space->getDomainSize()->getRe();
	ST a2 = space->getDomainSize()->getAlpha2()*iRe;

	//std::cout <<"a2: " <<a2 <<"\n";
	//std::cout <<"iRe: " <<iRe <<"\n";
	// computing zero mode of z
	// set paramteters
	auto para = Teuchos::parameterList();
	para->set<ST>("mulI", a2 );
	para->set<ST>("mulC", 1. );
	para->set<ST>("mulL", iRe);
	modeOp->setParameter(para);
	modeInv->setParameter(para);

	// initializtion
	{
		Pimpact::VectorField<SpaceT> wind(space);
		wind(Pimpact::F::U).init(1.);
		wind(Pimpact::F::V).init(1.);
		////wind(Pimpact::F::W).init(1.);

		zeroOp->assignField(wind);
		zeroInv->assignField(wind);
	}

	auto initFunC = [](ST x, ST y) ->ST { return(std::pow((y-0.5), 2)); };
	auto initFunS = [](ST x, ST y) ->ST { return(std::pow((x-0.5), 1)); };
	auto deriFunC = [=](ST y) ->ST { return(2.*(y-0.5)/ly - iRe*2./ly/ly); };
	auto deriFunS = [=](ST x) ->ST { return(1./lx); };

	x.getCField()(Pimpact::F::U).initFromFunction(
			[=](ST x, ST y, ST z) ->ST { return(initFunC(x, y)); });
	x.getCField()(Pimpact::F::V).initFromFunction(
			[=](ST x, ST y, ST z) ->ST { return(initFunC(x, y)); });
	x.getSField()(Pimpact::F::U).initFromFunction(
			[=](ST x, ST y, ST z) ->ST { return(initFunS(x, y)); });
	x.getSField()(Pimpact::F::V).initFromFunction(
			[=](ST x, ST y, ST z) ->ST { return(initFunS(x, y)); });

	sol = x;
	if(write) x.write(10);

	// solution init
	rhs.getCField()(Pimpact::F::U).initFromFunction(
	    	[=](ST x, ST y, ST z) ->ST {
					if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
						( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
						( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
						( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
						( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
						( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
						return(initFunC(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
					else
						return(a2*initFunS(x, y) + deriFunC(y)); });

	rhs.getCField()(Pimpact::F::V).initFromFunction(
	    	[=](ST x, ST y, ST z) ->ST {
					if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
						( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
						( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
						( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
						( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
						( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
						return(initFunC(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
					else
						return(a2*initFunS(x, y) + deriFunC(y)); });

	rhs.getSField()(Pimpact::F::U).initFromFunction(
				[=](ST x, ST y, ST z) ->ST {
					if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
						( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
						( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
						( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
						( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
						( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
						return(initFunS(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
					else
						return(-a2*initFunC(x, y) +deriFunS(x)); });

	rhs.getSField()(Pimpact::F::V).initFromFunction(
				[=](ST x, ST y, ST z) ->ST {
					if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
						( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
						( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
						( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
						( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
						( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
						return(initFunS(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
					else
						return(-a2*initFunC(x, y) + deriFunS(x)); });

	//if(write) rhs.write(30);

	modeOp->apply(sol, y);
	//if(write) y.write(20);
	if(2==print) y.print();
	if(3==print) rhs.print();

	err.add(1., y, -1., rhs);
	//if(write) err.write(0);
	if(1==print) err.print( );

	ST error = err.norm(Pimpact::ENorm::Inf)/rhs.norm(Pimpact::ENorm::Inf);
	std::cout <<"\nresidual: " <<error <<"\n";
	if(1==domain)
		TEST_EQUALITY(error<(1./nx/ny), true);

	x.init();
	modeInv->apply(rhs, x);

	//if(write) x.write(20);

	err.add(1., sol, -1., x);
	if(write) err.write(0);
	if(print) err.print( );

	error = err.norm(Pimpact::ENorm::Inf);
	std::cout <<"\nerror: " <<error <<"\n";
	if(1==domain)
		TEST_EQUALITY(error<1.e-6, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ModeSolver, ModeNonlinearOp, D2)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ModeSolver, ModeNonlinearOp, D3)




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiHarmonicSolver, MultiHarmonicDiagOp, SpaceT) {

	pl->set<bool>("spectral in time", true);

	ST pi2 = 2.*std::acos(-1.);

	lx = pi2 ;
	ly = pi2 ;
	lz = pi2 ;

	nf = 4;

	setParameter(SpaceT::sdim);
	Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(pl);


	Pimpact::MultiHarmonicField<Pimpact::VectorField<SpaceT> > x(space);
	Pimpact::MultiHarmonicField<Pimpact::VectorField<SpaceT> > y(space);
	Pimpact::MultiHarmonicField<Pimpact::VectorField<SpaceT> > rhs(space);
	Pimpact::MultiHarmonicField<Pimpact::VectorField<SpaceT> > sol(space);
	Pimpact::MultiHarmonicField<Pimpact::VectorField<SpaceT> > err(space);

	Teuchos::RCP<Teuchos::ParameterList > param = Teuchos::parameterList();// = Pimpact::createLinSolverParameter(solvName, 1.e-6);

	std::string solvName = "TFQMR";
	//std::string solvName = "GMRES";
	param->get<std::string>("Solver name", "TFQMR");

	param->sublist("Solver").set("Maximum Iterations", nMax);
	param->sublist("Solver").set("Convergence Tolerance", 1.e-6);
	//param->set("Flexible Gmres", false);
	
	auto zeroOp = Pimpact::create<ConvDiffOpT>(space);
	
	auto zeroInv = Pimpact::createInverseOp(zeroOp, param);

	auto modeOp =
		Teuchos::rcp(new Pimpact::ModeNonlinearOp<ConvDiffOpT<SpaceT> >(zeroOp));

	param->get<std::string>("Solver name", "GMRES");
	param->sublist("Solver").set("Output Style", Belos::Brief);
	param->sublist("Solver").set("Output Frequency", -1);
	param->sublist("Solver").set("Verbosity",			        
			Belos::Errors +
			Belos::Warnings +
			Belos::IterationDetails +
			Belos::OrthoDetails +
			Belos::FinalSummary +
			//Belos::TimingDetails +
			Belos::StatusTestDetails +
			Belos::Debug);

	auto modeInv = Pimpact::createInverseOp(modeOp, param);


	using FSpaceT = SpaceT;
	using CSpaceT = Pimpact::Space<ST, OT, SpaceT::sdim, SpaceT::dimension, 2>;

	using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT, CSpaceT>;
	auto mgSpaces = Pimpact::createMGSpaces<CS>(space, maxGrids);

  auto mgConvDiff = Pimpact::createMultiGrid<
    Pimpact::VectorField,
    TransVF,
    RestrVF,
    InterVF,
    ConvDiffOpT,
    ConvDiffOpT,
    //ConvDiffSORT,
    ConvDiffJT,
    //ConvDiffJT
    //POP2
    //POP3
    MOP
      > (mgSpaces, zeroOp, Teuchos::parameterList()) ;

	//if(0==space->rankST())
	//mgConvDiff->print();

	zeroInv->setRightPrec(Pimpact::createMultiOperatorBase(mgConvDiff));
	modeInv->setRightPrec(Pimpact::createMultiOperatorBase(Pimpact::create<Pimpact::ModePrec>(zeroInv)));
	//modeInv->setLeftPrec(Pimpact::createMultiOperatorBase(Pimpact::create<Pimpact::ModePrec>(zeroInv)));

	auto opV2Vprec = Pimpact::createMultiHarmonicDiagOp(zeroInv, modeInv);


	//std::cout <<"a2: " <<a2 <<"\n";
	//std::cout <<"iRe: " <<iRe <<"\n";
	// computing zero mode of z
	// set paramteters

	// initializtion
	{
		Pimpact::MultiHarmonicField<Pimpact::VectorField<SpaceT> > wind(space);
		wind.get0Field()(Pimpact::F::U).init(1.);
		wind.get0Field()(Pimpact::F::V).init(1.);
		////wind(Pimpact::F::W).init(1.);

		opV2Vprec->assignField(wind);
	}

	if(0==space()->si(Pimpact::F::U, 3)) {
		ST iRe = 1./space->getDomainSize()->getRe();

		auto initFunC = [](ST x, ST y) ->ST { return(std::pow((y-0.5), 2)); };
		auto initFunS = [](ST x, ST y) ->ST { return(std::pow((x-0.5), 1)); };
		auto deriFunC = [=](ST y) ->ST { return(2.*(y-0.5)/ly - iRe*2./ly/ly); };
		auto deriFunS = [=](ST x) ->ST { return(1./lx); };

		x.get0Field()(Pimpact::F::U).initFromFunction(
				[=](ST x, ST y, ST z) ->ST { return(initFunC(x, y)); });
		x.get0Field()(Pimpact::F::V).initFromFunction(
				[=](ST x, ST y, ST z) ->ST { return(initFunS(x, y)); });
		rhs.get0Field()(Pimpact::F::U).initFromFunction(
				[=](ST x, ST y, ST z) ->ST {
				if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
					( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
					( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
					( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
					( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
					( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
				return(initFunC(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
				else
				return(deriFunC(y)); });

		rhs.get0Field()(Pimpact::F::V).initFromFunction(
				[=](ST x, ST y, ST z) ->ST {
				if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
					( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
					( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
					( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
					( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
					( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
				return(initFunS(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
				else
				return(deriFunS(y)); });
	}

	for(OT i=std::max(space()->si(Pimpact::F::U, 3), 1); i<=space()->ei(Pimpact::F::U, 3); ++i) {
		ST iRe = 1./space->getDomainSize()->getRe();
		ST a2 = space->getDomainSize()->getAlpha2()*iRe*i;

		auto initFunC = [](ST x, ST y) ->ST { return(std::pow((y-0.5), 2)); };
		auto initFunS = [](ST x, ST y) ->ST { return(std::pow((x-0.5), 1)); };
		auto deriFunC = [=](ST y) ->ST { return(2.*(y-0.5)/ly - iRe*2./ly/ly); };
		auto deriFunS = [=](ST x) ->ST { return(1./lx); };

		x.getCField(i)(Pimpact::F::U).initFromFunction(
				[=](ST x, ST y, ST z) ->ST { return(initFunC(x, y)); });
		x.getCField(i)(Pimpact::F::V).initFromFunction(
				[=](ST x, ST y, ST z) ->ST { return(initFunC(x, y)); });
		x.getSField(i)(Pimpact::F::U).initFromFunction(
				[=](ST x, ST y, ST z) ->ST { return(initFunS(x, y)); });
		x.getSField(i)(Pimpact::F::V).initFromFunction(
				[=](ST x, ST y, ST z) ->ST { return(initFunS(x, y)); });

		rhs.getCField(i)(Pimpact::F::U).initFromFunction(
				[=](ST x, ST y, ST z) ->ST {
				if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
					( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
					( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
					( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
					( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
					( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
				return(initFunC(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
				else
				return(a2*initFunS(x, y) + deriFunC(y)); });

		rhs.getCField(i)(Pimpact::F::V).initFromFunction(
				[=](ST x, ST y, ST z) ->ST {
				if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
					( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
					( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
					( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
					( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
					( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
				return(initFunC(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
				else
				return(a2*initFunS(x, y) + deriFunC(y)); });

		rhs.getSField(i)(Pimpact::F::U).initFromFunction(
				[=](ST x, ST y, ST z) ->ST {
				if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
					( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
					( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
					( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
					( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
					( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
				return(initFunS(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
				else
				return(-a2*initFunC(x, y) +deriFunS(x)); });

		rhs.getSField(i)(Pimpact::F::V).initFromFunction(
				[=](ST x, ST y, ST z) ->ST {
				if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(0)>0) ||
					( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(0)>0) ||
					( (y  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(1)>0) ||
					( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(1)>0) ||
					( (z  )<= Teuchos::ScalarTraits<ST>::eps() && space->bcl(2)>0) ||
					( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && space->bcu(2)>0))
				return(initFunS(std::min(std::max(x, 0.), 1.), std::min(std::max(y, 0.), 1.)));
				else
				return(-a2*initFunC(x, y) + deriFunS(x)); });
	}

	sol = x;
	if(write) x.write(10);

	// solution init

	//if(write) rhs.write(30);

	//op->apply(sol, y);
	////if(write) y.write(20);
	//if(2==print) y.print();
	//if(3==print) rhs.print();

	//err.add(1., y, -1., rhs);
	////if(write) err.write(0);
	//if(1==print) err.print( );

	//ST error = err.norm(Pimpact::ENorm::Inf)/rhs.norm(Pimpact::ENorm::Inf);
	//std::cout <<"\nresidual: " <<error <<"\n";
	//if(1==domain)
		//TEST_EQUALITY(error<(1./nx/ny), true);

	x.init();
	opV2Vprec->apply(rhs, x);

	//if(write) x.write(20);

	err.add(1., sol, -1., x);
	if(write) err.write(0);
	if(print) err.print( );

	ST error = err.norm(Pimpact::ENorm::Inf);
	std::cout <<"\nerror: " <<error <<"\n";
	if(1==domain)
		TEST_EQUALITY(error<1.e-6, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiHarmonicSolver, MultiHarmonicDiagOp, D2)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiHarmonicSolver, MultiHarmonicDiagOp, D3)

} // namespace
