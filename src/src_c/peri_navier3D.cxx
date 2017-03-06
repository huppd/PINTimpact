#include <fstream>
#include <ostream>

#include <mpi.h>

#include <Teuchos_Array.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Tuple.hpp>
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTypes.hpp"

#include "NOX.H"

#include "BelosPimpactAdapter.hpp"
#include "NOX_Pimpact_Group.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_StatusTest.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_ModeNonlinearOp.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_MultiOpSmoother.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_RefinementStrategy.hpp"
#include "Pimpact_TransferCompoundOp.hpp"
#include "Pimpact_TransferMultiHarmonicOp.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_Space.hpp"

#include "Pimpact_DivGradNullSpace.hpp"




using S = double;
using O = int;

const int sd = 3;
//const int dNC = 4;
const int dNC = 3;

using SpaceT = Pimpact::Space<S,O,sd,4,dNC>;

using FSpaceT = Pimpact::Space<S,O,sd,4,dNC>;
using CSpaceT = Pimpact::Space<S,O,sd,4,2>;

//using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;
using CS = Pimpact::CoarsenStrategy<FSpaceT,CSpaceT>; // dirty fix: till gather isn't fixed

using VF = Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> >;
using SF = Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> >;

using MVF = Pimpact::MultiField<VF>;
using MSF = Pimpact::MultiField<SF>;

using CF = Pimpact::CompoundField< VF, SF>;
using MF = Pimpact::MultiField<CF>;

using BOp = Pimpact::OperatorBase<MF>;

using NV = NOX::Pimpact::Vector<MF>;


template<class T1,class T2>
using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T>
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T>
using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;


template<class T>
using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;

//template<class T> using PrecS = Pimpact::MultiOpSmoother< Pimpact::DivGradO2JSmoother<T> >;
//template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2SORSmoother >;
template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;
//template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::Chebyshev >;
template<class T> using POP2  = Pimpact::PrecInverseOp< T, ConvDiffJT >;
template<class T> using POP3  = Pimpact::PrecInverseOp< T, ConvDiffSORT >;







int main( int argi, char** argv ) {


	/////////////////////////////////////////// set up parameters ///////////////////////////
	// intialize MPI
	MPI_Init( &argi, &argv );

	Teuchos::CommandLineProcessor my_CLP;

	std::string xmlFilename = "parameter3D.xml";
	my_CLP.setOption("filename", &xmlFilename, "file name of the input xml paramerterlist");

	my_CLP.recogniseAllOptions(true);
	my_CLP.throwExceptions(true);

	my_CLP.parse(argi,argv);


	auto pl = Teuchos::getParametersFromXmlFile( xmlFilename );
	pl->print();

	/////////////////////////////////////////// end of set up parameters /////////////////////////


	///////////////////////////////////////////  set up initial stuff ////////////////////////////

	std::string initGuess = pl->sublist("Solver").get<std::string>( "initial guess", "zero" );

	int withoutput=pl->sublist("Solver").get<int>( "withoutput", 1 );

	int refinement     = pl->sublist("Solver").get<int>( "refinement level", 1     );
	S   refinementTol  = pl->sublist("Solver").get<S>(   "refinement tol",   1.e-6 );
	int refinementStep = pl->sublist("Solver").get<int>( "refinement step",  2     );

	Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(
			Teuchos::rcpFromRef( pl->sublist( "Space", true ) ) );


	// init vectors
	Teuchos::RCP<MF> x = wrapMultiField(
			createCompoundField(
				Teuchos::rcp( new VF(space,true) ),
				Teuchos::rcp( new SF(space) ) ) ) ;

	// init Fields
	x->getField(0).getVField().initField( pl->sublist("Base flow") );
	//if( withoutput ) x->write( 900 );



	/******************************************************************************************/
	for( int refine=0; refine<refinement; ++refine ) {

		auto opV2V = Pimpact::createMultiDtConvectionDiffusionOp( space );
		auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
		auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );

		{
			opV2S->apply( x->getField(0).getVField(), x->getField(0).getSField()  );
			S divergence = x->getField(0).getSField().norm();
			if( 0==space->rankST() )
				std::cout << "\ndiv(Base Flow): " << divergence << "\n\n";
			x->getField(0).getSField().init( 0. );
		}
		auto op = Pimpact::createCompoundOpWrap(
				opV2V,
				opS2V,
				opV2S );

		std::string rl = "";
		if( refinement>1 )
			rl = std::to_string( static_cast<long long>(refine) ); // long long needed on brutus(intel)

		auto fu = x->clone( Pimpact::ECopy::Shallow );

		{
			// to get the Dirichlet for the RHS (necessary interpolation) ugly
			// super ugly hack for BC::Dirichlet
			opV2V->apply( x->getField(0).getVField(), fu->getField(0).getVField() );
			fu->init( 0., Pimpact::B::N );
		}

		// Taylor Green Vortex
		std::string forceType = pl->sublist("Force").get<std::string>("force type","Dirichlet");
		if( "force"== forceType )
			fu->getField(0).getVField().initField( pl->sublist("Force"), Pimpact::Add::Y );
		else if( "Taylor-Green"==forceType ) {
			S pi2 = 2.*std::acos(-1.);
			S alpha2 = space->getDomainSize()->getAlpha2();
			S re = space->getDomainSize()->getRe();
			S A =  pl->sublist("Force").get<S>("A", 1.);
			S B =  pl->sublist("Force").get<S>("B",-1.);
			S C =  pl->sublist("Force").get<S>("C", 0.);
			S a =  pl->sublist("Force").get<S>("a", 1.);
			S b =  pl->sublist("Force").get<S>("b", 1.);
			S c =  pl->sublist("Force").get<S>("c", 1.);
			TEUCHOS_TEST_FOR_EXCEPT( std::abs( a*A + b*B + c*C )>1.e-16 );

			//fu->getField(0).getVField().get0Field()(Pimpact::F::U).initFromFunction(
					//[&pi2,&re]( S x, S y, S z ) ->S {  return(  -std::sin(2.*x*pi2)/4. ); } );
			//fu->getField(0).getVField().get0Field()(Pimpact::F::V).initFromFunction(
					//[&pi2,&re]( S x, S y, S z ) ->S {  return(  -std::sin(2.*y*pi2)/4. ); } );

			fu->getField(0).getVField().getCField(1)(Pimpact::F::U).initFromFunction(
					[&]( S x, S y, S z ) ->S { return( alpha2*A*std::cos(a*x*pi2)*std::sin(b*y*pi2)*std::sin(c*z*pi2)/re ); } );
			fu->getField(0).getVField().getCField(1)(Pimpact::F::V).initFromFunction(
					[&]( S x, S y, S z ) ->S { return( alpha2*B*std::sin(a*x*pi2)*std::cos(b*y*pi2)*std::sin(c*z*pi2)/re ); } );
			fu->getField(0).getVField().getCField(1)(Pimpact::F::W).initFromFunction(
					[&]( S x, S y, S z ) ->S { return( alpha2*C*std::sin(a*x*pi2)*std::sin(b*y*pi2)*std::cos(c*z*pi2)/re ); } );

			fu->getField(0).getVField().getSField(1)(Pimpact::F::U).initFromFunction(
					[&]( S x, S y, S z ) ->S { return( 3.*A*std::cos(a*x*pi2)*std::sin(b*y*pi2)*std::sin(c*z*pi2)/re ); } );
			fu->getField(0).getVField().getSField(1)(Pimpact::F::V).initFromFunction(
					[&]( S x, S y, S z ) ->S { return( 3.*B*std::sin(a*x*pi2)*std::cos(b*y*pi2)*std::sin(c*z*pi2)/re ); } );
			fu->getField(0).getVField().getSField(1)(Pimpact::F::W).initFromFunction(
					[&]( S x, S y, S z ) ->S { return( 3.*C*std::sin(a*x*pi2)*std::sin(b*y*pi2)*std::cos(c*z*pi2)/re ); } );

			//fu->getField(0).getVField().getCField(2)(Pimpact::F::U).initFromFunction(
					//[&pi2]( S x, S y, S z ) ->S { return( std::sin(2.*x*pi2)/4. ); } );
			//fu->getField(0).getVField().getCField(2)(Pimpact::F::V).initFromFunction(
					//[&pi2]( S x, S y, S z ) ->S { return( std::sin(2.*y*pi2)/4. ); } );
		}

		if( 0==refine ) {
			if( "zero"==initGuess )
				x->init( 0. );
			else if( "almost zero"==initGuess ) {
				x->random();
				x->scale(1.e-12);
			}
			else if( "random"==initGuess )
				x->random();
			else if( "exact"==initGuess ) {
				x->getField(0).getSField().get0Field().initField( Pimpact::Grad2D_inX, -2./space->getDomainSize()->getRe() );
			}
			x->getField(0).getVField().changed();
			x->getField(0).getSField().changed();
		}

		//if( withoutput ) fu->write( 90000 );

		pl->sublist("Picard Solver").sublist("Solver").set( "Output Stream", Pimpact::createOstream("Picard"+rl+".txt", space->rankST() ) );

		auto opInv = Pimpact::createInverseOp( op, Teuchos::rcpFromRef( pl->sublist("Picard Solver") ) );

		////--- nullspace 
		if( pl->sublist( "Picard Solver" ).get<bool>( "nullspace ortho", true ) ) {
			auto nullspace = x->clone( Pimpact::ECopy::Shallow );

			Pimpact::DivGradNullSpace<Pimpact::DivOp<SpaceT> > compNullspace;

			compNullspace.computeNullSpace( opV2S->getOperatorPtr(),
					nullspace->getField(0).getSField().get0Field(), true );

			nullspace->getField(0).getVField().get0Field()(Pimpact::F::U).initFromFunction(
					[&space]( S x, S y, S z ) -> S { return( ( (Pimpact::BC::Dirichlet==space->bcl(Pimpact::X)&&x<=0.)?1.:0.) + ( (Pimpact::BC::Dirichlet==space->bcu(Pimpact::X)&&1.<=x)?-1.:0.) ); } );
			nullspace->getField(0).getVField().get0Field()(Pimpact::F::V).initFromFunction(
					[&space]( S x, S y, S z ) -> S { return( ( (Pimpact::BC::Dirichlet==space->bcl(Pimpact::Y)&&y<=0.)?1.:0.) + ( (Pimpact::BC::Dirichlet==space->bcu(Pimpact::Y)&&1.<=y)?-1.:0.) ); } );
			nullspace->getField(0).getVField().get0Field()(Pimpact::F::W).initFromFunction(
					[&space]( S x, S y, S z ) -> S { return( ( (Pimpact::BC::Dirichlet==space->bcl(Pimpact::Z)&&z<=0.)?1.:0.) + ( (Pimpact::BC::Dirichlet==space->bcu(Pimpact::Z)&&1.<=z)?-1.:0.) ); } );

			S blup = std::sqrt( 1./nullspace->dot( *nullspace ) );
			nullspace->scale( blup );

			for( int i=1; i<=space->nGlo(3);++i) {
				nullspace->getField(0).getSField().getCField(i) =
					nullspace->getField(0).getSField().get0Field();
				nullspace->getField(0).getSField().getSField(i) =
					nullspace->getField(0).getSField().get0Field();

				nullspace->getField(0).getVField().getCField(i) =
					nullspace->getField(0).getVField().get0Field();
				nullspace->getField(0).getVField().getSField(i) =
					nullspace->getField(0).getVField().get0Field();
			}


			//nullspace->write(999);
			auto out = Pimpact::createOstream( "nullspace.txt", space->rankST() );
			nullspace->print( *out );

			opInv->setNullspace( nullspace );
		}
		//// --- end nullspace

		/*** init preconditioner ****************************************************************************/

		std::string picardPrecString = pl->sublist("Picard Solver").get<std::string>( "preconditioner", "none" );

		if( "none" != picardPrecString ) { 

			// create Multi space
			auto mgSpaces = Pimpact::createMGSpaces<CS>( space, pl->sublist("Multi Grid").get<int>("maxGrids") );

			/////////////////////////////////////////begin of opv2v//////////////////////////////////////

			pl->sublist("MH_ConvDiff").sublist("Solver").set(
					"Output Stream",
					Pimpact::createOstream( opV2V->getLabel()+rl+".txt", space->rankST() ) );
			auto opV2Vinv = Pimpact::createInverseOp( opV2V, Teuchos::rcpFromRef(
						pl->sublist("MH_ConvDiff") ) );

			std::string mhConvDiffPrecString = pl->sublist("MH_ConvDiff").get<std::string>( "preconditioner", "right" );

			if( "none" != mhConvDiffPrecString ) { 
				// creat H0-inv prec
				auto zeroOp = Pimpact::create<ConvDiffOpT>( space );

				pl->sublist("ConvDiff").sublist("Solver").set( "Output Stream",
						Pimpact::createOstream( zeroOp->getLabel()+rl+".txt", space->rankST() ) );

				auto zeroInv = Pimpact::createMultiOpUnWrap( Pimpact::createInverseOp(
							zeroOp, Teuchos::rcpFromRef( pl->sublist("ConvDiff") ) ) );

				auto modeOp = Teuchos::rcp( new
						Pimpact::ModeNonlinearOp< ConvDiffOpT<SpaceT> >( zeroOp ) );

				pl->sublist("M_ConvDiff").sublist("Solver").set(
						"Output Stream",
						Pimpact::createOstream( modeOp->getLabel()+rl+".txt", space->rankST() ) );

				auto modeInv = Pimpact::createMultiOpUnWrap( Pimpact::createInverseOp(
							modeOp, Teuchos::rcpFromRef( pl->sublist("M_ConvDiff") ) ) );

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
					MOP
					//POP2
					//POP3
						> ( mgSpaces, Teuchos::rcpFromRef( pl->sublist("ConvDiff").sublist("Multi Grid") ) ) ;

				//if( 0==space->rankST() )
				//mgConvDiff->print();

				std::string convDiffPrecString = pl->sublist("ConvDiff").get<std::string>( "preconditioner", "right" );
				if( "right" == convDiffPrecString ) 
					zeroInv->getOperatorPtr()->setRightPrec( Pimpact::createMultiOperatorBase(mgConvDiff) );
				modeInv->getOperatorPtr()->setRightPrec( Pimpact::createMultiOperatorBase( Pimpact::create<Pimpact::EddyPrec>(zeroInv) ) );

				// create Hinv prec
				Teuchos::RCP<Pimpact::OperatorBase<MVF> > opV2Vprec = 
					Pimpact::createMultiOperatorBase(
							Pimpact::createMultiHarmonicDiagOp(
								zeroInv, modeInv ) );

				if( "right" == mhConvDiffPrecString ) 
					opV2Vinv->setRightPrec( opV2Vprec );
				if( "left" == mhConvDiffPrecString ) 
					opV2Vinv->setLeftPrec( opV2Vprec );
			}


			/////////////////////////////////////////end of opv2v//////////////////////////////////////
			////--- inverse DivGrad

			pl->sublist("DivGrad").sublist("Solver").set( "Output Stream",
					Pimpact::createOstream( "DivGrad"+rl+".txt", space->rankST() ) );

			auto divGradOp =
				Pimpact::createDivGradOp(
						opV2S->getOperatorPtr(),
						opS2V->getOperatorPtr() );

			auto divGradInv2 =
				Pimpact::createInverseOp( 
						divGradOp,
						Teuchos::rcpFromRef( pl->sublist("DivGrad") ) );

			////--- nullspace 
			if( pl->sublist("DivGrad").get<bool>("nullspace ortho",true) ) {
				auto nullspace = Pimpact::wrapMultiField( x->getField(0).getSField().get0Field().clone() );
				auto zeros = nullspace->clone();

				Pimpact::DivGradNullSpace<Pimpact::DivOp<SpaceT> > compNullspace;// = Pimpact::create<Pimpact::DivGradNullSpace>( );

				compNullspace.computeNullSpace( divGradOp->getDivOp(), nullspace->getField(0) );
				//nullspace->init( 1. );
				//zeros->init( 0. );

				//auto divGradInvTT =
					//Pimpact::createInverseOp( 
							//Pimpact::createTranspose( divGradOp ),
							//Teuchos::rcpFromRef( pl->sublist("DivGrad^T") ) );

				//std::string divGradScalString =
					//pl->sublist("DivGrad^T").get<std::string>("scaling","none");
				//if( "none" != divGradScalString ) { 
					//if( "left" != divGradScalString )
						//divGradInvTT->setLeftPrec( Pimpact::createMultiOperatorBase(
									//Pimpact::createInvDiagonal( divGradOp ) ) );
					//if( "right" != divGradScalString )
						//divGradInvTT->setRightPrec( Pimpact::createMultiOperatorBase(
									//Pimpact::createInvDiagonal( divGradOp ) ) );
				//}

				//divGradInvTT->apply( *zeros, *nullspace );

				//nullspace->write(888);
				auto out = Pimpact::createOstream( "nullspace.txt", space->rankST() );
				nullspace->print( *out );

				//S blup = std::sqrt( 1./nullspace->dot(*nullspace) );
				//nullspace->scale( blup );
				divGradInv2->setNullspace( nullspace );
			}
			//// --- end nullspace

			std::string divGradScalString =
				pl->sublist("DivGrad").get<std::string>("scaling","none");

			if( "none" != divGradScalString ) { 
				if( "left" != divGradScalString )
					divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase(
								Pimpact::createInvDiagonal( divGradOp ) ) );
				if( "right" != divGradScalString )
					divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase(
								Pimpact::createInvDiagonal( divGradOp ) ) );
			}

			std::string divGradPrecString =
				pl->sublist("DivGrad").get<std::string>("preconditioner","none");

			if( "none" != divGradPrecString ) { // init multigrid divgrad

				auto mgDivGrad = Pimpact::createMultiGrid<
					Pimpact::ScalarField,
					Pimpact::TransferOp,
					Pimpact::RestrictionSFOp,
					Pimpact::InterpolationOp,
					Pimpact::DivGradOp,
					Pimpact::DivGradO2Op,
					Pimpact::DivGradO2JSmoother,
					//Pimpact::Chebyshev,
					//Pimpact::DivGradO2SORSmoother,
					MOP
					//Pimpact::Chebyshev
					//Pimpact::DivGradO2Inv
					//Pimpact::DivGradO2SORSmoother
					//Pimpact::DivGradO2JSmoother
						>( mgSpaces, Teuchos::rcpFromRef( pl->sublist("DivGrad").sublist("Multi Grid") ) );

				if( 0==space->rankST() )
					mgDivGrad->print();

				if( "right" == divGradPrecString ) 
					divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
				if( "left" == divGradPrecString ) 
					divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
			}

			auto divGradInv =
				Pimpact::createMultiHarmonicMultiOpWrap(
						divGradInv2 );

			// schurcomplement approximator ...
			auto opSchur =
				Pimpact::createTripleCompositionOp(
						opV2S,
						opV2V,
						opS2V );

			auto opS2Sinv =
				Pimpact::createTripleCompositionOp(
						divGradInv,
						opSchur,
						divGradInv );

			if( space->rankST()==0 )
				std::cout << opS2Sinv->getLabel() << "\n";

			auto invTriangOp =
				Pimpact::createInverseTriangularOp(
						Pimpact::createMultiOpUnWrap( opV2Vinv ),
						opS2V,
						opS2Sinv );

			if( "right" == picardPrecString ) 
				opInv->setRightPrec( Pimpact::createMultiOperatorBase( invTriangOp ) );
			if( "left" == picardPrecString ) 
				opInv->setLeftPrec( Pimpact::createMultiOperatorBase( invTriangOp ) );
		}
		/*** end of init preconditioner ****************************************************************************/

		auto inter = NOX::Pimpact::createInterface( fu, createMultiOpWrap(op), opInv );

		auto nx = NOX::Pimpact::createVector( x );

		auto bla = Teuchos::parameterList();

		auto group = NOX::Pimpact::createGroup( bla, inter, nx );

		// Set up the status tests
		auto statusTest =
			NOX::StatusTest::buildStatusTests( pl->sublist("NOX Solver").sublist("Status Test"), NOX::Utils() );

		// Create the solver
		Teuchos::RCP<NOX::Solver::Generic> solver =
			NOX::Solver::buildSolver( group, statusTest, Teuchos::rcpFromRef( pl->sublist("NOX Solver") ) );


		Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );
		if( 0==space->rankST() )
			std::cout << "\n\t--- Nf: "<< space->nGlo(3) <<"\tdof: "<<x->getLength()<<"\t---\n";

		// Solve the nonlinear system
		{
			Teuchos::TimeMonitor LocalTimer(*Teuchos::TimeMonitor::getNewCounter("Pimpact:: Solving Time"));
			try{
				solver->solve();
			}
			catch( std::logic_error & e ) {
				std::cout << e.what() << "\n";
			}
		}


		// Get the answer
		*group = solver->getSolutionGroup();

		x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr();
		if( withoutput ) {
			Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getField().level();
			Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getField().write( refine*1000 );
			//Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getField().getField(0).getVField().write( 500+refine*1000, true );
			//Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() ))->getField().write( refine*1000 );
		}

		{
			auto out = Pimpact::createOstream( "energy_dis"+rl+".txt", space->rankST() );

			*out << 0 << "\t" << x->getField(0).getVField().get0Field().norm() << "\t" << std::sqrt( static_cast<S>( x->getField(0).getVField().get0Field().getLength() ) ) << "\n";
			for( int i=1; i<=space->nGlo(3); ++i )
				*out << i << "\t" << x->getField(0).getVField().getField(i).norm() << "\t" << std::sqrt( static_cast<S>( x->getField(0).getVField().getField(i).getLength() ) ) << "\n";
		}

		// spectral refinement criterion
		if( refinement>1 ) {

			x->getField(0).getVField().exchange();
			S u_nf = x->getField(0).getVField().getField(space->nGlo(3)).norm();
			S u_1  = x->getField(0).getVField().getField(1             ).norm();
			S truncError = 1.;
			if( u_nf != u_1 ) // just in case u_1=u_nf=0
				truncError = u_nf / u_1 ;
			if( truncError < refinementTol ) {
				if( 0==space->rankST() )
					std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " < " << refinementTol << "\n\n"; 
				break;
			}	
			else
				if( 0==space->rankST() )
					std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " >= " << refinementTol << "\n\n"; 

			auto spaceF =
				Pimpact::RefinementStrategy<SpaceT>::createRefinedSpace(
						space, Teuchos::tuple<int>( 0, 0, 0, refinementStep ) );

			auto refineOp =
				Teuchos::rcp(
						new Pimpact::TransferCompoundOp<
						Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp<SpaceT> > >,
						Pimpact::TransferMultiHarmonicOp< Pimpact::InterpolationOp<SpaceT> >
						>( space, spaceF ) );
			// refineOp->print();

			// init Fields for fine Boundary conditions
			auto xf = Pimpact::create<CF>( spaceF );
			xf->getVField().initField( pl->sublist("Base flow") );

			auto temp = Pimpact::create<CF>( spaceF );

			refineOp->apply( x->getField(0), *temp );

			xf->add( 1., *temp, 0., *temp );

			x = Pimpact::wrapMultiField( xf );
			space = spaceF;
		}
	} // end of for( int refine=0; refine<refinement; ++refine ) {
	/******************************************************************************************/

	Teuchos::TimeMonitor::summarize();

	if( 0==space->rankST() )
		Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

	MPI_Finalize();
	return( 0 );
}
