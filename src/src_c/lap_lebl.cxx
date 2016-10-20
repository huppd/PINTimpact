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
#include "Pimpact_FieldFactory.hpp"
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
#include "Pimpact_Types.hpp"
#include "Pimpact_Space.hpp"




using S = double;
using O = int;
const int dNC = 4;
//const int dNC = 2;

using SpaceT = Pimpact::Space<S,O,4,dNC>;

using FSpaceT = SpaceT;
using CSpaceT = Pimpact::Space<S,O,4,2>;

using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;

using VF = Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> >;
using SF = Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> >;

using MVF = Pimpact::MultiField<VF>;
using MSF = Pimpact::MultiField<SF>;

using CF = Pimpact::CompoundField< VF, SF>;

template<class T>
using POP = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;




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
	auto space = Pimpact::createSpace<S,O,4,dNC>( Teuchos::rcpFromRef(
				pl->sublist( "Space", true ) ) );

	// ------------ init vectors
	auto x = Pimpact::create<CF>( space );
	auto f = x->clone( Pimpact::ShallowCopy );

	// ------------ init Fields
	x->getVFieldPtr()->initField( pl->sublist("Base flow") );
	//x->getVFieldPtr()->write();


	// ------------ construct Operators
	auto opV2V = Pimpact::createMultiDtConvectionDiffusionOp( space );
	auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
	auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );
	opV2V->assignField( x->getVField() );

	// ------------ check velocity field for divergence freeness
	opV2S->apply( x->getVField(), f->getSField() );
	if( 0==space->rankST() )
		std::cout << "\ndiv(U) = " << f->getSFieldPtr()->norm() << "\n";

	// ------------ generate RHS
	opV2V->apply( x->getVField(), f->getVField() );
	f->getVFieldPtr()->extrapolateBC();
	opV2S->apply( f->getVField(), f->getSField() );

	//f->getSFieldPtr()->write( 99 );

	// ------------ construct inverse DivGrad
	auto divGradOp = Pimpact::createDivGradOp(
			opV2S->getOperatorPtr(), opS2V->getOperatorPtr() );

	auto divGradInv2 = Pimpact::createInverseOp( 
				divGradOp,
				Teuchos::rcpFromRef( pl->sublist("DivGrad") ) );

	// ------------ construct inverse DivGrad^T
	auto divGradInvTT =
		Pimpact::createInverseOp( 
				Pimpact::createTranspose( divGradOp ),
				Teuchos::rcpFromRef( pl->sublist("DivGrad^T") ) );

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

		// create Multi space
		auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>( space,
				pl->sublist("Multi Grid").get<int>("maxGrids") );

		auto mgDivGrad = Pimpact::createMultiGrid<
			Pimpact::ScalarField,
			Pimpact::TransferOp,
			Pimpact::RestrictionHWOp,
			Pimpact::InterpolationOp,
			Pimpact::DivGradOp,
			//Pimpact::DivGradOp,
			Pimpact::DivGradO2Op,
			Pimpact::DivGradO2JSmoother,
			//Pimpact::Chebyshev,
			//Pimpact::DivGradO2SORSmoother,
			//POP
			//Pimpact::Chebyshev
			//Pimpact::DivGradO2Inv
			//Pimpact::DivGradO2SORSmoother
			Pimpact::DivGradO2JSmoother
				>( mgSpaces, Teuchos::rcpFromRef( pl->sublist("DivGrad").sublist("Multi Grid") ) );

		if( 0==space->rankST() )
			mgDivGrad->print();

		if( "right" == divGradPrecString ) 
			divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
		if( "left" == divGradPrecString ) 
			divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
	}


	auto divGradInv  = Pimpact::createMultiHarmonicMultiOpWrap( divGradInv2 );
	auto divGradInvT = Pimpact::createMultiHarmonicMultiOpWrap( divGradInvTT );

	// solve nullspace IMPACT style
	auto nullspace = f->getSFieldPtr()->clone();
	{
		nullspace->init( 1. );
		auto zeros = nullspace->clone();
		divGradInvT->getOperatorPtr()->getOperatorPtr()->apply(
				*Pimpact::createMultiField( nullspace->get0FieldPtr() ),
				*Pimpact::createMultiField( zeros->get0FieldPtr() ));

		zeros->write( 666 );
		zeros->init( 0. );
		//zeros->random();
		divGradInvT->apply( *zeros, *nullspace );
		nullspace->write(888);
	}

	////// --------------------- solve with plaine rhs
	divGradInv->apply( f->getSField(), x->getSField() );
	//x->getSField().init(0.);
	//divGradInv->apply( *nullspace, x->getSField() );
	//x->getSField().write(999);


	// compute Residual I
	{
		auto res =  Pimpact::createMultiField( f->getSFieldPtr()->get0FieldPtr()->clone() );
		divGradInv->getOperatorPtr()->getLinearProblem()->getProblem()->getOperator()->apply(
				*Pimpact::createMultiField( x->getSFieldPtr()->get0FieldPtr() ),
				*res );
		res->getFieldPtr( 0 )->add( 1., f->getSFieldPtr()->get0Field(), -1., res->getField( 0 ) );
		std::cout << "||res||: " << res->norm() << "\n";
		res->write( 777 );
	}

	x->level();
	x->getSFieldPtr()->write();

	//// --------------------- solve with nullspaced rhs
	//auto xl = x->getSFieldPtr()->clone();
	//xl->init(0.);
	//auto f_ = f->getSFieldPtr()->clone( Pimpact::DeepCopy );
	//f_->get0FieldPtr()->add(
			//1., f->getSFieldPtr()->get0Field(),
			//-nullspace->get0FieldPtr()->dot( f->getSFieldPtr()->get0Field() )/nullspace->get0FieldPtr()->dot( nullspace->get0Field() ), nullspace->get0Field() );
	//std::cout << "(rhs,nullspace): " << f_->get0FieldPtr()->dot( nullspace->get0Field() ) << "\n";
	////divGradInv->apply( *f_, *xl );
	//divGradInv->apply( *nullspace, *xl );

	//xl->level();
	//xl->write( 10 );
	//// diff between solution
	//x->getSFieldPtr()->add( 1., x->getSField(), -1., *xl );
	//x->getSFieldPtr()->write(666);
	//std::cout << "||x-xl||: " << x->getSFieldPtr()->norm() << "\n";

	// compute resudial II
	//opS2V->apply( xl->getSField(), x->getVField() );
	//x->getVField().add( 1., xl->getVField(), 1., f->getVField() );
	//std::cout << "residual: " << x->getVFieldPtr()->norm() << "\n";

	Teuchos::TimeMonitor::summarize();

	Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

	MPI_Finalize();
	return( 0 );

}
