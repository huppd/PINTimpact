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

using SpaceT = Pimpact::Space<S,O,4,4>;

using FSpaceT = Pimpact::Space<S,O,4,4>;
using CSpaceT = Pimpact::Space<S,O,4,2>;

using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;

using VF = Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> >;
using SF = Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> >;

using MVF = Pimpact::MultiField<VF>;
using MSF = Pimpact::MultiField<SF>;

using CF = Pimpact::CompoundField< VF, SF>;


template<class T1,class T2>
using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T>
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionHWOp<T> >;
template<class T>
using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;


template<class T>
using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;
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
	auto space =
		Pimpact::createSpace<S,O,4,4>( Teuchos::rcpFromRef( pl->sublist( "Space", true ) ) );


	// init vectors
	auto x = Pimpact::create<CF>( space );

	// init Fields
	x->getVFieldPtr()->initField( pl->sublist("Base flow") );

	std::string rl = "";

	auto f = x->clone( Pimpact::ShallowCopy );


	auto opV2V = Pimpact::createMultiDtConvectionDiffusionOp( space );
	auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
	auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );
	opV2V->assignField( x->getVField() );

	// create Multi space
	auto mgSpaces =
		Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>(
				space,
				pl->sublist("Multi Grid").get<int>("maxGrids") );



	// --- inverse DivGrad
	//pl->sublist("DivGrad").sublist("Solver").set(
			//"Output Stream",
			//Pimpact::createOstream( "DivGrad"+rl+".txt", space->rankST() ) );

	auto divGradInv2 =
		Pimpact::createInverseOp( 
				Pimpact::createDivGradOp(
					opV2S->getOperatorPtr(),
					opS2V->getOperatorPtr() ),
					Teuchos::rcpFromRef( pl->sublist("DivGrad") ));


	auto mgDivGrad = 
		Pimpact::createMultiGrid<
		Pimpact::ScalarField,
		Pimpact::TransferOp,
		Pimpact::RestrictionHWOp,
		Pimpact::InterpolationOp,
		Pimpact::DivGradOp,
		Pimpact::DivGradO2Op,
		Pimpact::DivGradO2JSmoother,
		//Pimpact::DivGradO2SORSmoother,
		//POP
		Pimpact::DivGradO2Inv
			>( mgSpaces, Teuchos::rcpFromRef( pl->sublist("DivGrad").sublist("Multi Grid") ) );

	if( 0==space->rankST() )
		mgDivGrad->print();

	divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );

	auto divGradInv =
		Pimpact::createMultiHarmonicMultiOpWrap(
				divGradInv2
				);

	opV2V->apply( x->getVField(), f->getVField() );
	opV2S->apply( f->getVField(), f->getSField() );

	f->getSFieldPtr()->write( 99 );
	//x->getSFieldPtr()->random();
	x->getSFieldPtr()->get0FieldPtr()->initField( Pimpact::Grad2D_inX, -2./space->getDomainSize()->getRe() );
	//for( int i=0; i<10; ++i )
		//mgDivGrad->apply( f->getSFieldPtr()->get0Field(), x->getSFieldPtr()->get0Field() );
	divGradInv->apply( f->getSField(), x->getSField() );

	x->getSFieldPtr()->write();

	Teuchos::TimeMonitor::summarize();

	Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

	MPI_Finalize();
	return( 0 );

}
