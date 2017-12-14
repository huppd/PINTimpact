#include <fstream>
#include <ostream>

#include <mpi.h>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTypes.hpp"

#include "NOX.H"
#include "NOX_PrePostOperator_Vector.H"

#include "BelosPimpactAdapter.hpp"
#include "NOX_Pimpact_Group.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_PrePostError.hpp"
#include "NOX_Pimpact_PrePostEnergy.hpp"
#include "NOX_Pimpact_PrePostSpectrum.hpp"
#include "NOX_Pimpact_PrePostWriter.hpp"
#include "NOX_Pimpact_PrePostWriteRestart.hpp"

#include "Pimpact_AnalysisTools.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_DivGradProjector.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_ModeNonlinearOp.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_MultiOpSmoother.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_PicardProjector.hpp"
#include "Pimpact_RefinementStrategy.hpp"
#include "Pimpact_TransferCompoundOp.hpp"
#include "Pimpact_TransferMultiHarmonicOp.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_Space.hpp"

#include "Pimpact_DivGradNullSpace.hpp"




using ST = double;
using OT = int;

const int sd = 3;
const int dNC = 4;
//const int dNC = 3;
//const int dNC = 2;

using SpaceT = Pimpact::Space<ST,OT,sd,4,dNC>;


using VF = Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> >;
using SF = Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> >;

using MVF = Pimpact::MultiField<VF>;
using MSF = Pimpact::MultiField<SF>;

using CF = Pimpact::CompoundField< VF, SF>;
using MF = Pimpact::MultiField<CF>;



int main( int argi, char** argv ) {

  // intialize MPI
  MPI_Init( &argi, &argv );

  {
    /////////////////////////////////////////// set up parameters ///////////////////////////
    Teuchos::CommandLineProcessor my_CLP;

    std::string xmlFilename = "parameter3D.xml";
    my_CLP.setOption("filename", &xmlFilename, "file name of the input xml parameterlist");
    int restart = -1;
    my_CLP.setOption("restart", &restart, "number of restart");

    my_CLP.recogniseAllOptions(true);
    my_CLP.throwExceptions(true);

    my_CLP.parse(argi,argv);

    Teuchos::RCP<Teuchos::ParameterList> pl;
    pl = Teuchos::getParametersFromXmlFile( "parameterOut.xml" );
    ////////////////////////////////////////// end of set up parameters /////////////////////////


    //////////////////////////////////////////  set up initial stuff ////////////////////////////
    Teuchos::RCP<const SpaceT> space =
      Pimpact::create<SpaceT>( Teuchos::sublist( pl, "Space", true ) );


    if( 0==space->rankST() ) std::cout << "initial field\n";

    // init vectors
    Teuchos::RCP<MF> x =
      wrapMultiField( createCompoundField(
            Teuchos::rcp( new VF(space,true) ),
            Teuchos::rcp( new SF(space) ) ) ) ;

    // init Fields
    x->getField(0).getVField().initField( pl->sublist("Base flow") );

    auto base = x->getField(0).getVField().get0Field().clone(Pimpact::ECopy::Deep);

    x->getField(0).read( restart );
    /*********************************************************************************/
    Pimpact::ECoord dir = Pimpact::ECoord::Y;
    ST gamma = 10.;
    std::string prefix = "energy_";

    // compute glob energy in y-dir
    if( 0==space->si(::Pimpact::F::U,3) ) {
      auto vel =  x->getField(0).getVField().get0Field().clone( ::Pimpact::ECopy::Deep );
      vel->add( 1., *vel, -1., *base );

      auto out = Pimpact::createOstream( prefix + "0.txt",
          space->getProcGrid()->getRankBar(dir));

      Pimpact::computeHEEnergyDir( *vel, *out, gamma );
    }

    for( OT i=std::max(space->si(::Pimpact::F::U,3),1); i<=space->ei(::Pimpact::F::U,3); ++i ) {
      {
        auto out = Pimpact::createOstream( prefix + "C"+std::to_string(i) + ".txt",
            space->getProcGrid()->getRankBar(dir) );

        Pimpact::computeHEEnergyDir(
            x->getField(0).getVField().getCField(i), *out, gamma );
      }
      {
        auto out = Pimpact::createOstream( prefix + "S"+std::to_string(i) + ".txt",
            space->getProcGrid()->getRankBar(dir) );

        Pimpact::computeHEEnergyDir(
            x->getField(0).getVField().getSField(i), *out, gamma );
      }
    }

    x->write();


    //pl->sublist("NOX spectrum").set<int>( "refinement", refine );
    //prePostOperators->pushBack( 
    //Teuchos::rcp(new NOX::Pimpact::PrePostSpectrum<NV>(Teuchos::sublist(pl, "NOX spectrum"))));

    //noxSolverPara->sublist("Solver Options").set<Teuchos::RCP<NOX::Abstract::PrePostOperator>>(
    //"User Defined Pre/Post Operator", prePostOperators);

    //Teuchos::RCP<NOX::Solver::Generic> solver =
    //NOX::Solver::buildSolver( group, statusTest, noxSolverPara );


    //if( 0==space->rankST() )
    //std::cout << "\n\t--- Nf: "<< space->nGlo(3) <<"\tdof: "<<x->getLength()<<"\t---\n";

    //// write ParameterList for restart
    //if( 0==space->rankST() ) {
    //pl->sublist("Space").set<OT>( "nf", space->nGlo(3) );
    //pl->sublist("NOX Solver").sublist("Solver Options").remove("Status Test Check Type"); // dirty fix probably, should be fixed in NOX
    //pl->sublist("NOX Solver").sublist("Solver Options").remove("User Defined Merit Function"); // dirty fix probably, should be fixed in NOX, just needed for restart
    //Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );
    //}
    //// Solve the nonlinear system
    //NOX::StatusTest::StatusType status = NOX::StatusTest::StatusType::Failed;
    //{
    //Teuchos::TimeMonitor LocalTimer(*Teuchos::TimeMonitor::getNewCounter("Pimpact:: Solving Time"));
    //try {
    //status = solver->solve();
    //} catch( std::logic_error & e ) {
    //std::cout << e.what() << "\n";
    //}
    //}


    //// Get the answer
    //*group = solver->getSolutionGroup();

    //x = Teuchos::rcp_const_cast<NV>(
    //Teuchos::rcp_dynamic_cast<const NV>(group->getXPtr()))->getFieldPtr();


    //// spectral refinement
    //if( maxRefinement>1 ) {
    //ST u_1, u_nf;

    ////space->print();
    //if( 0<space->nGlo(3) and space()->si(Pimpact::F::U,3)<=1 and 1<=space()->ei(Pimpact::F::U,3) )
    //u_1  = x->getField(0).getVField().getField(1).norm( Pimpact::ENorm::L2 );
    //if( 0<space->nGlo(3) and space()->si(Pimpact::F::U,3)<=space->nGlo(3) and space->nGlo(3)<=space()->ei(Pimpact::F::U,3) )
    //u_nf = x->getField(0).getVField().getField(space->nGlo(3)).norm( Pimpact::ENorm::L2 );


    //int rank_1 = 0;
    //if( 1==space->getProcGrid()->getNP(3) )
    //rank_1 = 0;
    //else if( 0==(space->nGlo(3)+1)%space->getProcGrid()->getNP(3) )
    //rank_1 = 1;
    //int rank_nf = space->getProcGrid()->getNP(3)-1;

    //// nice nonblocking version
    //MPI_Request req_1, req_nf;  

    //MPI_Ibcast(
    //&u_1,                                // buffer	starting address of buffer (choice)
    //1,                                   // number of entries in buffer (non-negative integer)
    //MPI_DOUBLE,                          // data type of buffer (handle)
    //rank_1,                              // rank of broadcast root (integer)
    //space->getProcGrid()->getCommBar(3), // communicator (handle)
    //&req_1);                             // communication request
    //MPI_Ibcast(
    //&u_nf,                               // buffer	starting address of buffer (choice)
    //1,                                   // number of entries in buffer (non-negative integer)
    //MPI_DOUBLE,                          // data type of buffer (handle)
    //rank_nf,                             // rank of broadcast root (integer)
    //space->getProcGrid()->getCommBar(3), // ccommunicator (handle) ommunicator (handle)
    //&req_nf);                            // communication request

    //MPI_Wait(&req_1, MPI_STATUS_IGNORE); 
    //MPI_Wait(&req_nf, MPI_STATUS_IGNORE); 

    //ST truncError = 1.;
    //if( u_nf != u_1 ) // just in case u_1=u_nf=0
    //truncError = u_nf / u_1 ;
    //if( NOX::StatusTest::StatusType::Converged == status && truncError < refinementTol ) {
    //if( 0==space->rankST() )
    //std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " < " << refinementTol << "\n\n";
    //break;
    //} else if( 0==space->rankST() )
    //std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " >= " << refinementTol << "\n\n";

    //auto spaceF =
    //Pimpact::RefinementStrategy<SpaceT>::createRefinedSpace(
    //space, Teuchos::tuple<int>( 0, 0, 0, refinementStep ) );

    //auto refineOp =
    //Teuchos::rcp( new Pimpact::TransferCompoundOp<
    //Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp<SpaceT> > >,
    //Pimpact::TransferMultiHarmonicOp< Pimpact::InterpolationOp<SpaceT> >
    //>( space, spaceF ) );
    ////refineOp->print();

    //// init Fields for fine Boundary conditions
    //Teuchos::RCP<CF> xf =
    //createCompoundField(
    //Teuchos::rcp( new VF(spaceF,true) ),
    //Teuchos::rcp( new SF(spaceF) ) );
    ////Teuchos::RCP<CF> fuf =
    ////createCompoundField(
    ////Teuchos::rcp( new VF(spaceF,true) ),
    ////Teuchos::rcp( new SF(spaceF) ) );
    ////xf->getVField().initField( pl->sublist("Base flow") );

    ////Teuchos::RCP<CF> temp =
    ////createCompoundField(
    ////Teuchos::rcp( new VF(spaceF,true) ),
    ////Teuchos::rcp( new SF(spaceF) ) );

    ////refineOp->apply( x->getField(0), *temp );
    //refineOp->apply( x->getField(0), *xf );
    ////refineOp->apply( fu->getField(0), *fuf );

    ////xf->add( 1., *temp, 0., *temp, Pimpact::B::N );

    //x = Pimpact::wrapMultiField( xf );
    ////fu = Pimpact::wrapMultiField( fuf );
    //space = spaceF;
    //}
    //prePostOperators->clear();

    //} // end of for( int refine=0; refine<maxRefinement; ++refine ) {
    //[>****************************************************************************************<]

    //Teuchos::TimeMonitor::summarize();

    //if( 0==space->rankST() ) {
    //pl->sublist("NOX Solver").sublist("Solver Options").remove("Status Test Check Type"); // dirty fix probably, will be fixed in NOX
    //Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );
    //}
  }

  MPI_Finalize();
  return 0;
}
