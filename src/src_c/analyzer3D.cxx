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

using SpaceT = Pimpact::Space<ST, OT, sd, 4, dNC>;


using VF = Pimpact::MultiHarmonicField<Pimpact::VectorField<SpaceT> >;
using SF = Pimpact::MultiHarmonicField<Pimpact::ScalarField<SpaceT> >;

using MVF = Pimpact::MultiField<VF>;
using MSF = Pimpact::MultiField<SF>;

using CF = Pimpact::CompoundField<VF, SF>;
using MF = Pimpact::MultiField<CF>;



int main(int argi, char** argv) {

  // intialize MPI
  MPI_Init(&argi, &argv);

  {
    /////////////////////////////////////////// set up parameters ///////////////////////////
    Teuchos::CommandLineProcessor my_CLP;

    std::string xmlFilename = "parameter3D.xml";
    my_CLP.setOption("filename", &xmlFilename, "file name of the input xml parameterlist");
    int restart = -1;
    my_CLP.setOption("restart", &restart, "number of restart");

    my_CLP.recogniseAllOptions(true);
    my_CLP.throwExceptions(true);

    my_CLP.parse(argi, argv);

    Teuchos::RCP<Teuchos::ParameterList> pl;
    pl = Teuchos::getParametersFromXmlFile("parameterOut.xml");
    ////////////////////////////////////////// end of set up parameters /////////////////////////


    ////////////////////////////////////////// set up initial stuff ////////////////////////////
    Teuchos::RCP<const SpaceT> space =
      Pimpact::create<SpaceT>(Teuchos::sublist(pl, "Space", true));


    if(0 == space->rankST())
      std::cout <<"initial field\n";

    // init vectors
    //Teuchos::RCP<MF> x = Pimpact::wrapMultiField(Pimpact::create<CF>(space));
    Teuchos::RCP<MF> x = Teuchos::rcp(new MF(space));

    // init Fields
    x->getField(0).getVField().initField(pl->sublist("Base flow"));

    auto base = x->getField(0).getVField().get0Field().clone(Pimpact::ECopy::Deep);

    x->getField(0).read(restart);
    /*********************************************************************************/
    Pimpact::ECoord dir = Pimpact::ECoord::Y;
    ST gamma = 10.;
    std::string prefix = "energy_";

    // compute glob energy in y-dir
    if(0==space->si(::Pimpact::F::U, 3)) {
      auto vel = x->getField(0).getVField().get0Field().clone(::Pimpact::ECopy::Deep);
      vel->add(1., *vel, -1., *base);

      auto out = Pimpact::createOstream(prefix + "0.txt",
          space->getProcGrid()->getRankBar(dir));

      Pimpact::computeHEEnergyDir(*vel, *out, gamma);
    }

    for(OT i=std::max(space->si(::Pimpact::F::U, 3), 1); i<=space->ei(::Pimpact::F::U, 3); ++i) {
      {
        auto out = Pimpact::createOstream(prefix + "C"+std::to_string(i) + ".txt",
            space->getProcGrid()->getRankBar(dir));

        Pimpact::computeHEEnergyDir(
            x->getField(0).getVField().getCField(i), *out, gamma);
      }
      {
        auto out = Pimpact::createOstream(prefix + "S"+std::to_string(i) + ".txt",
            space->getProcGrid()->getRankBar(dir));

        Pimpact::computeHEEnergyDir(
            x->getField(0).getVField().getSField(i), *out, gamma);
      }
    }

    x->write();
    Pimpact::writeMHLambda2evol(x->getField(0).getVField(), 10000);
    Pimpact::writeMHLambda2(x->getField(0).getVField(), 20000);
  }

  MPI_Finalize();
  return 0;
}
