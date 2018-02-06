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
#include "NOX_Pimpact_PrePostError.hpp"
#include "NOX_Pimpact_PrePostEnergy.hpp"
#include "NOX_Pimpact_PrePostSpectrum.hpp"
#include "NOX_Pimpact_PrePostWriter.hpp"
#include "NOX_Pimpact_PrePostWriteRestart.hpp"
#include "NOX_Pimpact_RefinementTest.hpp"
#include "NOX_Pimpact_Vector.hpp"

#include "Pimpact_AnalysisTools.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_DivGradProjector.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_ModeNonlinearOp.hpp"
#include "Pimpact_ModeSmoother.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_MultiOpSmoother.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_PicardProjector.hpp"
#include "Pimpact_RefinementStrategy.hpp"
#include "Pimpact_TransferModeOp.hpp"
#include "Pimpact_TransferCompoundOp.hpp"
#include "Pimpact_TransferMultiHarmonicOp.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_Grid.hpp"

#include "Pimpact_DivGradNullSpace.hpp"




// Grid types
using ST = double;
using OT = int;

const int sd = 3;
const int dNC = 4;
//const int dNC = 3;
//const int dNC = 2;

using GridT = Pimpact::Grid<ST, OT, sd, 4, dNC>;

using FGridT = GridT;
using CGridT = Pimpact::Grid<ST, OT, sd, 4, 2>;

using MGGridsT = Pimpact::MGGrids<FGridT, CGridT>;

using CS = Pimpact::CoarsenStrategyGlobal<FGridT, CGridT>;
//using CS = Pimpact::CoarsenStrategy<FGridT, CGridT>;


// Field types
using VF = Pimpact::MultiHarmonicField<Pimpact::VectorField<GridT> >;
using SF = Pimpact::MultiHarmonicField<Pimpact::ScalarField<GridT> >;

template<class T>
using ModeField = Pimpact::ModeField<Pimpact::VectorField<T> >;

using MVF = Pimpact::MultiField<VF>;
using MSF = Pimpact::MultiField<SF>;

using CF = Pimpact::CompoundField<VF, SF>;
using MF = Pimpact::MultiField<CF>;


// Operator types
using OpV2VT = Pimpact::MultiDtConvectionDiffusionOp<GridT>;
using OpV2ST = Pimpact::MultiHarmonicOpWrap<Pimpact::DivOp<GridT> >;
using OpS2VT = Pimpact::MultiHarmonicOpWrap<Pimpact::GradOp<GridT> >;

using OpT = Pimpact::CompoundOpWrap<OpV2VT, OpS2VT, OpV2ST>;
using IOpT = Pimpact::InverseOp<OpT, Pimpact::PicardProjector>;

using BOp = Pimpact::OperatorBase<MF>;

template<class T>
using ModeOp = Pimpact::ModeNonlinearOp<ConvDiffOpT<T> >;


// Smoother types
template<class T>
using MOP = Pimpact::InverseOp<T>;

template<class OpT>
using ModeSmootherT = Pimpact::ModeSmoother<OpT>;


// transfer types
template<class T1, class T2>
using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1, T2> >;
template<class T>
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T>
using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

template<class T1, class T2>
using ModeTransVF = Pimpact::TransferModeOp<Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1, T2> > >;
template<class T>
using ModeRestrVF = Pimpact::TransferModeOp<Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> > >;
template<class T>
using ModeInterVF = Pimpact::TransferModeOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> > >;


// NOX types
using NV = NOX::Pimpact::Vector<MF>;

using InterfaceT = NOX::Pimpact::Interface<MF, Pimpact::MultiOpWrap<OpT>,
      Pimpact::MultiOpWrap<IOpT> >;




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
    int nf_restart = -1;
    my_CLP.setOption("nf_restart", &nf_restart, "number of frequencies nf before restart");

    my_CLP.recogniseAllOptions(true);
    my_CLP.throwExceptions(true);

    my_CLP.parse(argi, argv);

    Teuchos::RCP<Teuchos::ParameterList> pl;
    if(restart==-1)
      pl = Teuchos::getParametersFromXmlFile(xmlFilename);
    else
      pl = Teuchos::getParametersFromXmlFile("parameterOut.xml");
    ////////////////////////////////////////// end of set up parameters /////////////////////////


    //////////////////////////////////////////  set up initial stuff ////////////////////////////

    pl->sublist("NOX Solver").sublist("Solver Options").set<std::string>("Status Test Check Type", "Complete"); // dirty fix probably, should be fixed in NOX

    std::string initGuess = pl->sublist("Solver").get<std::string>("initial guess", "zero");

    int withoutput=pl->sublist("Solver").get<int>("withoutput", 1);

    int maxRefinement     = pl->sublist("Solver").get<int>("max refinement", 1);
    ST  refinementTol  = pl->sublist("Solver").get<ST>("refinement tol",  1.e-6);
    int refinementStep = pl->sublist("Solver").get<int>("refinement step",  2  );

    Teuchos::RCP<const GridT> grid =
      Pimpact::create<GridT>(Teuchos::sublist(pl, "Grid", true));


    if(0==grid->rankST()) std::cout << "initial field\n";

    // init vectors
    Teuchos::RCP<MF> x = Pimpact::wrapMultiField(Pimpact::create<CF>(grid));


    if(restart!=-1) {
      x->getField(0).getVField().read(restart, nf_restart);
      x->getField(0).getSField().read(restart, nf_restart);
    }
    /*********************************************************************************/
    for(int refine=0; refine<maxRefinement; ++refine) {

      if(0==grid->rankST()) std::cout << "create operator\n";

      Teuchos::RCP<OpV2VT> opV2V = Pimpact::createMultiDtConvectionDiffusionOp(grid);

      Teuchos::RCP<OpS2VT> opS2V =
        Pimpact::createMultiHarmonicOpWrap(Pimpact::create<Pimpact::GradOp>(grid));

      Teuchos::RCP<OpV2ST> opV2S =
        Pimpact::createMultiHarmonicOpWrap(Pimpact::create<Pimpact::DivOp>(grid));

      Teuchos::RCP<OpT> op = Pimpact::createCompoundOpWrap(opV2V, opS2V, opV2S);

      std::string rl = "";
      if(maxRefinement>1)
        rl = std::to_string(static_cast<long long>(refine)); // long long needed on brutus(intel)

      if(0==grid->rankST()) std::cout << "create RHS:\n";
      Teuchos::RCP<MF> fu = x->clone(Pimpact::ECopy::Shallow);
      Teuchos::RCP<MF> sol = fu->clone(Pimpact::ECopy::Shallow);

      {
        if(0==grid->rankST()) std::cout << "\tBC interpolation\n";
        {
          Teuchos::RCP<VF> temp = x->getField(0).getVField().clone(Pimpact::ECopy::Shallow);
          temp->initField(pl->sublist("Base flow"));
          // to get the Dirichlet for the RHS (necessary interpolation) ugly
          // super ugly hack for BC::Dirichlet
          opV2V->apply(*temp, fu->getField(0).getVField());
          fu->init(0., Pimpact::B::N);
        }

        if(0==grid->rankST()) std::cout << "\tforcing\n";
        // Taylor-Green Vortex
        std::string forceType = pl->sublist("Force").get<std::string>("force type", "Dirichlet");
        if("force"== forceType)
          fu->getField(0).getVField().initField(pl->sublist("Force"), Pimpact::Add::Y);
        else if("Taylor-Green"==forceType) {
          ST pi2 = 2.*std::acos(-1.);
          ST alpha2 = grid->getDomainSize()->getAlpha2();
          ST re = grid->getDomainSize()->getRe();
          ST A =  pl->sublist("Force").get<ST>("A", 0.5);
          ST B =  pl->sublist("Force").get<ST>("B", -0.5);
          ST a =  pl->sublist("Force").get<ST>("a", 1.);
          ST b =  pl->sublist("Force").get<ST>("b", 1.);
          TEUCHOS_TEST_FOR_EXCEPT(std::abs(a*A + b*B)>1.e-16);

          // --- init RHS ---
          if(0==grid->si(Pimpact::F::U, 3)) { 
            fu->getField(0).getVField().get0Field()(Pimpact::F::U).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return A*(a*a+b*b)*std::cos(a*x*pi2)*std::sin(b*y*pi2)/re; });
            fu->getField(0).getVField().get0Field()(Pimpact::F::V).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return B*(a*a+b*b)*std::sin(a*x*pi2)*std::cos(b*y*pi2)/re; });
          }

          if(1>=grid->si(Pimpact::F::U, 3) && 1<=grid->ei(Pimpact::F::U, 3)) {
            fu->getField(0).getVField().getCField(1)(Pimpact::F::U).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return alpha2*A*std::cos(a*x*pi2)*std::sin(b*y*pi2)/re; });
            fu->getField(0).getVField().getCField(1)(Pimpact::F::V).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return alpha2*B*std::sin(a*x*pi2)*std::cos(b*y*pi2)/re; });

            fu->getField(0).getVField().getSField(1)(Pimpact::F::U).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return A*(a*a+b*b)*std::cos(a*x*pi2)*std::sin(b*y*pi2)/re; });
            fu->getField(0).getVField().getSField(1)(Pimpact::F::V).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return B*(a*a+b*b)*std::sin(a*x*pi2)*std::cos(b*y*pi2)/re; });
          }

          // --- init solution ---
          if(0==grid->si(Pimpact::F::U, 3)) {
            sol->getField(0).getVField().get0Field()(Pimpact::F::U).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return A*std::cos(a*x*pi2)*std::sin(b*y*pi2); });
            sol->getField(0).getVField().get0Field()(Pimpact::F::V).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return B*std::sin(a*x*pi2)*std::cos(b*y*pi2); });
          }

          if(1>=grid->si(Pimpact::F::U, 3) && 1<=grid->ei(Pimpact::F::U, 3)) {
            sol->getField(0).getVField().getSField(1)(Pimpact::F::U).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return A*std::cos(a*x*pi2)*std::sin(b*y*pi2); });
            sol->getField(0).getVField().getSField(1)(Pimpact::F::V).initFromFunction(
                [&](ST x, ST y, ST z) ->ST { return B*std::sin(a*x*pi2)*std::cos(b*y*pi2); });
          }
        }   

        if(0==grid->rankST()) std::cout << "set initial conditions\n";
        if(0==refine && restart==-1) {
          if("zero"==initGuess)
            x->init(0.);
          else if("almost zero"==initGuess) {
            x->random();
            x->scale(1.e-32);
          } else if("random"==initGuess) {
            x->random();
          }
          else if("exitor"==initGuess) {
            for(OT i=std::max(grid->si(Pimpact::F::U, 3), 1); i<=grid->ei(Pimpact::F::U, 3); ++i) {
              x->getField(0).getVField().getField(i).random();
              x->getField(0).getVField().getField(i).scale(0.1);
            }
          }
          else if("exact"==initGuess || "disturbed"==initGuess) {
            if("disturbed"==initGuess) {
              x->getField(0).getVField().random();
              x->getField(0).getVField().add(1.e-9, x->getField(0).getVField(), 1., sol->getField(0).getVField());
            } else
              x->getField(0).getVField() = sol->getField(0).getVField();

            ST pi2 = 2.*std::acos(-1.);
            //ST alpha2 = grid->getDomainSize()->getAlpha2();
            //ST re = grid->getDomainSize()->getRe();
            ST A =  pl->sublist("Force").get<ST>("A", 0.5);
            ST B =  pl->sublist("Force").get<ST>("B", -0.5);
            ST a =  pl->sublist("Force").get<ST>("a", 1.);
            ST b =  pl->sublist("Force").get<ST>("b", 1.);

            if(0==grid->si(Pimpact::F::U, 3)) {
              x->getField(0).getSField().get0Field().initFromFunction(
                  [&](ST x, ST y, ST z) ->ST {
                  return -3./8.*(A*A*std::cos(2.*a*pi2*x) + B*B*std::cos(2.*b*pi2*y)); });
            }
            if(1>=grid->si(Pimpact::F::U, 3) && 1<=grid->ei(Pimpact::F::U, 3)) {
              x->getField(0).getSField().getSField(1).initFromFunction(
                  [&](ST x, ST y, ST z) ->ST {
                  return -1./2.*(A*A*std::cos(2.*a*pi2*x) + B*B*std::cos(2.*b*pi2*y)); });
            }
            if(2>=grid->si(Pimpact::F::U, 3) && 2<=grid->ei(Pimpact::F::U, 3)) {
              x->getField(0).getSField().getCField(2).initFromFunction(
                  [&](ST x, ST y, ST z) ->ST {
                  return +1./8.*(A*A*std::cos(2.*a*pi2*x) + B*B*std::cos(2.*b*pi2*y)); });
            }
          }
          else
            x->getField(0).getVField().initField(pl->sublist("Base flow"));
          x->getField(0).getVField().changed();
          x->getField(0).getSField().changed();
        }
      }

      if(0==grid->rankST()) std::cout << "\tdiv test\n";
      {
        Teuchos::RCP<SF> tempField = x->getField(0).getSField().clone();
        opV2S->apply(x->getField(0).getVField(), *tempField);
        ST divergence = tempField->norm(Pimpact::ENorm::L2);
        if(0==grid->rankST())
          std::cout << "\n\tdiv(Base Flow): " << divergence << "\n\n";
      }


      pl->sublist("Picard Solver").sublist("Solver").set<Teuchos::RCP<std::ostream> >(
          "Output Stream",
          Pimpact::createOstream("Picard"+rl+".txt", withoutput?grid->rankST():-1,
            restart));

      auto opInv = Pimpact::createInverseOp<Pimpact::PicardProjector>(
            op, Teuchos::sublist(pl, "Picard Solver"));

      /*** init scaling *****************************************************************/
      //Teuchos::RCP<MF> scaleField = x->clone(Pimpact::ECopy::Shallow);
      //const ST pi = 4.*std::atan(1.);
      //const ST width = 0.925;
      //const ST eps = 1.e-6;

      //auto scalefunc =
        //[=](ST x, ST y, ST z) ->ST { return (y<=width)?1.:((1.-eps)*std::cos(
              //pi*(y-width)/(1.-width))/2. + 0.5+eps/2.); };

      //if(0==grid->si(Pimpact::F::U, 3)) { 
        //for(Pimpact::F i=Pimpact::F::U; i<GridT::sdim; ++i)
          //scaleField->getField(0).getVField().get0Field()(i).initFromFunction(scalefunc);
        //scaleField->getField(0).getSField().get0Field().initFromFunction(scalefunc);
      //}
      //for(OT i=std::max(grid()->si(Pimpact::F::U, 3), 1); i<=grid()->ei(Pimpact::F::U, 3); ++i) {
        //for(Pimpact::F f=Pimpact::F::U; f<GridT::sdim; ++f) {
          //scaleField->getField(0).getVField().getCField(i)(f).initFromFunction(scalefunc);
          //scaleField->getField(0).getVField().getSField(i)(f).initFromFunction(scalefunc);
        //}
        //scaleField->getField(0).getSField().getCField(i).initFromFunction(scalefunc);
        //scaleField->getField(0).getSField().getSField(i).initFromFunction(scalefunc);
      //}

      /*** init preconditioner **********************************************************/

      std::string picardPrecString =
        pl->sublist("Picard Solver").get<std::string>("preconditioner", "none");

      if("none" != picardPrecString) {
        if(0==grid->rankST()) std::cout << "\tinit Picard preconditioner\n";

        // create Multi grid
        Teuchos::RCP<const MGGridsT> mgGrids = Pimpact::createMGGrids<CS>(
            grid, pl->sublist("Multi Grid").get<int>("maxGrids"));

        ///////////////////////////////////////////begin of opv2v///////////////////////////////////
        //// creat F0-inv prec
        if(0==grid->rankST()) std::cout << "\tinit ConvDiff preconditioner\n";
        auto zeroOp = Pimpact::create<ConvDiffOpT>(grid);

        pl->sublist("ConvDiff").sublist("Solver").set<Teuchos::RCP<std::ostream> >(
            "Output Stream",
            Pimpact::createOstream(zeroOp->getLabel()+rl+".txt",
              withoutput?grid->rankST():-1, restart));

        auto zeroInv = Pimpact::createInverseOp(
            zeroOp, Teuchos::sublist(pl, "ConvDiff"));

        if(0==grid->rankST()) std::cout << "\tinit mgConvDiff preconditioner\n";
        auto mgConvDiff = Pimpact::createMultiGrid<
          Pimpact::VectorField,
          TransVF,
          RestrVF,
          InterVF,
          ConvDiffOpT,
          ConvDiffOpT,
          ConvDiffSORT,
          //ConvDiffJT,
          ConvDiffSORT > (
              mgGrids,
              zeroOp,
              Teuchos::sublist(Teuchos::sublist(pl, "ConvDiff"), "Multi Grid")) ;

        if(0==grid->rankST())
          mgConvDiff->print();

        std::string convDiffPrecString =
          pl->sublist("ConvDiff").get<std::string>("preconditioner", "right");
        if("right" == convDiffPrecString)
          zeroInv->setRightPrec(Pimpact::createMultiOperatorBase(mgConvDiff));
        if("left" == convDiffPrecString)
          zeroInv->setLeftPrec(Pimpact::createMultiOperatorBase(mgConvDiff));

        //// creat FMode-inv prec
        if(0==grid->rankST()) std::cout << "\tinit ModeConvDiff preconditioner\n";
        auto modeOp = Teuchos::rcp(new Pimpact::ModeNonlinearOp<ConvDiffOpT<GridT> >(
              zeroOp));

        pl->sublist("M_ConvDiff").sublist("Solver").set<Teuchos::RCP<std::ostream> >(
            "Output Stream",
            Pimpact::createOstream(modeOp->getLabel()+rl+".txt",
              withoutput?grid->rankST():-1, restart));

        auto modeInv = Pimpact::createInverseOp(modeOp, Teuchos::sublist(pl,
              "M_ConvDiff"));


        Teuchos::RCP<Pimpact::OperatorBase<
          Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<GridT> > > > >
          modePrec;

        //auto
        if(pl->sublist("M_ConvDiff").get<std::string>("preconditioner type", "block"
            )=="block")
          modePrec = Pimpact::createMultiOperatorBase(
              Pimpact::create<Pimpact::ModePrec>(
                mgConvDiff,
                //zeroInv,
                Teuchos::sublist(Teuchos::sublist(pl, "M_ConvDiff"), "Mode prec")));
        else { // precondtioner type =="multi grid"

          auto mgMConvDiff = Pimpact::createMultiGrid<
            ModeField,
            ModeTransVF,
            ModeRestrVF,
            ModeInterVF,
            ModeOp,
            ModeOp,
            ModeSmootherT,
            //ModeSmootherT
            //MOP,
            MOP 
              > (mgGrids, modeOp, Teuchos::sublist(Teuchos::sublist(pl, "M_ConvDiff"), "Multi Grid")) ;

          if(0==grid->rankST())
            mgMConvDiff->print();
          modePrec = Pimpact::createMultiOperatorBase(mgMConvDiff);
        }

        std::string modeConvDiffPrecString = pl->sublist("M_ConvDiff").get<std::string>(
            "preconditioner", "right");

        if("right" == modeConvDiffPrecString)
          modeInv->setRightPrec(modePrec);

        if("left" == modeConvDiffPrecString)
          modeInv->setLeftPrec(modePrec);

        // create Finv prec
        auto opV2Vprec = Pimpact::createMultiHarmonicDiagOp(zeroInv, modeInv);
        //auto opV2Vprec = Pimpact::createMultiHarmonicDiagOp(mgConvDiff, modeInv);

        /////////////////////////////////////////end of opv2v////////////////////////////
        ////--- inverse DivGrad
        if(0==grid->rankST()) std::cout << "\tinit DivGrad preconditioner\n";

        pl->sublist("DivGrad").sublist("Solver").set<Teuchos::RCP<std::ostream> >(
            "Output Stream",
            Pimpact::createOstream("DivGrad"+rl+".txt", withoutput?grid->rankST():-1,
              restart));

        auto divGradOp =
          Pimpact::createDivGradOp(
            opV2S->getOperatorPtr(),
            opS2V->getOperatorPtr());

        auto divGradInv2 =
          Pimpact::createInverseOp<Pimpact::DGProjector>(
            divGradOp,
            Teuchos::sublist(pl, "DivGrad"));

        std::string divGradScalString =
          pl->sublist("DivGrad").get<std::string>("scaling", "none");

        if("none" != divGradScalString) {
          if("left" != divGradScalString)
            divGradInv2->setLeftPrec(Pimpact::createMultiOperatorBase(
                                        Pimpact::createInvDiagonal(divGradOp)));
          if("right" != divGradScalString)
            divGradInv2->setRightPrec(Pimpact::createMultiOperatorBase(
                                         Pimpact::createInvDiagonal(divGradOp)));
        }

        std::string divGradPrecString =
            pl->sublist("DivGrad").get<std::string>("preconditioner", "none");

        if("none" != divGradPrecString) { // init multigrid divgrad

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
            //MOP
            //Pimpact::Chebyshev
            //Pimpact::DivGradO2Inv
            //Pimpact::DivGradO2SORSmoother
            Pimpact::DivGradO2JSmoother
              //Pimpact::DivGradO2Inv
              >(mgGrids,
                  divGradOp,
                  Teuchos::sublist(Teuchos::sublist(pl, "DivGrad"), "Multi Grid"));

          if(0==grid->rankST())
            mgDivGrad->print();

          if("right" == divGradPrecString)
            divGradInv2->setRightPrec(Pimpact::createMultiOperatorBase(mgDivGrad));
          if("left" == divGradPrecString)
            divGradInv2->setLeftPrec(Pimpact::createMultiOperatorBase(mgDivGrad));
        }

        auto divGradInv =
          Pimpact::createMultiHarmonicMultiOpWrap(
            divGradInv2);

        // schurcomplement approximator ...
        auto opSchur =
          Pimpact::createTripleCompositionOp(
            opV2S,
            opV2V,
            opS2V);

        auto opS2Sinv =
          Pimpact::createTripleCompositionOp(
            divGradInv,
            opSchur,
            divGradInv);

        //if(grid->rankST()==0)
        //std::cout << opS2Sinv->getLabel() << "\n";

        auto invTriangOp =
          Pimpact::createInverseTriangularOp(
            opV2Vprec,
            opS2V,
            opS2Sinv);

        if("right" == picardPrecString)
          opInv->setRightPrec(Pimpact::createMultiOperatorBase(invTriangOp));
        if("left" == picardPrecString)
          opInv->setLeftPrec(Pimpact::createMultiOperatorBase(invTriangOp));
      }
      //** end of init preconditioner ***************************************************

      auto inter = NOX::Pimpact::createInterface(
          fu,
          Pimpact::createMultiOpWrap(op),
          Pimpact::createMultiOpWrap(opInv));
          //,
          //scaleField);

      auto nx = NOX::Pimpact::createVector(x);

      Teuchos::RCP<NOX::Abstract::Group> group =
        NOX::Pimpact::createGroup(Teuchos::parameterList(), inter, nx);

      { // setting up refinement stopping cirtion
        Teuchos::RCP<std::ostream> refOut = Pimpact::createOstream(
            "refinementTest.txt", withoutput?grid->rankST():-1, restart);

        Teuchos::RCP<NOX::StatusTest::Generic> refinementTest =
          Teuchos::rcp(new NOX::Pimpact::RefinementTest<InterfaceT>(
                pl->sublist("Solver").get<double>("refinement residual tol", 1.),
                pl->sublist("NOX Status Test").sublist("Test 0").get<double>("Tolerance", 1.e-6),
                refOut));

        pl->sublist("NOX Status Test").sublist("Test 4").set("Test Type", "User Defined");
        pl->sublist("NOX Status Test").sublist("Test 4").set("User Status Test", refinementTest);
      }

      // Set up the status tests
      Teuchos::RCP<NOX::StatusTest::Generic> statusTest =
        NOX::StatusTest::buildStatusTests(pl->sublist("NOX Status Test"), NOX::Utils());

      // Create the solver
      Teuchos::RCP<Teuchos::ParameterList> noxSolverPara =
        Teuchos::sublist(pl, "NOX Solver");

      //pl->sublist("Printing").remove("MyPID");
      noxSolverPara->sublist("Printing").set<Teuchos::RCP<std::ostream> >(
          "Output Stream",
          Pimpact::createOstream("nonlinear"+rl+".txt", grid->rankST(), restart));

      // NOX PrePostOperators
      Teuchos::RCP<NOX::PrePostOperatorVector> prePostOperators =
        Teuchos::rcp(new NOX::PrePostOperatorVector());

      prePostOperators->pushBack(
          Teuchos::rcp(
            new NOX::Pimpact::PrePostErrorCompute<NV>(Teuchos::sublist(pl, "NOX error"), sol)));

      //pl->sublist("NOX energy").set<int>("refinement", refine);
      //prePostOperators->pushBack(
          //Teuchos::rcp(new NOX::Pimpact::PrePostEnergyCompute<NV>(
              //Teuchos::sublist(pl, "NOX energy"), base)));

      prePostOperators->pushBack(
          Teuchos::rcp(new NOX::Pimpact::PrePostWriter<NV>(Teuchos::sublist(pl, "NOX write"))));

      prePostOperators->pushBack(
          Teuchos::rcp(new NOX::Pimpact::PrePostWriteRestart<NV>(Teuchos::sublist(pl, "NOX write restart"))));

      pl->sublist("NOX spectrum").set<int>("refinement", refine);
      prePostOperators->pushBack(
          Teuchos::rcp(new NOX::Pimpact::PrePostSpectrum<NV>(Teuchos::sublist(pl, "NOX spectrum"))));

      noxSolverPara->sublist("Solver Options").set<Teuchos::RCP<NOX::Abstract::PrePostOperator>>(
          "User Defined Pre/Post Operator", prePostOperators);

      Teuchos::RCP<NOX::Solver::Generic> solver =
        NOX::Solver::buildSolver(group, statusTest, noxSolverPara);


      if(0==grid->rankST())
        std::cout << "\n\t--- Nf: "<< grid->nGlo(3) << "\tdof: "<< x->getLength()<< "\t---\n";

      // write ParameterList for restart
      if(0==grid->rankST()) {
        pl->sublist("Grid").set<OT>("nf", grid->nGlo(3));
        pl->sublist("NOX Solver").sublist("Solver Options").remove("Status Test Check Type"); // dirty fix probably, should be fixed in NOX
        pl->sublist("NOX Solver").sublist("Solver Options").remove("User Defined Merit Function"); // dirty fix probably, should be fixed in NOX, just needed for restart
        Teuchos::writeParameterListToXmlFile(*pl, "parameterOut.xml");
      }
      // Solve the nonlinear system
      NOX::StatusTest::StatusType status = NOX::StatusTest::StatusType::Failed;
      {
        Teuchos::TimeMonitor LocalTimer(*Teuchos::TimeMonitor::getNewCounter("Pimpact:: Solving Time"));
        try {
          status = solver->solve();
        } catch(std::logic_error & e) {
          std::cout << e.what() << "\n";
        }
      }


      // Get the answer
      *group = solver->getSolutionGroup();

      x = Teuchos::rcp_const_cast<NV>(
          Teuchos::rcp_dynamic_cast<const NV>(group->getXPtr()))->getFieldPtr();


      // spectral refinement
      if(maxRefinement>1) {
          
        ST truncError = Pimpact::truncErrorEstimate(x->getField(0).getVField());

        if(NOX::StatusTest::StatusType::Converged == status && truncError <refinementTol) {
          if(0==grid->rankST())
            std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " <" << refinementTol << "\n\n";
          break;
        } else if(0==grid->rankST())
          std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " >= " << refinementTol << "\n\n";

        auto gridF =
          Pimpact::RefinementStrategy<GridT>::createRefinedGrid(
            grid, Teuchos::tuple<int>(0, 0, 0, refinementStep));

        auto refineOp =
          Teuchos::rcp(new Pimpact::TransferCompoundOp<
              Pimpact::TransferMultiHarmonicOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<GridT> > >,
              Pimpact::TransferMultiHarmonicOp<Pimpact::InterpolationOp<GridT> >
              >(grid, gridF));
        //refineOp->print();

        // init Fields for fine Boundary conditions
        Teuchos::RCP<CF> xf = Pimpact::create<CF>(grid);


        //refineOp->apply(x->getField(0), *temp);
        refineOp->apply(x->getField(0), *xf);
        //refineOp->apply(fu->getField(0), *fuf);

        //xf->add(1., *temp, 0., *temp, Pimpact::B::N);

        x = Pimpact::wrapMultiField(xf);
        //fu = Pimpact::wrapMultiField(fuf);
        grid = gridF;
      }
      prePostOperators->clear();

    } // end of for(int refine=0; refine<maxRefinement; ++refine) {
    /******************************************************************************************/

    Teuchos::TimeMonitor::summarize();

    if(0==grid->rankST()) {
      //pl->sublist("NOX Solver").sublist("Solver Options").remove("Status Test Check Type"); // dirty fix probably, will be fixed in NOX
      Teuchos::writeParameterListToXmlFile(*pl, "parameterOut.xml");
    }
  }
  MPI_Finalize();
  return 0;
}
