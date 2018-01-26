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

const int sd = 2;
const int dNC = 4;
//const int dNC = 3;
//const int dNC = 2;

using SpaceT = Pimpact::Space<ST,OT,sd,4,dNC>;

using FSpaceT = SpaceT;
using CSpaceT = Pimpact::Space<ST,OT,sd,4,2>;

using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;
//using CS = Pimpact::CoarsenStrategy<FSpaceT,CSpaceT>;

using VF = Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> >;
using SF = Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> >;

using MVF = Pimpact::MultiField<VF>;
using MSF = Pimpact::MultiField<SF>;

using CF = Pimpact::CompoundField< VF, SF>;
using MF = Pimpact::MultiField<CF>;


using OpV2VT = Pimpact::MultiDtConvectionDiffusionOp<SpaceT>;
using OpV2ST = Pimpact::MultiHarmonicOpWrap< Pimpact::DivOp<SpaceT> >;
using OpS2VT = Pimpact::MultiHarmonicOpWrap< Pimpact::GradOp<SpaceT> >;

using OpT = Pimpact::CompoundOpWrap<OpV2VT, OpS2VT, OpV2ST>;
using IOpT = Pimpact::InverseOp<OpT, Pimpact::PicardProjector>;

using BOp = Pimpact::OperatorBase<MF>;

using NV = NOX::Pimpact::Vector<MF>;


template<class T1,class T2>
using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T>
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T>
using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;


template<class T>
using MOP = Pimpact::InverseOp<T>;

//template<class T> using PrecS = Pimpact::MultiOpSmoother< Pimpact::DivGradO2JSmoother<T> >;
//template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2SORSmoother >;
template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;
//template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::Chebyshev >;
template<class T> using POP2  = Pimpact::PrecInverseOp< T, ConvDiffJT >;
template<class T> using POP3  = Pimpact::PrecInverseOp< T, ConvDiffSORT >;


using InterfaceT = NOX::Pimpact::Interface< MF, Pimpact::MultiOpWrap<OpT>,
      Pimpact::MultiOpWrap<IOpT> >;





int main( int argi, char** argv ) {

  // intialize MPI
  MPI_Init( &argi, &argv );

  {
    /////////////////////////////////////////// set up parameters ///////////////////////////
    Teuchos::CommandLineProcessor my_CLP;

    std::string xmlFilename = "parameter3D.xml";
    my_CLP.setOption("filename", &xmlFilename, "file name of the input xml parameterlist");

    my_CLP.recogniseAllOptions(true);
    my_CLP.throwExceptions(true);

    my_CLP.parse(argi,argv);

    auto pl = Teuchos::getParametersFromXmlFile( xmlFilename );
    //pl->print();

    ////////////////////////////////////////// end of set up parameters /////////////////////////


    //////////////////////////////////////////  set up initial stuff ////////////////////////////

    std::string initGuess = pl->sublist("Solver").get<std::string>( "initial guess", "zero" );

    int withoutput=pl->sublist("Solver").get<int>( "withoutput", 1 );

    int maxRefinement     = pl->sublist("Solver").get<int>("max refinement", 1);
    ST  refinementTol  = pl->sublist("Solver").get<ST>(  "refinement tol",  1.e-6 );
    int refinementStep = pl->sublist("Solver").get<int>( "refinement step",  2     );

    Teuchos::RCP<const SpaceT> space =
      Pimpact::create<SpaceT>( Teuchos::sublist( pl, "Space", true ) );


    if( 0==space->rankST() ) std::cout << "initial field\n";

    // init vectors
    Teuchos::RCP<MF> x = Pimpact::wrapMultiField( Pimpact::create<CF>(space) );

    // init Fields
    x->getField(0).getVField().initField( pl->sublist("Base flow") );

    auto base = x->getField(0).getVField().get0Field().clone(Pimpact::ECopy::Deep);
    /*********************************************************************************/
    for( int refine=0; refine<maxRefinement; ++refine ) {

      if( 0==space->rankST() ) std::cout << "create operator\n";

      Teuchos::RCP<OpV2VT> opV2V = Pimpact::createMultiDtConvectionDiffusionOp( space );

      Teuchos::RCP<OpS2VT> opS2V =
        Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );

      Teuchos::RCP<OpV2ST> opV2S =
        Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );

      Teuchos::RCP<OpT> op = Pimpact::createCompoundOpWrap( opV2V, opS2V, opV2S );


      if( 0==space->rankST() ) std::cout << "\tdiv test\n";
      {
        auto tempField = x->getField(0).getSField().clone();
        opV2S->apply( x->getField(0).getVField(), *tempField );
        ST divergence = tempField->norm(Pimpact::ENorm::L2);
        if( 0==space->rankST() )
          std::cout << "\n\tdiv(Base Flow): " << divergence << "\n\n";
      }

      std::string rl = "";
      if( maxRefinement>1 )
        rl = std::to_string( static_cast<long long>(refine) ); // long long needed on brutus(intel)

      if( 0==space->rankST() ) std::cout << "create RHS:\n";
      Teuchos::RCP<MF> fu = x->clone( Pimpact::ECopy::Shallow );
      auto sol = fu->clone( Pimpact::ECopy::Shallow);

      {
        if( 0==space->rankST() ) std::cout << "\tBC interpolation\n";
        {
          auto temp = x->getField(0).getVField().clone(Pimpact::ECopy::Shallow);
          temp->initField( pl->sublist("Base flow") );
          // to get the Dirichlet for the RHS (necessary interpolation) ugly
          // super ugly hack for BC::Dirichlet
          opV2V->apply( *temp, fu->getField(0).getVField() );
          fu->init( 0., Pimpact::B::N );
        }

        if( 0==space->rankST() ) std::cout << "\tforcing\n";
        // Taylor-Green Vortex
        std::string forceType = pl->sublist("Force").get<std::string>("force type","Dirichlet");
        if( "force"== forceType )
          fu->getField(0).getVField().initField( pl->sublist("Force"), Pimpact::Add::Y );
        else if( "Taylor-Green"==forceType ) {
          ST pi2 = 2.*std::acos(-1.);
          ST alpha2 = space->getDomainSize()->getAlpha2();
          ST re = space->getDomainSize()->getRe();
          ST A =  pl->sublist("Force").get<ST>("A", 0.5);
          ST B =  pl->sublist("Force").get<ST>("B",-0.5);
          ST a =  pl->sublist("Force").get<ST>("a", 1.);
          ST b =  pl->sublist("Force").get<ST>("b", 1.);
          TEUCHOS_TEST_FOR_EXCEPT( std::abs( a*A + b*B )>1.e-16 );

          // --- init RHS ---
          if( 0==space->si(Pimpact::F::U,3) ) {
            fu->getField(0).getVField().get0Field()(Pimpact::F::U).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return A*(a*a+b*b)*std::sin(a*x*pi2)*std::cos(b*y*pi2)/re; } );
            fu->getField(0).getVField().get0Field()(Pimpact::F::V).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return B*(a*a+b*b)*std::cos(a*x*pi2)*std::sin(b*y*pi2)/re; } );
          }

          if( 1>=space->si(Pimpact::F::U,3) && 1<=space->ei(Pimpact::F::U,3) ) {
            fu->getField(0).getVField().getCField(1)(Pimpact::F::U).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return alpha2*A*std::sin(a*x*pi2)*std::cos(b*y*pi2)/re; } );
            fu->getField(0).getVField().getCField(1)(Pimpact::F::V).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return alpha2*B*std::cos(a*x*pi2)*std::sin(b*y*pi2)/re; } );

            fu->getField(0).getVField().getSField(1)(Pimpact::F::U).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return A*(a*a+b*b)*std::sin(a*x*pi2)*std::cos(b*y*pi2)/re; } );
            fu->getField(0).getVField().getSField(1)(Pimpact::F::V).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return B*(a*a+b*b)*std::cos(a*x*pi2)*std::sin(b*y*pi2)/re; } );
          }

          // --- init solution ---
          if( 0==space->si(Pimpact::F::U,3) ) {
            sol->getField(0).getVField().get0Field()(Pimpact::F::U).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return A*std::sin(a*x*pi2)*std::cos(b*y*pi2); } );
            sol->getField(0).getVField().get0Field()(Pimpact::F::V).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return B*std::cos(a*x*pi2)*std::sin(b*y*pi2); } );
          }

          if( 1>=space->si(Pimpact::F::U,3) && 1<=space->ei(Pimpact::F::U,3) ) {
            sol->getField(0).getVField().getSField(1)(Pimpact::F::U).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return A*std::sin(a*x*pi2)*std::cos(b*y*pi2); } );
            sol->getField(0).getVField().getSField(1)(Pimpact::F::V).initFromFunction(
                [&]( ST x, ST y, ST z ) ->ST { return B*std::cos(a*x*pi2)*std::sin(b*y*pi2); } );
          }
        }

        if( 0==space->rankST() ) std::cout << "set initial conditions\n";
        if( 0==refine ) {
          if( "zero"==initGuess )
            x->init( 0. );
          else if( "almost zero"==initGuess ) {
            x->random();
            x->scale(1.e-32);
          } else if( "random"==initGuess ) {
            x->random();
          }
          else if( "exitor"==initGuess ) {
            for( OT i=std::max(space->si(Pimpact::F::U,3),1); i<=space->ei(Pimpact::F::U,3); ++i ) {
              x->getField(0).getVField().getField(i).random();
              x->getField(0).getVField().getField(i).scale(0.1);
            }
          }
          else if( "exact"==initGuess || "disturbed"==initGuess ) {
            if( "disturbed"==initGuess ) {
              x->getField(0).getVField().random();
              x->getField(0).getVField().add( 1.e-9, x->getField(0).getVField(), 1., sol->getField(0).getVField() );
            } else
              x->getField(0).getVField() = sol->getField(0).getVField();

            ST pi2 = 2.*std::acos(-1.);
            //ST alpha2 = space->getDomainSize()->getAlpha2();
            //ST re = space->getDomainSize()->getRe();
            ST A =  pl->sublist("Force").get<ST>("A", 0.5);
            ST B =  pl->sublist("Force").get<ST>("B",-0.5);
            ST a =  pl->sublist("Force").get<ST>("a", 1.);
            ST b =  pl->sublist("Force").get<ST>("b", 1.);

            if( 0==space->si(Pimpact::F::U,3) ) {
              x->getField(0).getSField().get0Field().initFromFunction(
                  [&]( ST x, ST y, ST z ) ->ST {
                  return -3./8.*( A*A*std::cos(2.*a*pi2*x) + B*B*std::cos(2.*b*pi2*y) ); } );
            }
            if( 1>=space->si(Pimpact::F::U,3) && 1<=space->ei(Pimpact::F::U,3) ) {
              x->getField(0).getSField().getSField(1).initFromFunction(
                  [&]( ST x, ST y, ST z ) ->ST {
                  return -1./2.*( A*A*std::cos(2.*a*pi2*x) + B*B*std::cos(2.*b*pi2*y) ); } );
            }
            if( 2>=space->si(Pimpact::F::U,3) && 2<=space->ei(Pimpact::F::U,3) ) {
              x->getField(0).getSField().getCField(2).initFromFunction(
                  [&]( ST x, ST y, ST z ) ->ST {
                  return +1./8.*( A*A*std::cos(2.*a*pi2*x) + B*B*std::cos(2.*b*pi2*y) ); } );
            }
          }
          x->getField(0).getVField().changed();
          x->getField(0).getSField().changed();
        }
        //fu->write( 999 );
      }

      pl->sublist("Picard Solver").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
          "Output Stream",
          Pimpact::createOstream("Picard"+rl+".txt", withoutput?space->rankST():-1));

      Teuchos::RCP<IOpT> opInv = Pimpact::createInverseOp<Pimpact::PicardProjector>(
            op, Teuchos::sublist( pl, "Picard Solver" ) );

      /*** init preconditioner **********************************************************/

      std::string picardPrecString = pl->sublist("Picard Solver").get<std::string>( "preconditioner", "none" );

      if( "none" != picardPrecString ) {

        // create Multi space
        auto mgSpaces = Pimpact::createMGSpaces<CS>( space, pl->sublist("Multi Grid").get<int>("maxGrids") );

        ///////////////////////////////////////////begin of opv2v////////////////////////////////////
        pl->sublist("MH_ConvDiff").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
            "Output Stream",
            Pimpact::createOstream(opV2V->getLabel()+rl+".txt", withoutput?space->rankST():-1));

        auto opV2Vinv = Pimpact::createInverseOp( opV2V, Teuchos::rcpFromRef(
                          pl->sublist("MH_ConvDiff") ) );

        std::string mhConvDiffPrecString = pl->sublist("MH_ConvDiff").get<std::string>( "preconditioner", "right" );

        if( "none" != mhConvDiffPrecString ) {
          // creat H0-inv prec
          auto zeroOp = Pimpact::create<ConvDiffOpT>( space );

          pl->sublist("ConvDiff").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
              "Output Stream",
              Pimpact::createOstream(
                zeroOp->getLabel()+rl+".txt",
                withoutput?space->rankST():-1));

          auto zeroInv = Pimpact::createInverseOp(
                           zeroOp, Teuchos::sublist( pl, "ConvDiff" ) );

          auto modeOp = Teuchos::rcp(
                  new Pimpact::ModeNonlinearOp< ConvDiffOpT<SpaceT> >( zeroOp ) );

          pl->sublist("M_ConvDiff").sublist("Solver").set<Teuchos::RCP<std::ostream> >(
              "Output Stream",
              Pimpact::createOstream(
                modeOp->getLabel()+rl+".txt",
                withoutput?space->rankST():-1));

          auto modeInv = Pimpact::createInverseOp(
                  modeOp, Teuchos::sublist(pl, "M_ConvDiff") );

          auto mgConvDiff = Pimpact::createMultiGrid<
            Pimpact::VectorField,
            TransVF,
            RestrVF,
            InterVF,
            ConvDiffOpT,
            ConvDiffOpT,
            //ConvDiffSORT,
            ConvDiffJT,
            ConvDiffSORT
            //MOP
              > ( mgSpaces,
                  zeroOp,
                  Teuchos::sublist( Teuchos::sublist( pl, "ConvDiff"), "Multi Grid" ) ) ;

          if( 0==space->rankST() )
            mgConvDiff->print();

          std::string convDiffPrecString =
            pl->sublist("ConvDiff").get<std::string>( "preconditioner", "right" );
          if( "right" == convDiffPrecString )
            zeroInv->setRightPrec( Pimpact::createMultiOperatorBase(mgConvDiff) );
          if( "left" == convDiffPrecString )
            zeroInv->setLeftPrec( Pimpact::createMultiOperatorBase(mgConvDiff) );

          std::string modeConvDiffPrecString =
            pl->sublist("M_ConvDiff").get<std::string>( "preconditioner", "right" );

          auto modePrec =
            Pimpact::createMultiOperatorBase(
                Pimpact::create<Pimpact::ModePrec>(
                  //mgConvDiff,
                  zeroInv,
                  Teuchos::sublist(Teuchos::sublist(pl, "M_ConvDiff"), "Mode prec") ) );

          if("right" == modeConvDiffPrecString)
            modeInv->setRightPrec(modePrec);

          if("left" == modeConvDiffPrecString)
            modeInv->setLeftPrec(modePrec);

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


        /////////////////////////////////////////end of opv2v////////////////////////////////////
        ////--- inverse DivGrad

        pl->sublist("DivGrad").sublist("Solver").set<Teuchos::RCP<std::ostream> >(
            "Output Stream",
            Pimpact::createOstream( "DivGrad"+rl+".txt", withoutput?space->rankST():-1));

        auto divGradOp =
          Pimpact::createDivGradOp(
            opV2S->getOperatorPtr(),
            opS2V->getOperatorPtr() );

        auto divGradInv2 =
          Pimpact::createInverseOp<Pimpact::DGProjector>(
            divGradOp,
            Teuchos::sublist( pl, "DivGrad") );

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
              Pimpact::DivGradO2JSmoother >(
                  mgSpaces,
                  divGradOp,
                  Teuchos::sublist( Teuchos::sublist( pl, "DivGrad"), "Multi Grid") );

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

        //if( space->rankST()==0 )
        //std::cout << opS2Sinv->getLabel() << "\n";

        auto invTriangOp =
          Pimpact::createInverseTriangularOp(
            opV2Vinv,
            opS2V,
            opS2Sinv );

        if( "right" == picardPrecString )
          opInv->setRightPrec( Pimpact::createMultiOperatorBase( invTriangOp ) );
        if( "left" == picardPrecString )
          opInv->setLeftPrec( Pimpact::createMultiOperatorBase( invTriangOp ) );
      }
      //** end of init preconditioner ***********************************************************

      Teuchos::RCP<InterfaceT> inter =
        NOX::Pimpact::createInterface(
            fu, Pimpact::createMultiOpWrap(op), Pimpact::createMultiOpWrap(opInv) );

      auto nx = NOX::Pimpact::createVector( x );

      Teuchos::RCP<NOX::Abstract::Group> group =
        NOX::Pimpact::createGroup(Teuchos::parameterList(), inter, nx);

      { // setting up refinement stopping cirtion
        Teuchos::RCP<std::ostream> refOut = Pimpact::createOstream(
            "refinementTest.txt", space->rankST(), ((refine==0)?-1:1) );

        Teuchos::RCP<NOX::StatusTest::Generic> refinementTest =
          Teuchos::rcp( new NOX::Pimpact::RefinementTest<InterfaceT>(
                pl->sublist("Solver").get<double>( "refinement residual tol", 1.),
                pl->sublist("NOX Status Test").sublist("Test 0").get<double>("Tolerance", 1.e-6),
                refOut) );

        pl->sublist("NOX Status Test").sublist("Test 0").sublist("Test 4").set("Test Type", "User Defined");
        pl->sublist("NOX Status Test").sublist("Test 0").sublist("Test 4").set("User Status Test", refinementTest);
      }

      // Set up the status tests
      Teuchos::RCP<NOX::StatusTest::Generic> statusTest =
        NOX::StatusTest::buildStatusTests( pl->sublist("NOX Status Test"), NOX::Utils() );

      // Create the solver
      Teuchos::RCP<Teuchos::ParameterList> noxSolverPara =
        Teuchos::sublist(pl, "NOX Solver");

      //pl->sublist("Printing").remove("MyPID");
      noxSolverPara->sublist("Printing").set< Teuchos::RCP<std::ostream> >(
          "Output Stream", Pimpact::createOstream("nonlinear"+rl+".txt", space->rankST() ) );

      // NOX PrePostOperators
      Teuchos::RCP<NOX::PrePostOperatorVector> prePostOperators =
        Teuchos::rcp( new NOX::PrePostOperatorVector() );

      prePostOperators->pushBack( 
          Teuchos::rcp(
            new NOX::Pimpact::PrePostErrorCompute<NV>(Teuchos::sublist(pl, "NOX error"), sol)) );

      pl->sublist("NOX energy").set<int>( "refinement", refine );
      prePostOperators->pushBack( 
          Teuchos::rcp( new NOX::Pimpact::PrePostEnergyCompute<NV>(
              Teuchos::sublist(pl, "NOX energy"), base)));

      prePostOperators->pushBack( 
          Teuchos::rcp(new NOX::Pimpact::PrePostWriter<NV>( Teuchos::sublist(pl, "NOX write") ) ) );

      pl->sublist("NOX spectrum").set<int>( "refinement", refine );
      prePostOperators->pushBack( 
          Teuchos::rcp(new NOX::Pimpact::PrePostSpectrum<NV>(Teuchos::sublist(pl, "NOX spectrum"))));

      noxSolverPara->sublist("Solver Options").set<Teuchos::RCP<NOX::Abstract::PrePostOperator>>(
          "User Defined Pre/Post Operator", prePostOperators);

      Teuchos::RCP<NOX::Solver::Generic> solver =
        NOX::Solver::buildSolver( group, statusTest, noxSolverPara );


      if( 0==space->rankST() )
        std::cout << "\n\t--- Nf: "<< space->nGlo(3) <<"\tdof: "<<x->getLength()<<"\t---\n";

      // Solve the nonlinear system
      NOX::StatusTest::StatusType status = NOX::StatusTest::StatusType::Failed;
      {
        Teuchos::TimeMonitor LocalTimer(*Teuchos::TimeMonitor::getNewCounter("Pimpact:: Solving Time"));
        try {
          solver->solve();
        } catch( std::logic_error & e ) {
          std::cout << e.what() << "\n";
        }
      }


      // Get the answer
      *group = solver->getSolutionGroup();

      x = Teuchos::rcp_const_cast<NV>(
          Teuchos::rcp_dynamic_cast<const NV>(group->getXPtr()))->getFieldPtr();


      // spectral refinement
      if( maxRefinement>1 ) {

        ST truncError = Pimpact::truncErrorEstimate( x->getField(0).getVField() );

        if( NOX::StatusTest::StatusType::Converged == status && truncError < refinementTol ) {
          if( 0==space->rankST() )
            std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " < " << refinementTol << "\n\n";
          break;
        } else if( 0==space->rankST() )
          std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " >= " << refinementTol << "\n\n";

        auto spaceF =
          Pimpact::RefinementStrategy<SpaceT>::createRefinedSpace(
            space, Teuchos::tuple<int>( 0, 0, 0, refinementStep ) );

        auto refineOp =
          Teuchos::rcp( new Pimpact::TransferCompoundOp<
              Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp<SpaceT> > >,
              Pimpact::TransferMultiHarmonicOp< Pimpact::InterpolationOp<SpaceT> >
              >( space, spaceF ) );
        //refineOp->print();

        // init Fields for fine Boundary conditions
        Teuchos::RCP<CF> xf = Pimpact::create<CF>(spaceF);
        //xf->getVField().initField( pl->sublist("Base flow") );


        //refineOp->apply( x->getField(0), *temp );
        refineOp->apply( x->getField(0), *xf );
        //refineOp->apply( fu->getField(0), *fuf );

        //xf->add( 1., *temp, 0., *temp, Pimpact::B::N );

        x = Pimpact::wrapMultiField( xf );
        //fu = Pimpact::wrapMultiField( fuf );
        space = spaceF;
      }
      prePostOperators->clear();

    } // end of for( int refine=0; refine<maxRefinement; ++refine ) {
    /******************************************************************************************/

    Teuchos::TimeMonitor::summarize();

    if( 0==space->rankST() ) {
      //pl->sublist("NOX Solver").sublist("Solver Options").remove("Status Test Check Type"); // dirty fix probably, will be fixed in NOX
      Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );
    }
  }
  MPI_Finalize();
  return 0;
}
