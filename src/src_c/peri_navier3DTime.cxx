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
#include "NOX_Pimpact_StatusTest.hpp"
#include "NOX_Pimpact_PrePostError.hpp"
#include "NOX_Pimpact_PrePostWriter.hpp"

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
#include "Pimpact_RefinementStrategy.hpp"
#include "Pimpact_TransferCompoundOp.hpp"
#include "Pimpact_TransferMultiHarmonicOp.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_Space.hpp"

#include "Pimpact_DivGradNullSpace.hpp"


// full MG includes
#include "Pimpact_IntResCompoundOp.hpp"
#include "Pimpact_TransferCompoundOp.hpp"
#include "Pimpact_TransferTimeOp.hpp"
#include "Pimpact_TimeStokesBSmoother.hpp"

#include "Pimpact_InterpolationTimeOp.hpp"
#include "Pimpact_RestrictionTimeOp.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"



using ST = double;
using OT = int;

const int sd = 3;
const int dNC = 4;
const int method = 1;
//const int dNC = 3;
//const int dNC = 2;

using SpaceT = Pimpact::Space<ST,OT,sd,4,dNC>;

using FSpaceT = Pimpact::Space<ST,OT,sd,4,dNC>;
using CSpaceT = Pimpact::Space<ST,OT,sd,4,2  >;

//using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;
using CS = Pimpact::CoarsenStrategy<FSpaceT,CSpaceT>; // dirty fix: till gather isn't fixed

using VF = Pimpact::TimeField< Pimpact::VectorField<SpaceT> >;
using SF = Pimpact::TimeField< Pimpact::ScalarField<SpaceT> >;

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
using MOP = Pimpact::InverseOp<T>;

//template<class T> using PrecS = Pimpact::MultiOpSmoother< Pimpact::DivGradO2JSmoother<T> >;
//template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2SORSmoother >;
template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;
//template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::Chebyshev >;
template<class T> using POP2  = Pimpact::PrecInverseOp< T, ConvDiffJT >;
template<class T> using POP3  = Pimpact::PrecInverseOp< T, ConvDiffSORT >;



// full MG
template<class SpaceT>
using CVF = Pimpact::CompoundField<
            Pimpact::TimeField<Pimpact::VectorField<SpaceT> >,
            Pimpact::TimeField<Pimpact::ScalarField<SpaceT> > >;

template<class SpaceT>
using INT = Pimpact::IntResCompoundOp<
            Pimpact::InterpolationTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<SpaceT> > >,
            Pimpact::InterpolationTimeOp<                           Pimpact::InterpolationOp<SpaceT> > >;

template<class SpaceT>
using RES = Pimpact::IntResCompoundOp<
            Pimpact::RestrictionTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<SpaceT>
            > >, Pimpact::RestrictionTimeOp< Pimpact::RestrictionSFOp<SpaceT> > >;

template<class SpaceT1, class SpaceT2> using TCO = Pimpact::TransferCompoundOp<
    Pimpact::TransferTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<SpaceT1, SpaceT2> > >,
    Pimpact::TransferTimeOp<                           Pimpact::TransferOp<SpaceT1, SpaceT2> > >;

using CS4L = Pimpact::CoarsenStrategy<SpaceT,CSpaceT>;

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

    /////////////////////////////////////////// end of set up parameters /////////////////////////


    ///////////////////////////////////////////  set up initial stuff ////////////////////////////

    std::string initGuess = pl->sublist("Solver").get<std::string>( "initial guess", "zero" );

    int withoutput=pl->sublist("Solver").get<int>( "withoutput", 1 );


    Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(
                                         Teuchos::sublist( pl, "Space", true ) );


    if( 0==space->rankST() ) std::cout << "initial field\n";

    // init vectors
    Teuchos::RCP<MF> x = wrapMultiField(
                           createCompoundField(
                             Teuchos::rcp( new VF(space) ),
                             Teuchos::rcp( new SF(space) ) ) ) ;

    // init Fields
    //x->getField(0).getVField().initField( pl->sublist("Base flow") );

    /*********************************************************************************/

    if( 0==space->rankST() ) std::cout << "create operator\n";
    //auto opV2V = Pimpact::create<Pimpact::TimeDtConvectionDiffusionOp<SpaceT,method> >( space );
    //auto opS2V = Pimpact::createTimeOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
    //auto opV2S = Pimpact::createTimeOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );

    //auto op = Pimpact::createCompoundOpWrap(
                //opV2V,
                //opS2V,
                //opV2S );
    auto op = Pimpact::create<Pimpact::TimeNSOp<SpaceT>>( space );

    if( 0==space->rankST() ) std::cout << "create RHS:\n";
    if( 0==space->rankST() ) std::cout << "\tdiv test\n";
    //{
      //opV2S->apply( x->getField(0).getVField(), x->getField(0).getSField()  );
      //ST divergence = x->getField(0).getSField().norm();
      //if( 0==space->rankST() )
        //std::cout << "\n\tdiv(Base Flow): " << divergence << "\n\n";
      //x->getField(0).getSField().init( 0. );
    //}

    std::string rl = "";

    if( 0==space->rankST() ) std::cout << "\tcreate RHS\n";
    auto fu = x->clone( Pimpact::ECopy::Shallow );
    auto sol = fu->clone( Pimpact::ECopy::Shallow);

    if( 0==space->rankST() ) std::cout << "\tBC interpolation\n";
    {
      // to get the Dirichlet for the RHS (necessary interpolation) ugly
      // super ugly hack for BC::Dirichlet
      //opV2V->apply( x->getField(0).getVField(), fu->getField(0).getVField() );
      fu->init( 0., Pimpact::B::N );
    }

    if( 0==space->rankST() ) std::cout << "\tforcing\n";
    // Taylor Green Vortex
    std::string forceType = pl->sublist("Force").get<std::string>("force type","Dirichlet");
    if( "force"== forceType )
      ;
    //fu->getField(0).getVField().initField( pl->sublist("Force"), Pimpact::Add::Y );
    else if( "Taylor-Green"==forceType ) {
      ST pi2 = 2.*std::acos(-1.);
      ST alpha2 = space->getDomainSize()->getAlpha2();
      ST re = space->getDomainSize()->getRe();
      ST A =  pl->sublist("Force").get<ST>("A", 0.5);
      ST B =  pl->sublist("Force").get<ST>("B",-0.5);
      ST C =  pl->sublist("Force").get<ST>("C", 0.);
      ST a =  pl->sublist("Force").get<ST>("a", 1.);
      ST b =  pl->sublist("Force").get<ST>("b", 1.);
      ST c =  pl->sublist("Force").get<ST>("c", 0.);
      TEUCHOS_TEST_FOR_EXCEPT( std::abs( a*A + b*B + c*C )>1.e-16 );

      for( OT i=space->si(Pimpact::F::U,3); i<=space->ei(Pimpact::F::U,3); ++i ) {
        {
          // init RHS
          ST timei = space->getCoordinatesLocal()->getX(Pimpact::F::U,3,i);
          ST timeim =space->getCoordinatesLocal()->getX(Pimpact::F::U,3,i-1);
          //std::cout << time << "\n";
          //ST ctime = (std::cos( timei )+std::cos( timeim ))/2.;
          ST ctime = std::cos( (timei+timeim)/2. );
          ST stime = (std::sin( timei )+std::sin( timeim ))/2.;
          fu->getField(0).getVField()(i)(Pimpact::F::U).initFromFunction(
              [=]( ST x, ST y, ST z ) ->ST {
                return(
                  alpha2/re*A*        std::sin(a*pi2*x)*std::cos(b*pi2*y)*std::cos(c*pi2*z)*ctime	  	        // \alpha^2 dt u
                  +(a*a+b*b+c*c)*A/re*std::sin(a*pi2*x)*std::cos(b*pi2*y)*std::cos(c*pi2*z)*(1.+stime) ); } );// -\lap u

          fu->getField(0).getVField()(i)(Pimpact::F::V).initFromFunction(
              [=]( ST x, ST y, ST z ) ->ST {
                return(
                  alpha2/re*B*        std::cos(a*pi2*x)*std::sin(b*pi2*y)*std::cos(c*pi2*z)*ctime			        // \alpha^2 dt v
                  +(a*a+b*b+c*c)*B/re*std::cos(a*pi2*x)*std::sin(b*pi2*y)*std::cos(c*pi2*z)*(1.+stime) ); } );// -\lap u
        }
        // init sol
        {
          ST time = space->getCoordinatesLocal()->getX(Pimpact::F::U,3,i);
          //std::cout<< time << "\n";
          ST stime = std::sin( time );
          sol->getField(0).getVField()(i)(Pimpact::F::U).initFromFunction(
              [=]( ST x, ST y, ST z ) ->ST {
                return( A*std::sin(a*pi2*x)*std::cos(b*pi2*y)*std::cos(c*pi2*z)*(1.+stime) ); } );
          sol->getField(0).getVField()(i)(Pimpact::F::V).initFromFunction(
              [=]( ST x, ST y, ST z ) ->ST {
                return( B*std::cos(a*pi2*x)*std::sin(b*pi2*y)*std::cos(c*pi2*z)*(1.+stime) ); } );
        }
      }
      sol->getField(0).getVField().changed();
      fu->getField(0).getVField().changed();
    }


    if( 0==space->rankST() ) std::cout << "set initial conditions\n";
    if( "zero"==initGuess )
      x->init( 0. );
    else if( "almost zero"==initGuess ) {
      x->random();
      x->scale(1.e-12);
    } else if( "random"==initGuess )
      x->random();
    else if( "exact"==initGuess || "disturbed"==initGuess ) {
      if( "disturbed"==initGuess ) {
        x->getField(0).getVField().random();
        x->getField(0).getVField().add( 1.e-9, x->getField(0).getVField(), 1., sol->getField(0).getVField() );
      } else
        x->getField(0).getVField() = sol->getField(0).getVField();

      ST pi2 = 2.*std::acos(-1.);
      ST re = space->getDomainSize()->getRe();
      ST A =  pl->sublist("Force").get<ST>("A", 0.5);
      ST B =  pl->sublist("Force").get<ST>("B",-0.5);
      //ST C =  pl->sublist("Force").get<ST>("C", 0.);
      ST a =  pl->sublist("Force").get<ST>("a", 1.);
      ST b =  pl->sublist("Force").get<ST>("b", 1.);
      //ST c =  pl->sublist("Force").get<ST>("c", 1.);

      for( OT i=space->si(Pimpact::F::U,3); i<=space->ei(Pimpact::F::U,3); ++i ) {
        // init
        ST time0 = space->getCoordinatesLocal()->getX(Pimpact::F::U,3,i  );
        ST time1 = space->getCoordinatesLocal()->getX(Pimpact::F::U,3,i-1);
        ST stime0 = std::sin( time0 );
        ST stime1 = std::sin( time1 );
        ST s2time = ( std::pow(1.+stime0,2) + std::pow(1.+stime1,2) )/2.;

        x->getField(0).getSField()(i).initFromFunction(
        [=]( ST x, ST y, ST z ) ->ST {
          return( -1./4.*( A*A*std::cos(2.*a*pi2*x) + B*B*std::cos(2.*b*pi2*y) )*s2time );
        } );
      }
    }

    x->getField(0).getVField().changed();
    x->getField(0).getSField().changed();

    if( withoutput )
      pl->sublist("Picard Solver").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
        "Output Stream", Pimpact::createOstream("Picard.txt", space->rankST() ) );
    else
      pl->sublist("Picard Solver").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
        "Output Stream", Teuchos::rcp( new Teuchos::oblackholestream ) );

    auto opInv = Pimpact::createInverseOp( op, Teuchos::sublist( pl, "Picard Solver" ) );


    /*** init preconditioner *****************************************************************/

    std::string picardPrecString = pl->sublist("Picard Solver").get<std::string>( "preconditioner", "none" );

    //if( false ) {
    //if( "none" != picardPrecString ) {

      //// create Multi space
      //auto mgSpaces = Pimpact::createMGSpaces<CS>( space, pl->sublist("Multi Grid").get<int>("maxGrids") );

      ///////////////////////////////////////////////begin of opv2v////////////////////////////////////

      //if( withoutput )
        //pl->sublist("TimeConvDiff").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
          //"Output Stream", Pimpact::createOstream( opV2V->getLabel()+rl+".txt",
              //space->rankST() ) );
      //else
        //pl->sublist("TimeConvDiff").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
          //"Output Stream", Teuchos::rcp( new Teuchos::oblackholestream ) );

      //auto opV2Vinv = Pimpact::createInverseOp( opV2V,
                      //Teuchos::rcpFromRef( pl->sublist("TimeConvDiff") ) );

      //std::string timeConvDiffPrecString = pl->sublist("TimeConvDiff").get<std::string>( "preconditioner", "none" );

      //if( "none" != timeConvDiffPrecString ) {
        //// creat H0-inv prec
        //auto zeroOp = Pimpact::create<ConvDiffOpT>( space );

        //if( withoutput )
          //pl->sublist("ConvDiff").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
            //"Output Stream",
            //Pimpact::createOstream( zeroOp->getLabel()+rl+".txt",
                                    //space->rankST() ) );
        //else
          //pl->sublist("ConvDiff").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
            //"Output Stream", Teuchos::rcp( new Teuchos::oblackholestream ) );

        //auto zeroInv = Pimpact::createInverseOp(
                         //zeroOp, Teuchos::sublist( pl, "ConvDiff" ) );

        //auto mgConvDiff =
          //Pimpact::createMultiGrid<
          //Pimpact::VectorField,
          //TransVF,
          //RestrVF,
          //InterVF,
          //ConvDiffOpT,
          //ConvDiffOpT,
          ////ConvDiffSORT,
          //ConvDiffJT,
          //MOP
          ////POP2
          ////POP3
          //> ( mgSpaces, Teuchos::sublist( Teuchos::sublist( pl, "ConvDiff"), "Multi Grid" ) ) ;

        ////if( 0==space->rankST() )
        ////mgConvDiff->print();

        //std::string convDiffPrecString = pl->sublist("ConvDiff").get<std::string>( "preconditioner", "none" );
        //if( "right" == convDiffPrecString )
          //zeroInv->setRightPrec( Pimpact::createMultiOperatorBase(mgConvDiff) );
        //if( "left" == convDiffPrecString )
          //zeroInv->setLeftPrec( Pimpact::createMultiOperatorBase(mgConvDiff) );

        ////// create Hinv prec
        //Teuchos::RCP<Pimpact::OperatorBase<MVF> > opV2Vprec =
          //Pimpact::createMultiOperatorBase(
            //Pimpact::createTimeDtConvectionDiffusionPrec(zeroInv)
            ////Pimpact::create<Pimpact::TimeDtConvectionDiffusionPrec>(zeroInv)
            ////Teuchos::rcp( new Pimpact::TimeDtConvectionDiffusionPrec(zeroInv) )
          //);

        //if( "right" == timeConvDiffPrecString )
          //opV2Vinv->setRightPrec( opV2Vprec );
        //if( "left" == timeConvDiffPrecString )
          //opV2Vinv->setLeftPrec( opV2Vprec );
      //}
      /////////////////////////////////////////////end of opv2v//////////////////////////////////////
      ////////--- inverse DivGrad

      //if( withoutput )
        //pl->sublist("DivGrad").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
          //"Output Stream",
          //Pimpact::createOstream( "DivGrad"+rl+".txt", space->rankST() ) );
      //else
        //pl->sublist("DivGrad").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
          //"Output Stream", Teuchos::rcp( new Teuchos::oblackholestream ) );

      //auto divGradOp =
        //Pimpact::createDivGradOp(
          //opV2S->getOperatorPtr(),
          //opS2V->getOperatorPtr() );

      //auto divGradInv2 =
        //Pimpact::createInverseOp<Pimpact::DGProjector>(
          //divGradOp,
          //Teuchos::sublist( pl, "DivGrad") );

      //std::string divGradScalString =
        //pl->sublist("DivGrad").get<std::string>("scaling","none");

      //if( "none" != divGradScalString ) {
        //if( "left" != divGradScalString )
          //divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase(
                                      //Pimpact::createInvDiagonal( divGradOp ) ) );
        //if( "right" != divGradScalString )
          //divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase(
                                       //Pimpact::createInvDiagonal( divGradOp ) ) );
      //}

      //std::string divGradPrecString =
        //pl->sublist("DivGrad").get<std::string>("preconditioner","none");

      //if( "none" != divGradPrecString ) { // init multigrid divgrad

        //auto mgDivGrad = Pimpact::createMultiGrid<
                         //Pimpact::ScalarField,
                         //Pimpact::TransferOp,
                         //Pimpact::RestrictionSFOp,
                         //Pimpact::InterpolationOp,
                         //Pimpact::DivGradOp,
                         //Pimpact::DivGradO2Op,
                         //Pimpact::DivGradO2JSmoother,
                         ////Pimpact::Chebyshev,
                         ////Pimpact::DivGradO2SORSmoother,
                         ////MOP
                         ////Pimpact::Chebyshev
                         ////Pimpact::DivGradO2Inv
                         ////Pimpact::DivGradO2SORSmoother
                         //Pimpact::DivGradO2JSmoother
                         //>( mgSpaces, Teuchos::sublist( Teuchos::sublist( pl, "DivGrad"), "Multi Grid") );

        //if( 0==space->rankST() )
          //mgDivGrad->print();

        //if( "right" == divGradPrecString )
          //divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
        //if( "left" == divGradPrecString )
          //divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
      //}

      //auto divGradInv =
        //Pimpact::createTimeOpWrap(
          //divGradInv2 );

      //// schurcomplement approximator ...
      //auto opSchur =
        //Pimpact::createTripleCompositionOp(
          //opV2S,
          //opV2V,
          //opS2V );

      //auto opS2Sinv =
        //Pimpact::createTripleCompositionOp(
          //divGradInv,
          //opSchur,
          //divGradInv );

      ////if( space->rankST()==0 )
      ////std::cout << opS2Sinv->getLabel() << "\n";

      //auto invTriangOp =
        //Pimpact::createInverseTriangularOp(
          //opV2Vinv,
          //opS2V,
          //opS2Sinv );

      //if( "right" == picardPrecString )
        //opInv->setRightPrec( Pimpact::createMultiOperatorBase( invTriangOp ) );
      //if( "left" == picardPrecString )
        //opInv->setLeftPrec( Pimpact::createMultiOperatorBase( invTriangOp ) );
    //}
    ////if( false ) {
    if( "none" != picardPrecString ) {

      auto mgSpaces = Pimpact::createMGSpaces<CS4L>( space, 10);

      auto mgPL = Teuchos::parameterList();
      //mgPL->sublist("Smoother").set<std::string>("Solver name", "GMRES" );
      //mgPL->sublist("Smoother").set( "Solver",
      //*Pimpact::createLinSolverParameter( "GMRES", 1.e-16, -1,
      //Teuchos::rcp<std::ostream>( new Teuchos::oblackholestream() ), 4 ) );

      mgPL->sublist("Smoother").set<int>( "numIters", 4 );
      mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<int>( "Maximum Iterations", 1000 );
      mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "GMRES" );
      mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver" );
      mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<ST>("Convergence Tolerance" , 1.e-1 );

      auto mg = Pimpact::createMultiGrid<
                CVF,
                TCO,
                RES,
                INT,
                Pimpact::TimeNSOp,
                Pimpact::TimeNSOp,
								Pimpact::TimeNS4DBSmoother,
								//Pimpact::TimeStokesBSmoother,
                MOP > ( mgSpaces, mgPL );

      if( "right" == picardPrecString )
        opInv->setRightPrec( Pimpact::createMultiOperatorBase(mg) );
      //if( "left" == picardPrecString )
      //opInv->setLeftPrec( Pimpact::createMultiOperatorBase( invTriangOp ) );
    }
    // end of init preconditioner *************************************************************

    auto inter = NOX::Pimpact::createInterface(
        fu,
        Pimpact::createMultiOpWrap(op),
        Pimpact::createMultiOpWrap(opInv) );

    auto nx = NOX::Pimpact::createVector( x );

    Teuchos::RCP<NOX::Abstract::Group> group =
      NOX::Pimpact::createGroup(Teuchos::parameterList(), inter, nx);

    // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::Generic> statusTest =
      NOX::StatusTest::buildStatusTests( pl->sublist("NOX Status Test"), NOX::Utils() );

    // Create the solver
    Teuchos::RCP<Teuchos::ParameterList> noxSolverPara =
      Teuchos::sublist(pl, "NOX Solver");

    // NOX PrePostOperators
    Teuchos::RCP<NOX::PrePostOperatorVector> prePostOperators =
      Teuchos::rcp( new NOX::PrePostOperatorVector() );

    prePostOperators->pushBack( 
        Teuchos::rcp( new NOX::Pimpact::PrePostErrorCompute<NV>(Teuchos::sublist(pl, "NOX error"), sol)) );

    prePostOperators->pushBack( 
        Teuchos::rcp(new NOX::Pimpact::PrePostWriter<NV>( Teuchos::sublist(pl, "NOX write") ) ) );

    noxSolverPara->sublist("Solver Options").set<Teuchos::RCP<NOX::Abstract::PrePostOperator>>(
        "User Defined Pre/Post Operator", prePostOperators);

    Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver( group, statusTest, noxSolverPara );


    if( 0==space->rankST() )
      std::cout << "\n\t--- Nf: "<< space->nGlo(3) <<"\tdof: "<<x->getLength()<<"\t---\n";

    // Solve the nonlinear system
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

    x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr();
    if( withoutput ) {
      //Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getField().level();
      //Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getField().write();
      //Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getField().getField(0).getVField().write( 500+1000, true );
      //Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() ))->getField().write( 1000 );
    }

    ST solNorm = sol->getField(0).getVField().norm();
    sol->getField(0).getVField().add( 1., sol->getField(0).getVField(), -1., x->getField(0).getVField() );
    ST error = sol->getField(0).getVField().norm()/solNorm;
    if( 0==space->rankST() ) std::cout << "error: " << error << "\n";
    auto eStream = Pimpact::createOstream("error.txt", space->rankST() );
    *eStream << error << "\n";

    /******************************************************************************************/

    Teuchos::TimeMonitor::summarize();

    if( 0==space->rankST() ) {
      pl->sublist("NOX Solver").sublist("Solver Options").remove("Status Test Check Type"); // dirty fix probably, will be fixed in NOX
      Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );
    }

  }
  MPI_Finalize();
  return( 0 );
}
