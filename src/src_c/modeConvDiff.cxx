//#include <cmath>

#include <mpi.h>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Space.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_CoarsenStrategy.hpp"



using ST = double;
using OT = int;

const int sd = 2;
const int dNC = 4;
//const int dNC = 3;
//const int dNC = 2;

using SpaceT = Pimpact::Space<ST,OT,sd,4,dNC>;

using FSpaceT = Pimpact::Space<ST,OT,sd,4,dNC>;
using CSpaceT = Pimpact::Space<ST,OT,sd,4,2>;

//using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;
using CS = Pimpact::CoarsenStrategy<FSpaceT,CSpaceT>; // dirty fix: till gather isn't fixed

//using VF = Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> >;
//using SF = Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> >;

//using MVF = Pimpact::MultiField<VF>;
//using MSF = Pimpact::MultiField<SF>;

//using CF = Pimpact::CompoundField< VF, SF>;
//using MF = Pimpact::MultiField<CF>;

//using BOp = Pimpact::OperatorBase<MF>;

//using NV = NOX::Pimpact::Vector<MF>;


template<class T1,class T2>
using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T>
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T>
using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;


template<class T>
using MOP = Pimpact::InverseOp<T>;

////template<class T> using PrecS = Pimpact::MultiOpSmoother< Pimpact::DivGradO2JSmoother<T> >;
////template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2SORSmoother >;
//template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;
////template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::Chebyshev >;
//template<class T> using POP2  = Pimpact::PrecInverseOp< T, ConvDiffJT >;
//template<class T> using POP3  = Pimpact::PrecInverseOp< T, ConvDiffSORT >;




int main( int argi, char** argv ) {

  // intialize MPI
  MPI_Init( &argi, &argv );
  {
    /////////////////////////////////////////// set up parameters ///////////////////////////
    Teuchos::CommandLineProcessor my_CLP;

    std::string xmlFilename = "parameter3D.xml";
    my_CLP.setOption("filename", &xmlFilename, "file name of the input xml parameterlist");

    int realCase = 0;
    my_CLP.setOption("realCase", &realCase, "real case");

    my_CLP.recogniseAllOptions(true);
    my_CLP.throwExceptions(true);

    my_CLP.parse(argi,argv);

    auto pl = Teuchos::getParametersFromXmlFile( xmlFilename );
    ////////////////////////////////////////// end of set up parameters /////////////////////////

    int withoutput=pl->sublist("Solver").get<int>( "withoutput", 1 );

    Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(
        Teuchos::sublist( pl, "Space", true ) );

    Pimpact::ModeField< Pimpact::VectorField<SpaceT> > x( space );
    Pimpact::ModeField< Pimpact::VectorField<SpaceT> > y( space );
    Pimpact::ModeField< Pimpact::VectorField<SpaceT> > rhs( space );
    Pimpact::ModeField< Pimpact::VectorField<SpaceT> > sol( space );
    Pimpact::ModeField< Pimpact::VectorField<SpaceT> > err( space );

    ST error;

    Teuchos::RCP<ConvDiffOpT<SpaceT>> zeroOp = Pimpact::create<ConvDiffOpT>( space );

    auto modeOp = Teuchos::rcp(
        new Pimpact::ModeNonlinearOp< ConvDiffOpT<SpaceT> >( zeroOp ) );

    pl->sublist("M_ConvDiff").sublist("Solver").set< Teuchos::RCP<std::ostream> >(
        "Output Stream",
        Pimpact::createOstream( modeOp->getLabel()+".txt", space->rankST() ) );
    auto modeInv = Pimpact::createInverseOp( modeOp, Teuchos::sublist(pl, "M_ConvDiff") );

    auto mgSpaces =
      Pimpact::createMGSpaces<CS>( space, pl->sublist("Multi Grid").get<int>("maxGrids") );

    //pl->sublist("ConvDiff").sublist("Solver").set< Teuchos::RCP<std::ostream> >( "Output Stream",
        //Pimpact::createOstream( zeroOp->getLabel()+".txt", space->rankST() ) );

    //auto mgConvDiff =
      //Pimpact::create<ConvDiffSORT>(
          //zeroOp,
          //Teuchos::sublist(Teuchos::sublist(Teuchos::sublist(pl, "ConvDiff"), "Multi Grid"), "Coarse Grid Solver") ); 
    auto mgConvDiff =
      Pimpact::createMultiGrid<
      Pimpact::VectorField,
      TransVF,
      RestrVF,
      InterVF,
      ConvDiffOpT,
      ConvDiffOpT,
      ConvDiffJT,
      //ConvDiffSORT,
      ConvDiffSORT
      > ( mgSpaces, Teuchos::sublist( Teuchos::sublist( pl, "ConvDiff"), "Multi Grid" ) ) ;

    if( 0==space->rankST() )
      mgConvDiff->print();

    std::string modeConvDiffPrecString =
      pl->sublist("M_ConvDiff").get<std::string>( "preconditioner", "right" );

    auto modePrec =
      Pimpact::createMultiOperatorBase(
        Pimpact::create<Pimpact::EddyPrec>(
          //zeroInv,
          mgConvDiff,
          Teuchos::sublist(Teuchos::sublist(pl, "M_ConvDiff"), "Eddy prec") ) );

    if("right" == modeConvDiffPrecString)
      modeInv->setRightPrec(modePrec);

    if("left" == modeConvDiffPrecString)
      modeInv->setLeftPrec(modePrec);

    ST iRe = 1./space->getDomainSize()->getRe();
    ST a2 = space->getDomainSize()->getAlpha2()*iRe;

    ST lx = space->getDomainSize()->getSize(Pimpact::X) ;
    ST ly = space->getDomainSize()->getSize(Pimpact::Y) ;
    ST lz = space->getDomainSize()->getSize(Pimpact::Z) ;
    std::cout << "a2: " << a2 << "\n";
    std::cout << "iRe: " << iRe << "\n";
    std::cout << "lx: " << lx << "\n";
    std::cout << "ly: " << ly << "\n";
    std::cout << "lz: " << lz << "\n";
    //
    // computing zero mode of z
    // set paramteters
    auto para = Teuchos::parameterList();
    para->set<ST>( "mulI", a2  );
    para->set<ST>( "mulC", 1.  );
    para->set<ST>( "mulL", iRe );
    modeOp->setParameter( para );
    modeInv->setParameter( para );

    // initializtion
    {
      Pimpact::VectorField<SpaceT> wind( space );
      if( realCase )
        wind.initField( pl->sublist("Base flow").sublist("0 mode") );
      else {
        wind(Pimpact::F::U).init( 1. );
        wind(Pimpact::F::V).init( 1. );
        wind(Pimpact::F::W).init( 1. );
      }
      if( withoutput ) wind.write( 1000 );

      zeroOp->assignField( wind );
      mgConvDiff->assignField( wind );
    }

    if( realCase ) {
      rhs.getCField().initField( pl->sublist("Force").sublist("cos mode") );
      rhs.getSField().initField( pl->sublist("Force").sublist("sin mode") );
      if( withoutput ) rhs.write( 100 );
    }
    else{
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
      if( withoutput ) x.write( 10 );

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

      //if( withoutput ) rhs.write( 30 );

      modeOp->apply( sol, y );
      //if( withoutput ) y.write( 20 );
      //if( 2==print ) y.print();
      //if( 3==print ) rhs.print();

      err.add( 1., y, -1., rhs );
      //if( withoutput ) err.write( 0 );
      //if( 1==print ) err.print(   );

      error = err.norm(Pimpact::ENorm::Inf)/rhs.norm(Pimpact::ENorm::Inf);
      std::cout << "\nresidual: " << error << "\n";

      x.init();
    }
    modeInv->apply( rhs, x );

    if( withoutput ) x.write( 10 );

    if( realCase!=0 ) {
      err.add( 1., sol, -1., x );
      if( withoutput ) err.write( 0 );
      //if( print ) err.print(   );

      error = err.norm(Pimpact::ENorm::Inf);
      std::cout << "\nerror: " << error << "\n";
    }

    if( 0==space->rankST() ) {
      //pl->sublist("NOX Solver").sublist("Solver Options").remove("Status Test Check Type"); // dirty fix probably, will be fixed in NOX
      Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );
    }
    Teuchos::TimeMonitor::summarize();
  }
  MPI_Finalize();
  return( 0 );
}
