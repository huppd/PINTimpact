#include <mpi.h>

#include <ostream>
#include <fstream>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include "Teuchos_Range1D.hpp"
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"
#include "Teuchos_oblackholestream.hpp"

#include "pimpact.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinSolverParameter.hpp"
#include "BelosPimpactAdapter.hpp"

#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_LinearStokes.hpp"
#include "NOX_Pimpact_SimpleLinear.hpp"
#include "NOX_Pimpact_SimpleNonlinear.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Group.hpp"

#include "NOX_Pimpact_StatusTest.hpp"

#include "NOX.H"




auto CompTime = Teuchos::TimeMonitor::getNewCounter("Pimpact:: Solving Time");


int main(int argi, char** argv ) {

  typedef double S;
  typedef int O;

  typedef Pimpact::TimeField< Pimpact::VectorField<S,O,4> > VF;
  typedef Pimpact::TimeField< Pimpact::ScalarField<S,O,4> > SF;
  typedef Pimpact::CompoundField< VF, SF> CF;
  typedef Pimpact::MultiField<VF> MVF;
  typedef Pimpact::MultiField<SF> MSF;
  typedef Pimpact::MultiField<CF> MF;

  typedef Pimpact::ForcingOp<VF> Fo;
  typedef Pimpact::OperatorBase<MF> BOp;

  typedef NOX::Pimpact::Interface<MF> Inter;
  typedef NOX::Pimpact::Vector<typename Inter::Field> NV;

  // intialize MPI
  MPI_Init( &argi, &argv );

  //// get problem values form Comand line
  Teuchos::CommandLineProcessor my_CLP;

  // physical constants
  S re = 1.e0;
  my_CLP.setOption( "re", &re, "Reynolds number" );

  S alpha2 = 1.;
  my_CLP.setOption( "alpha2", &alpha2, "introduced frequency" );

  S rad = 0.1;
  my_CLP.setOption( "radius", &rad, "radius of disk or sig of dipol" );

  S amp = 0.1;
  my_CLP.setOption( "amp", &amp, "amplitude of elongation of disc" );

  S xm = 1./6.;
  my_CLP.setOption( "xm", &xm, "rotation of disc" );

  S ym = 0.5;
  my_CLP.setOption( "ym", &ym, "rotation of disc" );

  // flow type
  int flow = 5;
  my_CLP.setOption( "flow", &flow,
      "Flow type: 0=zero flow, 1=2D Poiseuille flow in x, 2=2D Poiseuille flow in y, 3=2D pulsatile flow in x, 4=2D pulsatile flow in, 5=2D streaming" );

  int forcing = 0;
  my_CLP.setOption( "force", &forcing,
      "forcing, ja?" );


  // domain type
  int domain = 2;
  my_CLP.setOption( "domain", &domain,
      "Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel" );

  // domain size
  S l1 = 6.;
  my_CLP.setOption( "lx", &l1, "length in x-direction" );

  S l2 = 2.;
  my_CLP.setOption( "ly", &l2, "length in y-direction" );

  S l3 = 1.;
  my_CLP.setOption( "lz", &l3, "length in z-direction" );

  int dim = 2;
  my_CLP.setOption( "dim", &dim, "dimension of problem" );

  // grid size
  O n1 = 33;
  my_CLP.setOption( "nx", &n1, "amount of grid points in x-direction: a*2**q+1" );

  O n2 = 33;
  my_CLP.setOption( "ny", &n2, "amount of grid points in y-direction: a*2**q+1" );

  O n3 = 2.;
  my_CLP.setOption( "nz", &n3, "amount of grid points in z-direction: a*2**q+1" );

  O nt = 32.;
  my_CLP.setOption( "nt", &nt, "amount of grid points in time-direction" );

  const bool cny = true;


  // processor grid size
  O np1 = 2;
  my_CLP.setOption( "npx", &np1, "amount of processors in x-direction" );

  O np2 = 1;
  my_CLP.setOption( "npy", &np2, "amount of processors in y-direction" );

  O np3 = 1.;
  my_CLP.setOption( "npz", &np3, "amount of processors in z-direction" );

  O np4 = 2.;
  my_CLP.setOption( "npt", &np4, "amount of processors in z-direction" );

  // solver stuff
  std::string linSolName = "GMRES";
  my_CLP.setOption( "linSolName", &linSolName, "name of the linear solver" );

  std::string nonLinSolName = "Newton";
  my_CLP.setOption( "nonLinSolName", &nonLinSolName , "name of the non linear solver" );

  std::string lineSearchName = "Backtrack";
  my_CLP.setOption( "linesearch", &lineSearchName, "name of the line search" );

  int fixType = 1.;
  my_CLP.setOption( "fixType", &fixType, "type of fixpoint iteration matrix: 1=Newton, 2=Piccard, 3=lin diag... " );

  bool isNewton = true;
  my_CLP.setOption( "newton","piccard", &isNewton, "type of fixpoint iteration matrix: 1=Newton, 2=Piccard, 3=lin diag... " );

  bool leftPrec = true;
  my_CLP.setOption( "leftPrec","rightPrec", &leftPrec, "type of fixpoint iteration matrix: 1=Newton, 2=Piccard, 3=lin diag... " );

  int precType = 0.;
  my_CLP.setOption( "precType", &precType, "type of preconditioners " );

  int maxIter = 10;
  my_CLP.setOption( "maxIter", &maxIter, "maximum iterations" );

  S tolBelos = 1.e-4;
  my_CLP.setOption( "tolBelos", &tolBelos, "tolerance for linear solver" );

  S tolNOX = 1.e-2;
  my_CLP.setOption( "tolNOX", &tolNOX, "tolerance for non-linear solver" );

  S tolSchur = 1.e-4;
  my_CLP.setOption( "tolSchur", &tolSchur, "tolerance for non-linear solver" );

  S tolPrec = 1.e-4;
  my_CLP.setOption( "tolPrec", &tolPrec, "tolerance for non-linear solver" );


  my_CLP.recogniseAllOptions(true);
  my_CLP.throwExceptions(true);

  my_CLP.parse(argi,argv);
  // end of parsing

  // starting with ininializing
  int rank = Pimpact::init_impact_pre();

  // outputs
  Teuchos::RCP<std::ostream> outPar;
  Teuchos::RCP<std::ostream> outLinSolve;
  Teuchos::RCP<std::ostream> outPrec;
  Teuchos::RCP<std::ostream> outSchur;
  //  Teuchos::RCP<std::ostream> outLap2;
  //  Teuchos::RCP<std::ostream> outSchur;

  if(rank==0) {
    outPar   = Teuchos::rcp( new std::ofstream("para_case.txt") );
    outLinSolve  = Teuchos::rcp( new std::ofstream("stats_linSolve.txt") );
    outPrec  = Teuchos::rcp( new std::ofstream("stats_solvPrec.txt") );
    outSchur  = Teuchos::rcp( new std::ofstream("stats_solvSchur.txt") );
    //    outSchur = Teuchos::rcp( new std::ofstream("stats_solvSchur.txt") );
  } else
    //    outPar = Teuchos::rcp( &blackhole, false) ;
    outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;


  *outPar << " \tCrankNicolson=" << cny << "\n";
  *outPar << " \tNewton=" << isNewton << "\n";
  *outPar << " \tflow=" << flow << "\n";
  *outPar << " \tforce=" << forcing << "\n";
  *outPar << " \tdomain=" << domain << "\n";
  *outPar << " \tre=" << re << "\n";
  *outPar << " \talpha2=" << alpha2 << "\n";


  // init space
  auto ds = Pimpact::createDomainSize<S>( l1, l2, l3 );
  ds->print( *outPar );

  auto bc = Pimpact::createBoudaryConditionsGlobal( Pimpact::EDomainType(domain) );
  bc->print( *outPar );

  auto gs = Pimpact::createGridSizeGlobal( n1, n2, n3, nt );
  gs->print( *outPar );
  //  gs->print();

  auto pgs = Pimpact::createProcGridSize<O>( np1, np2, np3, np4 );
  pgs->print( *outPar );

  auto lgs = Pimpact::createGridSizeLocal( gs, pgs );
  //  lgs->print();

  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
  }
  outPar = Teuchos::null;

  // init IMPACT
  Pimpact::init_impact_mid();

  auto pg = Pimpact::createProcGrid<O>( lgs, bc, pgs );
  //  pg->print();

  Pimpact::init_impact_postpost();

  auto fS = Pimpact::createFieldSpace<O,4>();

  auto iS = Pimpact::createScalarIndexSpace<O>();
  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto space = Pimpact::createSpace<O,4>(fS,iS,iIS,fIS,gs,lgs,pgs,pg);
//  space->print();


  // init vectors
  auto x = Pimpact::createMultiField(
      Pimpact::createCompoundField(
          Pimpact::createTimeField< Pimpact::VectorField<S,O,4> >( space ),
          Pimpact::createTimeField< Pimpact::ScalarField<S,O,4> >( space ) ) );


  // init Fields, init and rhs
  Pimpact::initVectorTimeField( x->getFieldPtr(0)->getVFieldPtr(), Pimpact::EFlowType(flow) );
  x->init( 0. );
  //  x->random();
  auto fu   = x->clone(Pimpact::ShallowCopy);
//  fu->init( 0. );

  Teuchos::RCP<VF> force=Teuchos::null;
  Teuchos::RCP<VF> forcem1=Teuchos::null;

  if( 0!=forcing ) {
    force =
        Pimpact::initVectorTimeField(
            x->getConstFieldPtr(0)->getConstVFieldPtr()->clone( Pimpact::ShallowCopy ),
            Pimpact::OscilatingDisc2D,
            xm*l1, ym*l2, rad, amp );
//    force->write(500);
    Pimpact::initVectorTimeField(
        fu->getFieldPtr(0)->getVFieldPtr(),
        Pimpact::OscilatingDisc2DVel,
        xm*l1, ym*l2, rad, amp );
    fu->getFieldPtr(0)->getVFieldPtr()->scale( *force );

    forcem1 = force->clone(Pimpact::ShallowCopy);
    forcem1->init( 1. );
    forcem1->add( 1., *forcem1, -1., *force );

    x->getFieldPtr(0)->getVFieldPtr()->scale( *force );

  }



  auto para = Pimpact::createLinSolverParameter( linSolName, tolBelos, 10 );
  //    auto para = Pimpact::createLinSolverParameter( linSolName, tol, -1 );
//  para->set( "Maximum Iterations", 30000 );
  para->set( "Output Stream", outLinSolve );


  Teuchos::RCP<Fo> forcingOp = Teuchos::null;
  Teuchos::RCP<Fo> forcingm1Op = Teuchos::null;
  if( 0!=forcing ) {
    forcingOp   = Pimpact::createForcingOp( force   );
    forcingm1Op = Pimpact::createForcingOp( forcem1 );
  }

  S pi = 4.*std::atan(1.);
  S idt = ((S)space->nGlo()[3])/2./pi;


  auto dt  = Pimpact::createDtTimeOp<S,O>( alpha2*idt/re );
  auto lap = Pimpact::createTimeOpWrap< Pimpact::HelmholtzOp<S,O,4>, cny >(
      Pimpact::createHelmholtzOp<S,O,4>( 0., 1./re ),
      Pimpact::createVectorField<S,O,4>( space ) );
  auto adv = Pimpact::createTimeOpWrap<Pimpact::Nonlinear<S,O,4>,cny>  (
      Pimpact::createNonlinear<S,O,4>(),
      Pimpact::createVectorField<S,O,4>( space ) );

//  auto opV2V =
//      Pimpact::createAdd3Op(
//          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
//          dt,
//          lap,
//          adv );
  auto opV2V =
      Pimpact::createAdd3Op(
          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
          Pimpact::createCompositionOp(
              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
              forcingm1Op,
              Pimpact::createAdd3Op(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
                  dt,
                  lap,
                  adv ) ),
           forcingOp );

//    auto opS2V = Pimpact::createTimeOpWrap< Pimpact::Grad<S,O,4> >();
  auto opS2V =
      Pimpact::createCompositionOp(
          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
          forcingm1Op,
//          Pimpact::createTimeOpWrap<Pimpact::Grad<S,O,4>,cny>(
          Pimpact::createTimeOpWrap<Pimpact::Grad<S,O,4> >(
              Pimpact::createGradOp<S,O,4>(),
              Pimpact::createVectorField<S,O,4>( space ) ) );

//  auto opV2S = Pimpact::createTimeOpWrap< Pimpact::Div<S,O,4>,cny >(
  auto opV2S = Pimpact::createTimeOpWrap< Pimpact::Div<S,O,4> >(
              Pimpact::createDivOp<S,O,4>(),
              Pimpact::createScalarField<S,O,4>( space ) );
//  auto opV2S = Pimpact::createTimeOpWrap< Pimpact::Div<S,O,4> >();

//  auto opV2S =
//      Pimpact::createCompositionOp(
//          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
//          Pimpact::createTimeOpWrap<Pimpact::Div<S,O,4>,cny>(
////          Pimpact::createTimeOpWrap<Pimpact::Div<S,O,4> >(
//              Pimpact::createDivOp<S,O,4>(),
//              Pimpact::createScalarField<S,O,4>( space ) ),
//          forcingm1Op
//              );

  auto op =
      Pimpact::createMultiOperatorBase<MF>(
          Pimpact::createCompoundOpWrap(
              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
              opV2V,
              opS2V,
              opV2S )
      );


  Teuchos::RCP<BOp> jop;
  {
//    auto opV2V =
//        Pimpact::createAdd3Op(
//            x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
//            dt,
//            lap,
//            Pimpact::createTimeNonlinearJacobian<S,O>(
//                Teuchos::null,
//                isNewton,
//                Pimpact::createVectorField<S,O,4>(space) ) );
    auto opV2V =
        Pimpact::createAdd3Op(
            x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
            Pimpact::createCompositionOp(
                x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
                forcingm1Op,
                Pimpact::createAdd3Op(
                    x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
                    dt,
                    lap,
                    Pimpact::createTimeNonlinearJacobian<S,O,cny>(
                        x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
                        isNewton,
                        Pimpact::createVectorField<S,O,4>( space )) ) ),
             forcingOp );

    jop =
        Pimpact::createMultiOperatorBase<MF>(
            Pimpact::createCompoundOpWrap(
                x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                opV2V,
                opS2V,
                opV2S ) );
  }


  // init lprec
  Teuchos::RCP<BOp> lprec = Teuchos::null;
  auto schurParams = Pimpact::createLinSolverParameter( linSolName , tolSchur );
  //    schurPara->set( "Verbosity", int( Belos::Errors) );
  schurParams->set( "Output Stream", outSchur );

  auto precParams = Pimpact::createLinSolverParameter( "CG", tolPrec );
  //    solverParams->set( "Verbosity", int( Belos::Errors) );
  precParams->set( "Output Stream", outPrec );

  if( 1==precType ) {
    para->set( "Maximum Iterations", 500 );
    if(0==rank) std::cout << "\n\t---\tprec Type(10): linear block Schur complement\t---\n";

    S pi = 4.*std::atan(1.);
    S deltat = 2.*pi/((S)space->nGlo()[3]);
    auto mglap =
        Pimpact::createMultiOperatorBase<MVF>(
            Pimpact::createTimeOpWrap(
                Pimpact::createMGVHelmholtzOp<S,O,4>( alpha2/re/deltat, 1./re, true ) ) );


    auto A = Pimpact::createMultiOperatorBase<MVF>(
        Pimpact::createTimeOpWrap(
            Pimpact::createHelmholtzOp<S,O,4>( alpha2/re/deltat, 1./re ) ) );

    auto prob2 =
        Pimpact::createLinearProblem<MVF>(
            A,
            Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone() ),
            Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone() ),
            precParams, "GMRES" );
//            Teuchos::null, "PCPG" );
    prob2->setLeftPrec( mglap );

    auto opV2Vinv = Pimpact::createInverseOperatorBase( prob2 );

    // opS2Sinv alias schurcomplement ...
    auto opSchur =
           Pimpact::createOperatorBase< MSF >(
               Pimpact::createTripleCompositionOp(
                   Pimpact::createMultiField(fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy)),
                   Pimpact::createMultiField(fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy)),
//                   fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                   Pimpact::createMultiOpWrap(opS2V),
                   opV2Vinv,
                   Pimpact::createMultiOpWrap(opV2S) ) );//

    auto lp_ = Pimpact::createLinearProblem< MSF >( opSchur, Teuchos::null, Teuchos::null, schurParams, "GMRES");
    /// simple prec
    auto lappre = Pimpact::createMultiOperatorBase< MSF >(
        Pimpact::createTimeOpWrap( Pimpact::createMGVDivGradOp<S,O,4>(true) ) );
    lp_->setLeftPrec(lappre);
    lp_->setRightPrec(lappre);
    ///
    auto opS2Sinv = Pimpact::createInverseOperatorBase( lp_ );

    auto invSchur = Pimpact::createInverseTriangularOp(
          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
          x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
          opV2Vinv,
          opS2V,
          opS2Sinv );

    lprec =
        Pimpact::createMultiOperatorBase<MF>(
            invSchur );
  }
  else
  if( 2==precType ) {
      para->set( "Maximum Iterations", 300 );

     if(0==rank) std::cout << "\n\t---\tprec Type(10): linear block commuter Schur complement\t---\n";

//     S pi = 4.*std::atan(1.);
//     S deltat = 2.*pi/((S)space->nGlo()[3]);
     auto helMGV =
         Pimpact::createMultiOperatorBase<MVF>(
             Pimpact::createTimeOpWrap<Pimpact::MGVHelmholtzOp<S,O,4>,false>(
                 Pimpact::createMGVHelmholtzOp<S,O,4>( alpha2*idt/re, 1./re, true ) ) );
//                 Pimpact::createMGVHelmholtzOp<S,O,4>( alpha2*idt/re, 0., true ) ) );


//     auto helOp = Pimpact::createMultiOperatorBase<MVF>(
//         Pimpact::createTimeOpWrap< Pimpact::HelmholtzOp<S,O,4>, false>(
//             Pimpact::createHelmholtzOp<S,O,4>( alpha2*idt/re, 1./re ) ) );
////             Pimpact::createHelmholtzOp<S,O,4>( alpha2*idt/re, 0. ) ) );
//////             Pimpact::createVectorField( space ) ) );
     auto helOp = Pimpact::createMultiOperatorBase<MVF>(opV2V);

     auto helProb =
         Pimpact::createLinearProblem<MVF>(
             helOp,
             Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy) ),
             Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy) ),
             precParams, "GMRES" );
//             precParams, "GCRODR" );
//             precParams, "TFQMR" );
//             precParams, "CG" );
 //            Teuchos::null, "PCPG" );
//     helProb->setLeftPrec( helMGV );

     auto opV2Vinv = Pimpact::createInverseOperatorBase( helProb );

     // opS2Sinv alias schurcomplement ...
     auto opSchur =
            Pimpact::createMultiOperatorBase< MSF >(
                Pimpact::createTripleCompositionOp(
                    fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                    fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                    opS2V,
                    opV2V,
                    opV2S)  );

     //--- inverse DivGrad
     auto divGradPrec =
         Pimpact::createMultiOperatorBase< Pimpact::MultiField<SF> >(
             Pimpact::createTimeOpWrap( Pimpact::createMGVDivGradOp<S,O,4>(true) ) );
     auto divGradProb =
         Pimpact::createLinearProblem< MSF >(
             Pimpact::createMultiOperatorBase< Pimpact::MultiField<SF> >(
                 Pimpact::createTimeOpWrap(
                     Pimpact::createDivGradOp<S,O,4>(/*x->getConstFieldPtr(0)->getConstVFieldPtr()->getConstField(0)*/) ) ),
             Teuchos::null,
             Teuchos::null,
             schurParams,
             "GMRES" );
     divGradProb->setRightPrec( divGradPrec );
     auto divGradInv = Pimpact::createInverseOperatorBase( divGradProb );

     auto opS2Sinv =
            Pimpact::createOperatorBase< MSF >(
                Pimpact::createTripleCompositionOp(
                    Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy) ),
                    Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy) ),
                    divGradInv,
                    opSchur,
                    divGradInv ) );

////     auto lp_ = Pimpact::createLinearProblem< Pimpact::MultiField<SF> >( opSchur, Teuchos::null, Teuchos::null, schurParams, "GMRES" );
////
////     auto opS2Sinv = Pimpact::createInverseOperatorBase( lp_ );

     auto invSchur = Pimpact::createInverseTriangularOp(
           x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
           x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
           opV2Vinv,
           opS2V,
           opS2Sinv );

     lprec =
         Pimpact::createMultiOperatorBase<MF>(
             invSchur );
  }
  else
  if( 3==precType ) {
     if(0==rank) std::cout << "\n\t---\tprec Type(10): linear block commuter Schur complement\t---\n";

     S pi = 4.*std::atan(1.);
     S deltat = 2.*pi/((S)space->nGlo()[3]);
     auto mglap =
         Pimpact::createMultiOperatorBase<MVF>(
             Pimpact::createTimeOpWrap(
                 Pimpact::createMGVHelmholtzOp<S,O,4>( alpha2/re/deltat, 1./re, true ) ) );


     auto A = Pimpact::createMultiOperatorBase<MVF>(
         Pimpact::createTimeOpWrap(
             Pimpact::createHelmholtzOp<S,O,4>( alpha2/re/deltat, 1./re ) ) );

     auto prob2 =
         Pimpact::createLinearProblem<MVF>(
             A,
             Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone() ),
             Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone() ),
//             precParams, "GMRES" );
             precParams, "CG" );
 //            Teuchos::null, "PCPG" );
     prob2->setLeftPrec( mglap );

     auto opV21Vinv = Pimpact::createInverseOperatorBase( prob2 );
    auto opV2V =
        Pimpact::createAdd3Op(
            x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
            Pimpact::createCompositionOp(
                x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
                forcingm1Op,
                Pimpact::createAdd3Op(
                    x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
                    dt,
                    lap,
                    Pimpact::createTimeNonlinearJacobian<S,O,cny>(
                        Teuchos::null,
                        isNewton,
                        Pimpact::createVectorField<S,O,4>( space )) ) ),
             forcingOp );

    auto lp_ =
        Pimpact::createLinearProblem<MVF>(
            Pimpact::createMultiOperatorBase<MVF>(opV2V),
            Teuchos::null,
            Teuchos::null,
            precParams,
//            "GCRODR" );
            linSolName );
    lp_->setRightPrec( opV21Vinv );

    auto opV2Vinv = Pimpact::createInverseOperatorBase<MVF>( lp_ );

     // opS2Sinv alias schurcomplement ...
     auto opSchur =
            Pimpact::createMultiOperatorBase< MSF >(
                Pimpact::createTripleCompositionOp(
                    fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                    fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                    opS2V,
                    opV2V,
                    opV2S)  );

     //--- inverse DivGrad
     auto divGradPrec =
         Pimpact::createMultiOperatorBase< Pimpact::MultiField<SF> >(
             Pimpact::createTimeOpWrap( Pimpact::createMGVDivGradOp<S,O,4>(true) ) );
     auto divGradProb =
         Pimpact::createLinearProblem< MSF >(
             Pimpact::createMultiOperatorBase< Pimpact::MultiField<SF> >(
                 Pimpact::createTimeOpWrap(
                     Pimpact::createDivGradOp<S,O,4>(/*x->getConstFieldPtr(0)->getConstVFieldPtr()->getConstField(0)*/) ) ),
             Teuchos::null,
             Teuchos::null,
             schurParams,
             "GMRES" );
     divGradProb->setRightPrec( divGradPrec );
     auto divGradInv = Pimpact::createInverseOperatorBase( divGradProb );

     auto opS2Sinv =
            Pimpact::createOperatorBase< MSF >(
                Pimpact::createTripleCompositionOp(
                    Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy) ),
                    Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy) ),
                    divGradInv,
                    opSchur,
                    divGradInv ) );

////     auto lp_ = Pimpact::createLinearProblem< Pimpact::MultiField<SF> >( opSchur, Teuchos::null, Teuchos::null, schurParams, "GMRES" );
////
////     auto opS2Sinv = Pimpact::createInverseOperatorBase( lp_ );

     auto invSchur = Pimpact::createInverseTriangularOp(
           x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
           x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
           opV2Vinv,
           opS2V,
           opS2Sinv );

     lprec =
         Pimpact::createMultiOperatorBase<MF>(
             invSchur );
  }
//  else if( 3==precType ) {
//
//    if(0==rank) std::cout << "\n\t---\tprec Type(10): full Schur complement\t---\n";
//
//    auto opV2V =
//        Pimpact::createAdd3Op(
//            x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
//            Pimpact::createCompositionOp(
//                x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
//                forcingm1Op,
//                Pimpact::createAdd3Op(
//                    x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
//                    dt,
//                    lap,
//                    Pimpact::createTimeNonlinearJacobian<S,O,cny>(
//                        Teuchos::null,
//                        isNewton,
//                        Pimpact::createVectorField<S,O,4>( space )) ) ),
//             forcingOp );
//
//    auto lp_ =
//        Pimpact::createLinearProblem<MVF>(
//            Pimpact::createMultiOperatorBase<MVF>(opV2V),
//            Teuchos::null,
//            Teuchos::null,
//            precParams,
////            "GCRODR" );
//            linSolName );
//
//    auto opV2Vinv = Pimpact::createInverseOperatorBase<MVF>( lp_ );
//
//    auto invSchur = Pimpact::createInverseSchurOp(
//        x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
//        x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
//        opV2Vinv,
//        opS2V,
//        opV2S,
//        schurParams );
//
//    lprec =
//        Pimpact::createMultiOperatorBase<MF>(
//            invSchur );
//  }
//  else if( 4==precType ) {
//     if(0==rank) std::cout << "\n\t---\tprec Type(10): linear block diag\t---\n";
//
//     S pi = 4.*std::atan(1.);
//     S deltat = 2.*pi/((S)space->nGlo()[3]);
//     auto prec =
//         Pimpact::createMultiOperatorBase<MVF>(
//             Pimpact::createTimeOpWrap(
//                 Pimpact::createMGVHelmholtzOp<S,O,4>( alpha2/re/deltat, 1./re, true ) ) );
//
//
//     auto A = Pimpact::createMultiOperatorBase<MVF>(
//         Pimpact::createTimeOpWrap(
//             Pimpact::createHelmholtzOp<S,O,4>( alpha2/re/deltat, 1./re ) ) );
//
//     auto prob2 =
//         Pimpact::createLinearProblem<MVF>(
//             A,
//             Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone() ),
//             Pimpact::createMultiField( fu->getConstFieldPtr(0)->getConstVFieldPtr()->clone() ),
//             precParams, "CG" );
//     prob2->setLeftPrec( prec );
//
//     auto opV2Vinv = Pimpact::createInverseOperatorBase( prob2 );
//
//
//     auto invSchur = Pimpact::createInverseSchurOp(
//         x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
//         x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
//         opV2Vinv,
//         opS2V,
//         opV2S,
//         schurParams );
//
//     lprec =
//         Pimpact::createMultiOperatorBase<MF>(
//             invSchur );
//   }

//      if( 10==precType ) {
  //      if(0==rank) std::cout << "\n\t---\tprec Type(10): full Newton iteration Schur complement\t---\n";
  //      auto opV2V =
  //          Pimpact::createMultiOperatorBase<MVF>(
  //              Pimpact::createAdd3Op(
  //                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
  //                  Pimpact::createCompositionOp(
  //                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
  //                      forcingm1Op,
  //                      //                          Pimpact::createAdd3Op(
  //                      //                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
  //                      dtl ),//,
  //                      //                              Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
  //                      //                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(Pimpact::ShallowCopy),
  //                      //                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy), false ) ) ),
  //                      forcingOp ) );
  //      auto lp_ = Pimpact::createLinearProblem<MVF>(
  //          opV2V, Teuchos::null, Teuchos::null, Teuchos::parameterList(), linSolName );
  //      auto opV2Vinv = Pimpact::createInverseOperatorBase<MVF>( lp_ );
  //
  //      auto invSchur = Pimpact::createInverseSchurOp(
  //          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
  //          x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
  //          opV2Vinv,
  //          opS2V,
  //          opV2S);
  //
  //      lprec =
  //          Pimpact::createMultiOperatorBase<MF>(
  //              invSchur );
  //    }
  //    else if( 11==precType ) {
  //      if(0==rank) std::cout << "\n\t---\tprec Type(10): full Newton iteration Schur complement\t---\n";
  //
  //      auto prec = Pimpact::createMultiHarmonicMLEddy<S,O>( space, 20, nt, alpha2, 1./re );
  //
  //      auto opV2V =
  //          Pimpact::createMultiOperatorBase<MVF>(
  //              Pimpact::createAdd3Op(
  //                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
  //                  Pimpact::createCompositionOp(
  //                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
  //                      forcingm1Op,
  //                      //                          Pimpact::createAdd3Op(
  //                      //                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
  //                      dtl ),//,
  //                      //                              Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
  //                      //                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(Pimpact::ShallowCopy),
  //                      //                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy), false ) ) ),
  //                      forcingOp ) );
  //
  //      //     auto lp_ = Pimpact::createLinearProblem<MVF>(
  //      //           opV2V, Teuchos::null, Teuchos::null, Teuchos::parameterList(), linSolName );
  //      //     lp_->setLeftPrec(prec);
  //      //     auto opV2Vinv = Pimpact::createInverseOperatorBase<MVF>( lp_ );
  //
  //      auto invSchur = Pimpact::createInverseSchurOp(
  //          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
  //          x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
  //          prec,
  //          //         opV2Vinv,
  //          opS2V,
  //          opV2S);
  //
  //      lprec =
  //          Pimpact::createMultiOperatorBase<MF>(
  //              invSchur );
  //    }
  //    //  if( 3==precType ){
  //    //
  //    //    if(0==rank) std::cout << "\n\t---\tpreconditioner matrix(3): complex diagonal Newton iteration\t---\n";
  //    //
  //    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<DJMAdv,DtL>, Fo > > JOp;
  //    //
  //    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
  //    //        Pimpact::createMultiOpWrap(
  //    //            Pimpact::createAdd2Op<Pimpact::AddOp<DJMAdv,DtL>,Fo>(
  //    //                Pimpact::createAdd2Op<DJMAdv,DtL>(
  //    //                    Pimpact::createMultiHarmonicDiagNonlinearJacobian<S,O>(
  //    //                        temp->getFieldPtr(0), isNewton ),
  //    //                    Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
  //    //                    temp->getFieldPtr(0)->clone() ) ,
  //    //                Pimpact::createMultiHarmonicOpWrap(
  //    //                    Pimpact::createForcingOp<S,O>( force ) ) ,
  //    //                temp->getConstFieldPtr(0)->clone() ) )
  //    //    );
  //    //    auto prob =
  //    //          Pimpact::createLinearProblem<MVF>(
  //    //              precOp,
  //    //              fu,
  //    //              fu, para, "GMRES" );
  //    //
  //    //      lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
  //    ////      break;
  //    //  }
  //    //  else if( 4==precType ){
  //    //
  //    //    if(0==rank) std::cout << "\n\t---\tpreconditioner  matrix(4): complex diagonal Picard iteration\t---\n";
  //    //
  //    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<DJMAdv,DtL>, Fo > > JOp;
  //    //
  //    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
  //    //        Pimpact::createMultiOpWrap(
  //    //            Pimpact::createAdd2Op<Pimpact::AddOp<DJMAdv,DtL>,Fo>(
  //    //                Pimpact::createAdd2Op<DJMAdv,DtL>(
  //    //                    Pimpact::createMultiHarmonicDiagNonlinearJacobian<S,O>(
  //    //                        temp->getFieldPtr(0), false ),
  //    //                    Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
  //    //                    temp->getFieldPtr(0)->clone() ) ,
  //    //                Pimpact::createMultiHarmonicOpWrap(
  //    //                    Pimpact::createForcingOp<S,O>( force ) ) ,
  //    //                temp->getConstFieldPtr(0)->clone() ) )
  //    //    );
  //    //    auto prob =
  //    //          Pimpact::createLinearProblem<MVF>(
  //    //              precOp,
  //    //              fu,
  //    //              fu, para, "GMRES" );
  //    //
  //    //      lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
  //    //  }
  //    //  else if( 5==precType ){
  //    //
  //    //    if(0==rank) std::cout << "\n\t---\titeration matrix(5): semi-complex diagonal Newton iteration\t---\n";
  //    //
  //    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<JAdv,DtL>, Fo > > JOp;
  //    //
  //    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
  //    //        Pimpact::createMultiOpWrap(
  //    //            Pimpact::createAdd2Op<Pimpact::AddOp<JAdv,DtL>,Fo>(
  //    //                Pimpact::createAdd2Op<JAdv,DtL>(
  //    //                    Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
  //    //                        Pimpact::createNonlinearJacobian<S,O>( temp->getFieldPtr(0)->get0FieldPtr(), isNewton ) ),
  //    //                    Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
  //    //                    temp->getFieldPtr(0)->clone() ) ,
  //    //                Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
  //    //                temp->getConstFieldPtr(0)->clone() ) )
  //    //    );
  //    //    auto prob =
  //    //          Pimpact::createLinearProblem<MVF>(
  //    //              precOp,
  //    //              fu,
  //    //              fu, para, "GMRES" );
  //    //
  //    //      lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
  //    //  }
  //    //  else if( 6==precType ){
  //    //
  //    //    if(0==rank) std::cout << "\n\t---\titeration matrix(6): semi-complex diagonal Picard iteration\t---\n";
  //    //
  //    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<JAdv,DtL>, Fo > > JOp;
  //    //
  //    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
  //    //        Pimpact::createMultiOpWrap(
  //    //            Pimpact::createAdd2Op<Pimpact::AddOp<JAdv,DtL>,Fo>(
  //    //                Pimpact::createAdd2Op<JAdv,DtL>(
  //    //                    Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
  //    //                        Pimpact::createNonlinearJacobian<S,O>( temp->getFieldPtr(0)->get0FieldPtr(), false ) ),
  //    //                    Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
  //    //                    temp->getFieldPtr(0)->clone() ) ,
  //    //                Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
  //    //                temp->getConstFieldPtr(0)->clone() ) )
  //    //    );
  //    //    auto prob =
  //    //          Pimpact::createLinearProblem<MVF>(
  //    //              precOp,
  //    //              fu,
  //    //              fu, para, "GMRES" );
  //    //    lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
  //    //  }
  //    //  else if( 7==precType ){
  //    //
  //    //    if(0==rank) std::cout << "\n\t---\titeration matrix(7): real diagonal Newton iteration\t---\n";
  //    //
  //    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<JAdv,DtL>, Fo > > JOp;
  //    //
  //    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
  //    //        Pimpact::createMultiOpWrap(
  //    //            Pimpact::createAdd2Op<Pimpact::AddOp<JAdv,DtL>,Fo>(
  //    //                Pimpact::createAdd2Op<JAdv,DtL>(
  //    //                    Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
  //    //                        Pimpact::createNonlinearJacobian<S,O>( temp->getFieldPtr(0)->get0FieldPtr(), isNewton) ),
  //    //                    Pimpact::createMultiDtHelmholtz<S,O>( 0., 0., 1./re ),
  //    //                    temp->getFieldPtr(0)->clone() ) ,
  //    //                Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
  //    //                temp->getConstFieldPtr(0)->clone() ) )
  //    //    );
  //    //    auto prob =
  //    //          Pimpact::createLinearProblem<MVF>(
  //    //              precOp,
  //    //              fu,
  //    //              fu, para, "GMRES" );
  //    //    lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
  //    //  }
  //    //  else if( 8==precType ){
  //    //    if(0==rank) std::cout << "\n\t---\titeration matrix(8): real diagonal Picard iteration\t---\n";
  //    //
  //    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<JAdv,DtL>, Fo > > JOp;
  //    //
  //    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
  //    //        Pimpact::createMultiOpWrap(
  //    //            Pimpact::createAdd2Op<Pimpact::AddOp<JAdv,DtL>,Fo>(
  //    //                Pimpact::createAdd2Op<JAdv,DtL>(
  //    //                    Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
  //    //                        Pimpact::createNonlinearJacobian<S,O>( temp->getFieldPtr(0)->get0FieldPtr(), false) ),
  //    //                    Pimpact::createMultiDtHelmholtz<S,O>( 0., 0., 1./re ),
  //    //                    temp->getFieldPtr(0)->clone() ) ,
  //    //                Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
  //    //                temp->getConstFieldPtr(0)->clone() ) )
  //    //    );
  //    //    auto prob =
  //    //          Pimpact::createLinearProblem<MVF>(
  //    //              precOp,
  //    //              fu,
  //    //              fu, para, "GMRES" );
  //    //    lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
  //    //  }
  //    //  else if( 9==precType ) {
  //    //    if(0==rank) std::cout << "\n\t---\titeration matrix(9): linear terms\t---\n";
  //    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< DtL, Fo > > JOp;
  //    //
  //    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
  //    //           Pimpact::createMultiOpWrap(
  //    //               Pimpact::createAdd2Op<DtL,Fo>(
  //    ////                   Pimpact::createAdd2Op<JMAdv,DtL>(
  //    ////                       Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
  //    ////                           x->getConstFieldPtr(0)->getConst0FieldPtr()->clone(), x->getConstFieldPtr(0)->clone(), false ),
  //    //                   Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
  //    ////                           temp->getFieldPtr(0)->clone()  ,
  //    //                   Pimpact::createMultiHarmonicOpWrap(
  //    //                       Pimpact::createForcingOp<S,O>( force ) ) ,
  //    //                   temp->getFieldPtr(0)->clone() ) )
  //    //       );
  //    //    auto prob =
  //    //          Pimpact::createLinearProblem<MVF>(
  //    //              precOp,
  //    //              fu,
  //    //              fu, para, "GMRES" );
  //    //    lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
  //    //  }
  //    //  else {
  //    //    lprec = Teuchos::null;
  //    //  }
  //
  //
  //    if( 10!=fixType ) {
//  para->set( "Output Stream", outLinSolve );
  auto lp_ = Pimpact::createLinearProblem<MF>(
      jop, x->clone(), fu->clone(), para, linSolName );
  if( leftPrec )
    lp_->setLeftPrec( lprec );
  else
    lp_->setRightPrec( lprec );
  auto lp = Pimpact::createInverseOperatorBase<MF>( lp_ );

  jop=lp;
  //    }

  auto inter = NOX::Pimpact::createInterface<MF>( fu, op, jop );

  auto nx = NOX::Pimpact::createVector(x);

  auto bla = Teuchos::parameterList();

  auto group = NOX::Pimpact::createGroup<Inter>( bla, inter, nx );

  // Set up the status tests
  auto statusTest = NOX::Pimpact::createStatusTest( maxIter, tolNOX, tolBelos*1e-4 );

  // Create the list of solver parameters
  auto solverParametersPtr =
      NOX::Pimpact::createNOXSolverParameter( nonLinSolName, lineSearchName );

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver( group, statusTest, solverParametersPtr);

  if(0==rank) std::cout << "\n\t--- Nf: 0\tdof: "<<x->getLength(true)<<"\t---\n";
  // Solve the nonlinear system
  {
    Teuchos::TimeMonitor LocalTimer(*CompTime);
    //    NOX::StatusTest::StatusType status =
    solver->solve();
  }

  //    // Print the parameter list
  //    //    if( nt==nfs ) {
  //    //      if(rank==0) std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
  //    //      if(rank==0) std::cout << "\n" << status << "\n";
  //    //      if(rank==0) solver->getList().print(std::cout);
  //    //    }

  // Get the answer
  *group = solver->getSolutionGroup();

  // Get a summary from the time monitor.
  Teuchos::TimeMonitor::summarize();

  //    // Print the answer if(rank==0) std::cout << "\n" << "-- Final Solution From Solver --" << "\n";
  //
  Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->write();
  Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(1000);
  //
  //  x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr();
  //  x = Teuchos::rcp_const_cast<NV>( group->getXPtr() )->getFieldPtr();



  //  write solution
  //  x->write();

  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>(outPrec)->close();
    Teuchos::rcp_static_cast<std::ofstream>(outSchur)->close();
  }

  MPI_Finalize();
  return( 0 );

}
