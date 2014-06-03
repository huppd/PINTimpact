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
#include "Pimpact_GridSize.hpp"
#include "Pimpact_ProcGridSize.hpp"
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

  typedef Pimpact::MultiHarmonicField< Pimpact::VectorField<S,O> > VF;
  typedef Pimpact::MultiHarmonicField< Pimpact::ScalarField<S,O> > SF;
  typedef Pimpact::CompoundField< VF, SF> CF;
  typedef Pimpact::MultiField<VF> MVF;
  typedef Pimpact::MultiField<CF> MF;

  typedef Pimpact::MultiHarmonicOpWrap< Pimpact::ForcingOp<S,O> > Fo;
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

  //  S drift = 1.;
  //  my_CLP.setOption( "drift", &drift, "phase velocity" );

  S rad = 0.1;
  my_CLP.setOption( "radius", &rad, "radius of disk or sig of dipol" );

  S rot = 1.;
  my_CLP.setOption( "rotation", &rot, "rotation of disc" );

  S xm = 0.25;
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
  S l1 = 1.;
  my_CLP.setOption( "lx", &l1, "length in x-direction" );

  S l2 = 1.;
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

  O nf = 4.;
  my_CLP.setOption( "nf", &nf, "amount of grid points in f-direction" );

  O nfs = 1.;
  my_CLP.setOption( "nfs", &nfs, "start amount of grid points in f-direction" );

  O nfe = 1.;
  my_CLP.setOption( "nfe", &nfe, "end amount of grid points in f-direction" );

  // processor grid size
  O np1 = 2;
  my_CLP.setOption( "npx", &np1, "amount of processors in x-direction" );

  O np2 = 2;
  my_CLP.setOption( "npy", &np2, "amount of processors in y-direction" );

  O np3 = 1.;
  my_CLP.setOption( "npz", &np3, "amount of processors in z-direction" );

  // solver stuff
  std::string linSolName = "GMRES";
  my_CLP.setOption( "linSolName", &linSolName, "name of the linear solver" );

  std::string nonLinSolName = "Newton";
  my_CLP.setOption( "nonLinSolName", &nonLinSolName , "name of the non linear solver" );

  std::string lineSearchName = "Backtrack";
  my_CLP.setOption( "linesearch", &lineSearchName, "name of the line search" );

  int fixType = 1.;
  my_CLP.setOption( "fixType", &fixType, "type of fixpoint iteration matrix: 1=Newton, 2=Piccard, 3=lin diag... " );

  int precType = 1.;
  my_CLP.setOption( "precType", &precType, "type of preconditioners " );

  int maxIter = 10;
  my_CLP.setOption( "maxIter", &maxIter, "maximum iterations" );

  S tolBelos = 1.e-4;
  my_CLP.setOption( "tolBelos", &tolBelos, "tolerance for linear solver" );

  S tolNOX = 1.e-2;
  my_CLP.setOption( "tolNOX", &tolNOX, "tolerance for non-linear solver" );

  S tolNF = 1.e-2;
  my_CLP.setOption( "tolNF", &tolNF, "tolerance for non-linear solver" );


  my_CLP.recogniseAllOptions(true);
  my_CLP.throwExceptions(true);

  my_CLP.parse(argi,argv);
  // end of parsing

  if( nfs==nfe ) {
    nfs=nf;
    nfe=nf+1;
  }
  else {
    nf=nfe;
  }
  if( 7==flow ) {
    if( nfs<2 ) {
      nfs = 2 ;
      if( nfe <= nfs ) {
        ++nfe;
      }
    }
  }
  // starting with ininializing
  int rank = Pimpact::init_impact_pre();

  // outputs
  Teuchos::RCP<std::ostream> outPar;
  Teuchos::RCP<std::ostream> outLinSolve;
  Teuchos::RCP<std::ostream> outPrec;
  //  Teuchos::RCP<std::ostream> outLap2;
  //  Teuchos::RCP<std::ostream> outSchur;

  if(rank==0) {
    outPar   = Teuchos::rcp( new std::ofstream("para_case.txt") );
    outLinSolve  = Teuchos::rcp( new std::ofstream("stats_linSolve.txt") );
    outPrec  = Teuchos::rcp( new std::ofstream("stats_solvPrec.txt") );
    //    outSchur = Teuchos::rcp( new std::ofstream("stats_solvSchur.txt") );
  } else
    //    outPar = Teuchos::rcp( &blackhole, false) ;
    outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;

  *outPar << " \tflow=" << flow << "\n";
  *outPar << " \tforce=" << forcing << "\n";
  *outPar << " \tdomain=" << domain << "\n";
  *outPar << " \tre=" << re << "\n";
  *outPar << " \talpha2=" << alpha2 << "\n";

  auto ds = Pimpact::createDomainSize<S>( l1, l2, l3 );
  ds->set_Impact();
  ds->print( *outPar );

  //  auto bc = Pimpact::createBC( Pimpact::AllPeriodic );
  auto bc = Pimpact::createBoudaryConditionsGlobal( Pimpact::EDomainType(domain) );
  bc->set_Impact();

  auto gs = Pimpact::createGridSize<O>( n1, n2, n3 );
  gs->set_Impact();
  gs->print( *outPar );
  *outPar << " \tnf=" << nf << "\n";

  auto pgs = Pimpact::createProcGridSize<O>( np1, np2, np3 );
  pgs->set_Impact();
  pgs->print( *outPar );

  if(rank==0) {
    Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
  }
  outPar = Teuchos::null;

  // init IMPACT
  Pimpact::init_impact_post();


  // init Spaces
  auto fS = Pimpact::createFieldSpace<O>();
//  fS->print();

  auto iIS = Pimpact::createInnerFieldIndexSpaces<O>();
  auto fIS = Pimpact::createFullFieldIndexSpaces<O>();

  auto space = Pimpact::createSpace( fS, iIS, fIS );

  // init vectors
  auto x    = Pimpact::createMultiField( Pimpact::createCompoundField(
      Pimpact::createMultiHarmonicVectorField<S,O>( fS, iIS, fIS, nfs ),
      Pimpact::createMultiHarmonicScalarField<S,O>( fS, nfs )) );

  auto fu   = x->clone(Pimpact::ShallowCopy);
  //  fu->init( 0. );

  Teuchos::RCP<Pimpact::VectorField<S,O> > force=Teuchos::null;
  Teuchos::RCP<Pimpact::VectorField<S,O> > forcem1=Teuchos::null;
  if( 0!=forcing ) {
    force = x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(Pimpact::ShallowCopy);
    switch( Pimpact::EForceType(forcing) ) {
    case Pimpact::Dipol:
      force->initField( Pimpact::VPoint2D, rad );
      fu->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::VPoint2D, rad );
      break;
    case Pimpact::Disc:
      force->initField( Pimpact::Disc2D, xm*l1, ym*l2, rad );
      break;
    case Pimpact::RotatingDisc:
      force->initField( Pimpact::Disc2D, xm*l1, ym*l2, rad );
      fu->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::RotationDisc2D, xm*l1, ym*l2, rot );
      fu->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->scale( *force );
      break;
    case Pimpact::PseudoOscilatingDisc:
      force->initField( Pimpact::Disc2D, xm*l1, ym*l2, rad );
      fu->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->init( Teuchos::tuple(0., rot*rad, 0.) );
      fu->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->scale( *force );
      break;
    }
    forcem1 = force->clone(Pimpact::ShallowCopy);
    forcem1->init( 1. );
    forcem1->add( 1., *forcem1, -1., *force );
    //    force  ->write( 111 );
    //    forcem1->write( 222 );
  }

  // init Fields, init and rhs
  switch( Pimpact::EFlowType(flow) ) {
  case Pimpact::Zero2DFlow:
    break;
  case Pimpact::Poiseuille_inX:
    x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->initField( Pimpact::Poiseuille2D_inX );
    break;
  case Pimpact::Poiseuille_inY:
    x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->initField( Pimpact::Poiseuille2D_inY );
    break;
  case Pimpact::Pulsatile_inX:
    x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::Pulsatile2D_inXC, -1. );
    x->getFieldPtr(0)->getVFieldPtr()->getSFieldPtr(0)->initField( Pimpact::Pulsatile2D_inXS,  1. );
    break;
  case Pimpact::Pulsatile_inY:
    x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::Pulsatile2D_inYC, -1. );
    x->getFieldPtr(0)->getVFieldPtr()->getSFieldPtr(0)->initField( Pimpact::Pulsatile2D_inYS,  1. );
    break;
  case Pimpact::Streaming2DFlow:
    x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::Streaming2DS, 1. );
    break;
  case Pimpact::Streaming2DFlow2:
    x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::Streaming2DS, -1. );
    x->getFieldPtr(0)->getVFieldPtr()->getSFieldPtr(0)->initField( Pimpact::Streaming2DC,  1. );
    break;
  case Pimpact::Streaming2DFlow3:
    x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr( )->initField( Pimpact::Streaming2DS, -0.5 );
    x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(1)->initField( Pimpact::Streaming2DS, -0.5 );
    x->getFieldPtr(0)->getVFieldPtr()->getSFieldPtr(1)->initField( Pimpact::Streaming2DC,  0.5 );
    break;
  default:
    x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::Streaming2DS, 1. );
    break;
  }

    x->init( 0. );
  //  x->random();



  /******************************************************************************************/
  //  for( nf=nfs; nf<nfe; ++nf) {
  //  for( nf=nfs; nf<nfe; nf+=2 ) {
  for( nf=nfs; nf<nfe; nf*=2) {

    if( nf!=nfs ) {
      S toltemp = x->getConstFieldPtr(0)->getConstVFieldPtr()->getConstFieldPtr(x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes()-1)->norm()/std::sqrt(l1*l2/n1/n2);
      if(0==rank) std::cout << "\n\t--- ||u_Nf||: "<<toltemp<<"\t---\n";
      if( toltemp < tolNF ) {
        if(0==rank) std::cout << "\n\t--- Nf: "<<x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes()<<"\tdof: "<<x->getLength(true)<<"\t---\n";
        break;
      }
      do {
        x->getFieldPtr(0)->getVFieldPtr()->push_back();
        x->getFieldPtr(0)->getSFieldPtr()->push_back();
        fu->getFieldPtr(0)->getVFieldPtr()->push_back();
        fu->getFieldPtr(0)->getSFieldPtr()->push_back();
      } while( x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes() < nf );
            tolNOX /= 5;
//      tolNOX /= 10;
    }

    //    if(0==rank) std::cout << "\n\t--- Nf: "<<nf<<"\tdof: "<<x->getLength(true)<<"\t---\n";
    if(0==rank) std::cout << "\n\t--- Nf: "<<x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes()<<"\tdof: "<<x->getLength(true)<<"\t---\n";

    auto para = Pimpact::createLinSolverParameter( linSolName, tolBelos*l1*l2/n1/n2*(nfe-1)/nf, -1 );
    //    auto para = Pimpact::createLinSolverParameter( linSolName, tol, -1 );
    para->set( "Maximum Iterations", 3000 );
    para->set( "Implicit Residual Scaling", "Norm of RHS" );
    para->set( "Explicit Residual Scaling", "Norm of RHS" );
    //    para->set( "Output Stream", outLinSolve );

    auto dtl = Pimpact::createMultiDtHelmholtz<S,O>( alpha2/re, 1./re );

    Teuchos::RCP<Fo> forcingOp = Teuchos::null;
    Teuchos::RCP<Fo> forcingm1Op = Teuchos::null;
    if( 0!=forcing ) {
      forcingOp   = Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force   ) );
      forcingm1Op = Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( forcem1 ) );
    }


    auto opV2V =
        Pimpact::createAdd3Op(
            x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
            Pimpact::createCompositionOp(
                x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
                forcingm1Op,
                Pimpact::createAdd3Op(
                    x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
                    dtl,
                    Pimpact::createMultiHarmonicNonlinear<S,O>( /*x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone()*/ )
                )
            ),
            forcingOp
        );

    //    auto opS2V = Pimpact::createMultiHarmonicOpWrap< Pimpact::Grad<S,O> >();
    auto opS2V =
        Pimpact::createCompositionOp(
            x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
            forcingm1Op,
            Pimpact::createMultiHarmonicOpWrap(
                Pimpact::createGradOp<S,O>()
            )
        );
    auto opV2S = Pimpact::createMultiHarmonicOpWrap< Pimpact::Div<S,O> >();

    auto op =
        Pimpact::createMultiOperatorBase<MF>(
            Pimpact::createCompoundOpWrap(
                x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
                opV2V,
                opS2V,
                opV2S )
        );


    Teuchos::RCP<BOp> jop;

    /// this could be a factory in goes opS2V, opV2S, forcingOp, forcing1mOp, x,
    if( 1==fixType ) {
      if(0==rank) std::cout << "\n\t---\titeration matrix(1): full Newton iteration\t---\n";

      auto opV2V = Pimpact::createAdd3Op(
          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
          Pimpact::createCompositionOp(
              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
              forcingm1Op,
              Pimpact::createAdd3Op(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  dtl,
                  Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
                      //                    x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(Pimpact::ShallowCopy),
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy)  ) ) ),
                      forcingOp );
      jop =
          Pimpact::createMultiOperatorBase<MF>(
              Pimpact::createCompoundOpWrap(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  opV2V,
                  opS2V,
                  opV2S ) );
    }
    else if( 2==fixType ){

      if(0==rank) std::cout << "\n\t---\titeration matrix(2): full Picard iteration\t---\n";
      jop =
          Pimpact::createMultiOperatorBase<MF>(
              Pimpact::createCompoundOpWrap(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  Pimpact::createAdd3Op(
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      Pimpact::createCompositionOp(
                          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                          forcingm1Op,
                          Pimpact::createAdd3Op(
                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                              dtl,
                              Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
                                  //                                x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(Pimpact::ShallowCopy),
                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                                  false  ) ) ),
                                  forcingOp
                  ),
                  opS2V,
                  opV2S ) );
    }
    else if( 3==fixType ){

      if(0==rank) std::cout << "\n\t---\titeration matrix(3): complex diagonal Newton iteration\t---\n";
      jop =
          Pimpact::createMultiOperatorBase<MF>(
              Pimpact::createCompoundOpWrap(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  Pimpact::createAdd3Op(
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      Pimpact::createCompositionOp(
                          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                          forcingm1Op,
                          Pimpact::createAdd3Op(
                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                              dtl,
                              Pimpact::createMultiHarmonicDiagNonlinearJacobian<S,O>(
                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy), true )
                          )
                      ),
                      forcingOp
                  ),
                  opS2V,
                  opV2S ) );
    }
    else if( 4==fixType ){

      if(0==rank) std::cout << "\n\t---\titeration matrix(4): complex diagonal Picard iteration\t---\n";
      jop =
          Pimpact::createMultiOperatorBase<MF>(
              Pimpact::createCompoundOpWrap(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  Pimpact::createAdd3Op(
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      Pimpact::createCompositionOp(
                          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                          forcingm1Op,
                          Pimpact::createAdd3Op(
                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                              dtl,
                              Pimpact::createMultiHarmonicDiagNonlinearJacobian<S,O>(
                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy), false )
                          )
                      ),
                      forcingOp
                  ),
                  opS2V,
                  opV2S ) );

    }
    else if( 5==fixType ){

      if(0==rank) std::cout << "\n\t---\titeration matrix(5): semi-complex diagonal Newton iteration\t---\n";
      jop =
          Pimpact::createMultiOperatorBase<MF>(
              Pimpact::createCompoundOpWrap(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  Pimpact::createAdd3Op(
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      Pimpact::createCompositionOp(
                          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                          forcingm1Op,
                          Pimpact::createAdd3Op(
                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                              dtl,
                              Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
                                  Pimpact::createNonlinearJacobian<S,O>(
                                      x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(Pimpact::ShallowCopy), true ) )
                          )
                      ),
                      forcingOp
                  ),
                  opS2V,
                  opV2S ) );

    }
    else if( 6==fixType){

      if(0==rank) std::cout << "\n\t---\titeration matrix(6): semi-complex diagonal Picard iteration\t---\n";
      jop =
          Pimpact::createMultiOperatorBase<MF>(
              Pimpact::createCompoundOpWrap(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  Pimpact::createAdd3Op(
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      Pimpact::createCompositionOp(
                          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                          forcingm1Op,
                          Pimpact::createAdd3Op(
                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                              dtl,
                              Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
                                  Pimpact::createNonlinearJacobian<S,O>(
                                      x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(Pimpact::ShallowCopy), false ) )
                          )
                      ),
                      forcingOp
                  ),
                  opS2V,
                  opV2S ) );

    }
    else if( 9==fixType ) {

      if(0==rank) std::cout << "\n\t---\titeration matrix(9): linear terms\t---\n";
      jop =
          Pimpact::createMultiOperatorBase<MF>(
              Pimpact::createCompoundOpWrap(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  Pimpact::createAdd3Op(
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      Pimpact::createCompositionOp(
                          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                          forcingm1Op,
                          dtl
                      ),
                      forcingOp
                  ),
                  opS2V,
                  opV2S ) );
    }
    else if( 10==fixType ) {

      if(0==rank) std::cout << "\n\t---\titeration matrix(10): full Newton iteration Schur complement\t---\n";

      auto opV2V =
          Pimpact::createMultiOperatorBase<MVF>(
              Pimpact::createAdd3Op(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  Pimpact::createCompositionOp(
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      forcingm1Op,
                      Pimpact::createAdd3Op(
                          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                          dtl ,//),
                          Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy), false ) ) ),
                              forcingOp ) );
      auto lp_ = Pimpact::createLinearProblem<MVF>(
          opV2V, Teuchos::null, Teuchos::null, para, linSolName );
      auto opV2Vinv = Pimpact::createInverseOperatorBase<MVF>( lp_ );

      auto invSchur = Pimpact::createInverseSchurOp(
          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
          x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
          opV2Vinv,
          opS2V,
          opV2S);

      jop =
          Pimpact::createMultiOperatorBase<MF>(
              invSchur );
      //      jop =
      //          Pimpact::createMultiOperatorBase<MF>(
      //              Pimpact::createTripleCompositionOp(
      //                   x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
      //                   x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
      //                   op
      //                   opV2Vinv,


      //      jop =
      //          Pimpact::createMultiOperatorBase<MF>(
      //              Pimpact::createCompoundOpWrap(
      //                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(),
      //                  opV2V,
      //                  opS2V,
      //                  opV2S ) );

    }


    // init lprec
    Teuchos::RCP<BOp> lprec = Teuchos::null;
    para->set( "Output Stream", outPrec );

    if( 10==precType ) {
      if(0==rank) std::cout << "\n\t---\tprec Type(10): full Newton iteration Schur complement\t---\n";
      auto opV2V =
          Pimpact::createMultiOperatorBase<MVF>(
              Pimpact::createAdd3Op(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  Pimpact::createCompositionOp(
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      forcingm1Op,
                      //                          Pimpact::createAdd3Op(
                      //                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      dtl ),//,
                      //                              Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
                      //                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(Pimpact::ShallowCopy),
                      //                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy), false ) ) ),
                      forcingOp ) );
      auto lp_ = Pimpact::createLinearProblem<MVF>(
          opV2V, Teuchos::null, Teuchos::null, Teuchos::parameterList(), linSolName );
      auto opV2Vinv = Pimpact::createInverseOperatorBase<MVF>( lp_ );

      auto invSchur = Pimpact::createInverseSchurOp(
          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
          x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
          opV2Vinv,
          opS2V,
          opV2S);

      lprec =
          Pimpact::createMultiOperatorBase<MF>(
              invSchur );
    }
    else if( 11==precType ) {
      if(0==rank) std::cout << "\n\t---\tprec Type(10): full Newton iteration Schur complement\t---\n";

      auto prec = Pimpact::createMultiHarmonicMLEddy<S,O>( space, 20, nf, alpha2, 1./re );

      auto opV2V =
          Pimpact::createMultiOperatorBase<MVF>(
              Pimpact::createAdd3Op(
                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                  Pimpact::createCompositionOp(
                      x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      forcingm1Op,
                      //                          Pimpact::createAdd3Op(
                      //                              x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
                      dtl ),//,
                      //                              Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
                      //                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->getConst0FieldPtr()->clone(Pimpact::ShallowCopy),
                      //                                  x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy), false ) ) ),
                      forcingOp ) );

      //     auto lp_ = Pimpact::createLinearProblem<MVF>(
      //           opV2V, Teuchos::null, Teuchos::null, Teuchos::parameterList(), linSolName );
      //     lp_->setLeftPrec(prec);
      //     auto opV2Vinv = Pimpact::createInverseOperatorBase<MVF>( lp_ );

      auto invSchur = Pimpact::createInverseSchurOp(
          x->getConstFieldPtr(0)->getConstVFieldPtr()->clone(Pimpact::ShallowCopy),
          x->getConstFieldPtr(0)->getConstSFieldPtr()->clone(Pimpact::ShallowCopy),
          prec,
          //         opV2Vinv,
          opS2V,
          opV2S);

      lprec =
          Pimpact::createMultiOperatorBase<MF>(
              invSchur );
    }
    //  if( 3==precType ){
    //
    //    if(0==rank) std::cout << "\n\t---\tpreconditioner matrix(3): complex diagonal Newton iteration\t---\n";
    //
    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<DJMAdv,DtL>, Fo > > JOp;
    //
    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
    //        Pimpact::createMultiOpWrap(
    //            Pimpact::createAdd2Op<Pimpact::AddOp<DJMAdv,DtL>,Fo>(
    //                Pimpact::createAdd2Op<DJMAdv,DtL>(
    //                    Pimpact::createMultiHarmonicDiagNonlinearJacobian<S,O>(
    //                        temp->getFieldPtr(0), true ),
    //                    Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
    //                    temp->getFieldPtr(0)->clone() ) ,
    //                Pimpact::createMultiHarmonicOpWrap(
    //                    Pimpact::createForcingOp<S,O>( force ) ) ,
    //                temp->getConstFieldPtr(0)->clone() ) )
    //    );
    //    auto prob =
    //          Pimpact::createLinearProblem<MVF>(
    //              precOp,
    //              fu,
    //              fu, para, "GMRES" );
    //
    //      lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
    ////      break;
    //  }
    //  else if( 4==precType ){
    //
    //    if(0==rank) std::cout << "\n\t---\tpreconditioner  matrix(4): complex diagonal Picard iteration\t---\n";
    //
    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<DJMAdv,DtL>, Fo > > JOp;
    //
    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
    //        Pimpact::createMultiOpWrap(
    //            Pimpact::createAdd2Op<Pimpact::AddOp<DJMAdv,DtL>,Fo>(
    //                Pimpact::createAdd2Op<DJMAdv,DtL>(
    //                    Pimpact::createMultiHarmonicDiagNonlinearJacobian<S,O>(
    //                        temp->getFieldPtr(0), false ),
    //                    Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
    //                    temp->getFieldPtr(0)->clone() ) ,
    //                Pimpact::createMultiHarmonicOpWrap(
    //                    Pimpact::createForcingOp<S,O>( force ) ) ,
    //                temp->getConstFieldPtr(0)->clone() ) )
    //    );
    //    auto prob =
    //          Pimpact::createLinearProblem<MVF>(
    //              precOp,
    //              fu,
    //              fu, para, "GMRES" );
    //
    //      lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
    //  }
    //  else if( 5==precType ){
    //
    //    if(0==rank) std::cout << "\n\t---\titeration matrix(5): semi-complex diagonal Newton iteration\t---\n";
    //
    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<JAdv,DtL>, Fo > > JOp;
    //
    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
    //        Pimpact::createMultiOpWrap(
    //            Pimpact::createAdd2Op<Pimpact::AddOp<JAdv,DtL>,Fo>(
    //                Pimpact::createAdd2Op<JAdv,DtL>(
    //                    Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
    //                        Pimpact::createNonlinearJacobian<S,O>( temp->getFieldPtr(0)->get0FieldPtr(), true ) ),
    //                    Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
    //                    temp->getFieldPtr(0)->clone() ) ,
    //                Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
    //                temp->getConstFieldPtr(0)->clone() ) )
    //    );
    //    auto prob =
    //          Pimpact::createLinearProblem<MVF>(
    //              precOp,
    //              fu,
    //              fu, para, "GMRES" );
    //
    //      lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
    //  }
    //  else if( 6==precType ){
    //
    //    if(0==rank) std::cout << "\n\t---\titeration matrix(6): semi-complex diagonal Picard iteration\t---\n";
    //
    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<JAdv,DtL>, Fo > > JOp;
    //
    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
    //        Pimpact::createMultiOpWrap(
    //            Pimpact::createAdd2Op<Pimpact::AddOp<JAdv,DtL>,Fo>(
    //                Pimpact::createAdd2Op<JAdv,DtL>(
    //                    Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
    //                        Pimpact::createNonlinearJacobian<S,O>( temp->getFieldPtr(0)->get0FieldPtr(), false ) ),
    //                    Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
    //                    temp->getFieldPtr(0)->clone() ) ,
    //                Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
    //                temp->getConstFieldPtr(0)->clone() ) )
    //    );
    //    auto prob =
    //          Pimpact::createLinearProblem<MVF>(
    //              precOp,
    //              fu,
    //              fu, para, "GMRES" );
    //    lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
    //  }
    //  else if( 7==precType ){
    //
    //    if(0==rank) std::cout << "\n\t---\titeration matrix(7): real diagonal Newton iteration\t---\n";
    //
    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<JAdv,DtL>, Fo > > JOp;
    //
    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
    //        Pimpact::createMultiOpWrap(
    //            Pimpact::createAdd2Op<Pimpact::AddOp<JAdv,DtL>,Fo>(
    //                Pimpact::createAdd2Op<JAdv,DtL>(
    //                    Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
    //                        Pimpact::createNonlinearJacobian<S,O>( temp->getFieldPtr(0)->get0FieldPtr(), true) ),
    //                    Pimpact::createMultiDtHelmholtz<S,O>( 0., 0., 1./re ),
    //                    temp->getFieldPtr(0)->clone() ) ,
    //                Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
    //                temp->getConstFieldPtr(0)->clone() ) )
    //    );
    //    auto prob =
    //          Pimpact::createLinearProblem<MVF>(
    //              precOp,
    //              fu,
    //              fu, para, "GMRES" );
    //    lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
    //  }
    //  else if( 8==precType ){
    //    if(0==rank) std::cout << "\n\t---\titeration matrix(8): real diagonal Picard iteration\t---\n";
    //
    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< Pimpact::AddOp<JAdv,DtL>, Fo > > JOp;
    //
    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
    //        Pimpact::createMultiOpWrap(
    //            Pimpact::createAdd2Op<Pimpact::AddOp<JAdv,DtL>,Fo>(
    //                Pimpact::createAdd2Op<JAdv,DtL>(
    //                    Pimpact::createMultiHarmonicOpWrap<Pimpact::NonlinearJacobian<S,O> >(
    //                        Pimpact::createNonlinearJacobian<S,O>( temp->getFieldPtr(0)->get0FieldPtr(), false) ),
    //                    Pimpact::createMultiDtHelmholtz<S,O>( 0., 0., 1./re ),
    //                    temp->getFieldPtr(0)->clone() ) ,
    //                Pimpact::createMultiHarmonicOpWrap( Pimpact::createForcingOp<S,O>( force ) ) ,
    //                temp->getConstFieldPtr(0)->clone() ) )
    //    );
    //    auto prob =
    //          Pimpact::createLinearProblem<MVF>(
    //              precOp,
    //              fu,
    //              fu, para, "GMRES" );
    //    lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
    //  }
    //  else if( 9==precType ) {
    //    if(0==rank) std::cout << "\n\t---\titeration matrix(9): linear terms\t---\n";
    //    typedef Pimpact::MultiOpWrap< Pimpact::AddOp< DtL, Fo > > JOp;
    //
    //    auto precOp = Pimpact::createOperatorBase<MVF,JOp>(
    //           Pimpact::createMultiOpWrap(
    //               Pimpact::createAdd2Op<DtL,Fo>(
    ////                   Pimpact::createAdd2Op<JMAdv,DtL>(
    ////                       Pimpact::createMultiHarmonicNonlinearJacobian<S,O>(
    ////                           x->getConstFieldPtr(0)->getConst0FieldPtr()->clone(), x->getConstFieldPtr(0)->clone(), false ),
    //                   Pimpact::createMultiDtHelmholtz<S,O>( alpha2, 0., 1./re ),
    ////                           temp->getFieldPtr(0)->clone()  ,
    //                   Pimpact::createMultiHarmonicOpWrap(
    //                       Pimpact::createForcingOp<S,O>( force ) ) ,
    //                   temp->getFieldPtr(0)->clone() ) )
    //       );
    //    auto prob =
    //          Pimpact::createLinearProblem<MVF>(
    //              precOp,
    //              fu,
    //              fu, para, "GMRES" );
    //    lprec = Pimpact::createOperatorBase< MVF, Pimpact::InverseOperator<MVF> >( Pimpact::createInverseOperator<MVF>( prob ) );
    //  }
    //  else {
    //    lprec = Teuchos::null;
    //  }


    if( 10!=fixType ) {
      para->set( "Output Stream", outLinSolve );
      auto lp_ = Pimpact::createLinearProblem<MF>(
          jop, x->clone(), fu->clone(), para, linSolName );
      lp_->setLeftPrec( lprec );
      auto lp = Pimpact::createInverseOperatorBase<MF>( lp_ );

      jop=lp;
    }


    auto inter = NOX::Pimpact::createInterface<MF>( fu, op, jop );

    auto nx = NOX::Pimpact::createVector(x);

    auto bla = Teuchos::parameterList();

    auto group = NOX::Pimpact::createGroup<Inter>( bla, inter, nx );

    // Set up the status tests
    auto statusTest = NOX::Pimpact::createStatusTest( maxIter, tolNOX, tolBelos );

    // Create the list of solver parameters
    auto solverParametersPtr =
        NOX::Pimpact::createNOXSolverParameter( nonLinSolName, lineSearchName );

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
        NOX::Solver::buildSolver( group, statusTest, solverParametersPtr);

    // Solve the nonlinear system
    //    NOX::StatusTest::StatusType status =
    {
      Teuchos::TimeMonitor LocalTimer(*CompTime);
      solver->solve();
    }

    // Print the parameter list
    //    if( nf==nfs ) {
    //      if(rank==0) std::cout << "\n" << "-- Parameter List From Solver --" << "\n";
    //      if(rank==0) std::cout << "\n" << status << "\n";
    //      if(rank==0) solver->getList().print(std::cout);
    //    }

    // Get the answer
    *group = solver->getSolutionGroup();

    // Print the answer if(rank==0) std::cout << "\n" << "-- Final Solution From Solver --" << "\n";

    //    Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() )->getConstFieldPtr()->write(800);
    //    Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(900);

    x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr();
    //    x = Teuchos::rcp_const_cast<NV>( group->getXPtr() )->getFieldPtr();

//    Teuchos::TimeMonitor::summarize();
        x->write(nf*100);
//      x->write();
    //  x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() ))->getFieldPtr();
    //  x->write(900);
  }
  /******************************************************************************************/

  // Get a summary from the time monitor.
  Teuchos::TimeMonitor::summarize();

  //  write so
//  x->write();


  MPI_Finalize();
  return( 0 );

}
