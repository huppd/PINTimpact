#include "Pimpact_Operator.hpp"




#ifdef COMPILE_ETI

// BaseOperators
template class Pimpact::DivOp< Pimpact::Space<double,int,3,2> >;
template class Pimpact::DivOp< Pimpact::Space<double,int,3,4> >;
template class Pimpact::DivOp< Pimpact::Space<double,int,4,2> >;
template class Pimpact::DivOp< Pimpact::Space<double,int,4,4> >;
template class Pimpact::GradOp< Pimpact::Space<double,int,3,2> >;
template class Pimpact::GradOp< Pimpact::Space<double,int,3,4> >;
template class Pimpact::GradOp< Pimpact::Space<double,int,4,2> >;
template class Pimpact::GradOp< Pimpact::Space<double,int,4,4> >;
template class Pimpact::DivGradOp< Pimpact::Space<double,int,3,2> >;
template class Pimpact::DivGradOp< Pimpact::Space<double,int,3,4> >;
template class Pimpact::DivGradOp< Pimpact::Space<double,int,4,2> >;
template class Pimpact::DivGradOp< Pimpact::Space<double,int,4,4> >;
template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,2> >;
template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,4> >;
template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,2> >;
template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,4> >;
template class Pimpact::DivGradO2Inv< Pimpact::DivGradO2Op<Pimpact::Space<double,int,3,2> > >;
template class Pimpact::DivGradO2Inv< Pimpact::DivGradO2Op<Pimpact::Space<double,int,3,4> > >;
template class Pimpact::DivGradO2Inv< Pimpact::DivGradO2Op<Pimpact::Space<double,int,4,2> > >;
template class Pimpact::DivGradO2Inv< Pimpact::DivGradO2Op<Pimpact::Space<double,int,4,4> > >;
template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op<Pimpact::Space<double,int,3,2> > >;
template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op<Pimpact::Space<double,int,3,4> > >;
template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op<Pimpact::Space<double,int,4,2> > >;
template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op<Pimpact::Space<double,int,4,4> > >;
namespace Pimpact {
template class DivGradO2SORSmoother< DivGradO2Op<Space<double,int,3,2> > >;
template class DivGradO2SORSmoother< DivGradO2Op<Space<double,int,3,4> > >;
template class DivGradO2SORSmoother< DivGradO2Op<Space<double,int,4,2> > >;
template class DivGradO2SORSmoother< DivGradO2Op<Space<double,int,4,4> > >;
template class HelmholtzOp< Space<double,int,3,2> >;
template class HelmholtzOp< Space<double,int,3,4> >;
template class HelmholtzOp< Space<double,int,4,2> >;
template class HelmholtzOp< Space<double,int,4,4> >;
template class InterpolateS2V< Space<double,int,3,2> >;
template class InterpolateS2V< Space<double,int,3,4> >;
template class InterpolateS2V< Space<double,int,4,2> >;
template class InterpolateS2V< Space<double,int,4,4> >;
template class InterpolateV2S<double,int,3,2>;
template class InterpolateV2S<double,int,3,4>;
template class InterpolateV2S<double,int,4,2>;
template class InterpolateV2S<double,int,4,4>;
//	template class TransferOp< Space<double,int,3,4>, Space<double,int,3,2> >;
//	template class TransferOp< Space<double,int,4,4>, Space<double,int,4,2> >;
template class ConvectionSOp< Space<double,int,3,2> >;
template class ConvectionSOp< Space<double,int,3,4> >;
template class ConvectionSOp< Space<double,int,4,2> >;
template class ConvectionSOp< Space<double,int,4,4> >;
template class ConvectionDiffusionSOp< Space<double,int,3,2> >;
template class ConvectionDiffusionSOp< Space<double,int,3,4> >;
template class ConvectionDiffusionSOp< Space<double,int,4,2> >;
template class ConvectionDiffusionSOp< Space<double,int,4,4> >;
template class ConvectionDiffusionJSmoother< ConvectionDiffusionSOp< Space<double,int,3,2> > >;
template class ConvectionDiffusionJSmoother< ConvectionDiffusionSOp< Space<double,int,3,4> > >;
template class ConvectionDiffusionJSmoother< ConvectionDiffusionSOp< Space<double,int,4,2> > >;
template class ConvectionDiffusionJSmoother< ConvectionDiffusionSOp< Space<double,int,4,4> > >;
template class ConvectionDiffusionSORSmoother< ConvectionDiffusionSOp< Space<double,int,3,2> > >;
template class ConvectionDiffusionSORSmoother< ConvectionDiffusionSOp< Space<double,int,3,4> > >;
template class ConvectionDiffusionSORSmoother< ConvectionDiffusionSOp< Space<double,int,4,2> > >;
template class ConvectionDiffusionSORSmoother< ConvectionDiffusionSOp< Space<double,int,4,4> > >;
template class NonlinearOp< ConvectionDiffusionSOp< Space<double,int,3,2> > >;
template class NonlinearOp< ConvectionDiffusionSOp< Space<double,int,3,4> > >;
template class NonlinearOp< ConvectionDiffusionSOp< Space<double,int,4,2> > >;
template class NonlinearOp< ConvectionDiffusionSOp< Space<double,int,4,4> > >;
template class NonlinearOp< ConvectionSOp< Space<double,int,3,2> > >;
template class NonlinearOp< ConvectionSOp< Space<double,int,3,4> > >;
template class NonlinearOp< ConvectionSOp< Space<double,int,4,2> > >;
template class NonlinearOp< ConvectionSOp< Space<double,int,4,4> > >;
template class NonlinearSmoother< NonlinearOp< ConvectionDiffusionSOp< Space<double,int,3,2> > >, ConvectionDiffusionSORSmoother>;
template class NonlinearSmoother< NonlinearOp< ConvectionDiffusionSOp< Space<double,int,3,4> > >, ConvectionDiffusionSORSmoother>;
template class NonlinearSmoother< NonlinearOp< ConvectionDiffusionSOp< Space<double,int,4,2> > >, ConvectionDiffusionSORSmoother>;
template class NonlinearSmoother< NonlinearOp< ConvectionDiffusionSOp< Space<double,int,4,4> > >, ConvectionDiffusionSORSmoother>;

} // end of namespace Pimpact

template class Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,2> >;
template class Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,4> >;
template class Pimpact::TimeNSOp< Pimpact::Space<double,int,4,2> >;
template class Pimpact::TimeNSOp< Pimpact::Space<double,int,4,4> >;
template class Pimpact::TimeStokesBSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,2> > >;
template class Pimpact::TimeStokesBSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,4> > >;
template class Pimpact::TimeNSBSmoother< Pimpact::TimeNSOp< Pimpact::Space<double,int,4,2> > >;
template class Pimpact::TimeNSBSmoother< Pimpact::TimeNSOp< Pimpact::Space<double,int,4,4> > >;
template class Pimpact::TimeNS4DBSmoother< Pimpact::TimeNSOp< Pimpact::Space<double,int,4,2> > >;
template class Pimpact::TimeNS4DBSmoother< Pimpact::TimeNSOp< Pimpact::Space<double,int,4,4> > >;
template class Pimpact::TimeStokesLSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,2> > >;
template class Pimpact::TimeStokesLSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,4> > >;
#endif
