#include "Pimpact_Operator.hpp"




#ifdef COMPILE_ETI
namespace Pimpact {
	// BaseOperators
	template class DivOp< Space<double,int,3,2> >;
	template class DivOp< Space<double,int,3,4> >;
	template class DivOp< Space<double,int,4,2> >;
	template class DivOp< Space<double,int,4,4> >;
	template class GradOp< Space<double,int,3,2> >;
	template class GradOp< Space<double,int,3,4> >;
	template class GradOp< Space<double,int,4,2> >;
	template class GradOp< Space<double,int,4,4> >;
	template class DivGradOp< Space<double,int,3,2> >;
	template class DivGradOp< Space<double,int,3,4> >;
	template class DivGradOp< Space<double,int,4,2> >;
	template class DivGradOp< Space<double,int,4,4> >;
	template class DivGradO2Op< Space<double,int,3,2> >;
	template class DivGradO2Op< Space<double,int,3,4> >;
	template class DivGradO2Op< Space<double,int,4,2> >;
	template class DivGradO2Op< Space<double,int,4,4> >;
	template class DivGradO2JSmoother< DivGradO2Op<Space<double,int,3,2> > >;
	template class DivGradO2JSmoother< DivGradO2Op<Space<double,int,3,4> > >;
	template class DivGradO2JSmoother< DivGradO2Op<Space<double,int,4,2> > >;
	template class DivGradO2JSmoother< DivGradO2Op<Space<double,int,4,4> > >;
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
template class Pimpact::TimeStokesBSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,2> > >;
template class Pimpact::TimeStokesBSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,4> > >;

#endif
