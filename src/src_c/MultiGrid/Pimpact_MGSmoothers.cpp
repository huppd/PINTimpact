#include "Pimpact_MGSmoothers.hpp"

#ifdef COMPILE_ETI
template class Pimpact::MGSmoothers< Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >, Pimpact::DivGradOp, Pimpact::DivGradO2Op >, Pimpact::DivGradO2JSmoother >;
template class Pimpact::MGSmoothers< Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >, Pimpact::DivGradOp, Pimpact::DivGradO2Op >, Pimpact::DivGradO2JSmoother >;
#endif
