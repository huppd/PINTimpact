#include "Pimpact_MGOperators.hpp"

#ifdef COMPILE_ETI
template class Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >, Pimpact::DivGradOp, Pimpact::DivGradOp >;
template class Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >, Pimpact::DivGradOp, Pimpact::DivGradOp >;
template class Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >, Pimpact::DivGradOp, Pimpact::DivGradO2Op >;
template class Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >, Pimpact::DivGradOp, Pimpact::DivGradO2Op >;
template class Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >, ConvDiffOpT,ConvDiffOpT >;
template class Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >, ConvDiffOpT,ConvDiffOpT >;
#endif
