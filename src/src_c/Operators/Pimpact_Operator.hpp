#pragma once
#ifndef PIMPACT_OPERATOR_HPP
#define PIMPACT_OPERATOR_HPP


#include "Pimpact_Fields.hpp" // ETI

/// \defgroup Operator Operators
#include "Pimpact_Add2Op.hpp"
#include "Pimpact_Add3Op.hpp"
#include "Pimpact_Chebyshev.hpp"
#include "Pimpact_CompositionOp.hpp"
#include "Pimpact_TripleCompositionOp.hpp"
#include "Pimpact_InverseOp.hpp"
#include "Pimpact_InvDiagonal.hpp"
#include "Pimpact_PrecInverseOp.hpp"
#include "Pimpact_ForcingOp.hpp"
#include "Pimpact_Transpose.hpp"


/// \defgroup BaseOperator Basic Operator
/// \ingroup Operator
/// Operator operate on basic fields, \c Pimpact::ScalarField and \c Pimpact::VectorField

/// \defgroup NonliearOperator  Nonlinear Operators
/// \ingroup BaseOperator
//// Basic Operator that are nonlinear, such that the can be wrapped by \c Pimpact::ConvectionVWrap
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_DivGradOp.hpp"
#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_DivGradO2Inv.hpp"
#include "Pimpact_DivGradO2JSmoother.hpp"
#include "Pimpact_DivGradO2LSmoother.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_InterpolateS2VOp.hpp"
#include "Pimpact_InterpolateV2SOp.hpp"
#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_ConvectionDiffusionJSmoother.hpp"
#include "Pimpact_ConvectionDiffusionSORSmoother.hpp"
#include "Pimpact_NonlinearOp.hpp"
#include "Pimpact_NonlinearSmoother.hpp"
#include "Pimpact_NonlinearVWrap.hpp"




/// \defgroup TimeOperator Time Operator
/// \ingroup Operator
/// Operator operate on time fields: \c Pimpact::TimeField
#include "Pimpact_DtTimeOp.hpp"
#include "Pimpact_TimeOpWrap.hpp"
#include "Pimpact_TimeDtConvectionDiffusionOp.hpp"


/// \defgroup ModeOperator Mode Operator
/// \ingroup Operator
/// Operator operate on mode fields: \c Pimpact::ModeField
#include "Pimpact_EddyPrec.hpp"
#include "Pimpact_ModeOpWrap.hpp"
#include "Pimpact_ModeNonlinearOp.hpp"

/// \defgroup MultiOperator Multi Operator
/// \ingroup Operator
/// Operator operate on multi fields: \c Pimpact::MultiField
#include "Pimpact_MultiOpWrap.hpp"
#include "Pimpact_MultiOpUnWrap.hpp"

/// \defgroup CompoundOperator Compound Operator
/// \ingroup Operator
/// Operator operate on compound fields: \c Pimpact::CompoundField
#include "Pimpact_CompoundOpWrap.hpp"
#include "Pimpact_InverseTriangularOp.hpp"
#include "Pimpact_TimeStokesOp.hpp"
#include "Pimpact_TimeStokesBSmoother.hpp"
#include "Pimpact_TimeStokesLSmoother.hpp"
#include "Pimpact_TimeNSOp.hpp"
#include "Pimpact_TimeNSBSmoother.hpp"
#include "Pimpact_TimeNS4DBSmoother.hpp"


/// \defgroup MultiHarmonicOperator Multi-harmonic Operator
/// \ingroup Operator
/// Operator operate on multi-harmonic fields: \c Pimpact::MultiHarmonicField
#include "Pimpact_MultiHarmonicOpWrap.hpp"
#include "Pimpact_MultiHarmonicMultiOpWrap.hpp"
#include "Pimpact_MultiHarmonicDiagOp.hpp"
#include "Pimpact_MultiDtConvectionDiffusionOp.hpp"

#include "Pimpact_OperatorFactory.hpp"

#endif // end of #ifndef PIMPACT_OPERATOR_HPP
