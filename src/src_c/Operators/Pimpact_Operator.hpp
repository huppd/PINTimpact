#pragma once
#ifndef PIMPACT_OPERATOR_HPP
#define PIMPACT_OPERATOR_HPP

#include "Pimpact_Fields.hpp" // ETI

/// \defgroup Operator Operators
/// \todo think about adding a compResidual( x, b, res ) method to Operators
#include "Pimpact_Add2Op.hpp"
#include "Pimpact_Add3Op.hpp"
#include "Pimpact_CompositionOp.hpp"
#include "Pimpact_TripleCompositionOp.hpp"
#include "Pimpact_InverseOperator.hpp"
#include "Pimpact_InverseOp.hpp"
#include "Pimpact_ForcingOp.hpp"


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
#include "Pimpact_DivGradO2JSmoother.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_InterpolateS2VOp.hpp"
#include "Pimpact_InterpolateV2SOp.hpp"
#include "Pimpact_TransferOp.hpp"
#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_ConvectionDiffusionJSmoother.hpp"
#include "Pimpact_ConvectionDiffusionSORSmoother.hpp"
#include "Pimpact_ConvectionVOp.hpp"
#include "Pimpact_ConvectionVSmoother.hpp"
#include "Pimpact_ConvectionVWrap.hpp"
#include "Pimpact_ConvectionJacobianOp.hpp"
#include "Pimpact_VectorFieldOpWrap.hpp"




/// \defgroup TimeOperator Time Operator
/// \ingroup Operator
/// Operator operate on time fields: \c Pimpact::TimeField
#include "Pimpact_DtTimeOp.hpp"
#include "Pimpact_TimeNonlinearJacobianOp.hpp"
#include "Pimpact_TimeOpWrap.hpp"
#include "Pimpact_TimeDtConvectionDiffusionOp.hpp"


/// \defgroup ModeOperator Mode Operator
/// \ingroup Operator
/// Operator operate on mode fields: \c Pimpact::ModeField
#include "Pimpact_DtModeOp.hpp"
#include "Pimpact_DtLapOp.hpp"
#include "Pimpact_EddyPrec.hpp"
#include "Pimpact_ModeOpWrap.hpp"

/// \defgroup MultiOperator Multi Operator
/// \ingroup Operator
/// Operator operate on multi fields: \c Pimpact::MultiField
#include "Pimpact_MultiOpWrap.hpp"
#include "Pimpact_MultiOpUnWrap.hpp"

/// \defgroup CompoundOperator Compound Operator
/// \ingroup Operator
/// Operator operate on compound fields: \c Pimpact::CompoundField
#include "Pimpact_CompoundOpWrap.hpp"
#include "Pimpact_InverseSchurOp.hpp"
#include "Pimpact_InverseTriangularOp.hpp"
#include "Pimpact_TimeStokesOp.hpp"
#include "Pimpact_TimeStokesBSmoother.hpp"
#include "Pimpact_TimeStokesLSmoother.hpp"
#include "Pimpact_TimeNSOp.hpp"
#include "Pimpact_TimeNSBSmoother.hpp"

/// \defgroup MultiHarmonicOperator Multi-harmonic Operator
/// \ingroup Operator
/// Operator operate on multi-harmonic fields: \c Pimpact::MultiHarmonicField
#include "Pimpact_MultiDtHelmholtzOp.hpp"
#include "Pimpact_MultiConvectionOp.hpp"
#include "Pimpact_MultiConvectionJacobianOp.hpp"
#include "Pimpact_MultiHarmonicOpWrap.hpp"
#include "Pimpact_MultiHarmonicMultiOpWrap.hpp"
#include "Pimpact_MultiHarmonicDiagOp.hpp"
#include "Pimpact_MultiDiagConvectionOp.hpp"
#include "Pimpact_MultiDtConvectionDiffusionOp.hpp"


#endif // end of #ifndef PIMPACT_OPERATOR_HPP
