#pragma once
#ifndef PIMPACT_OPERATOR_HPP
#define PIMPACT_OPERATOR_HPP

/// \defgroup Operator Operators

/// \defgroup BaseOperator Basic Operator
/// \ingroup Operator
/// Operator operate on basic fields, for example \c Pimpact::ScalarField and \c Pimpact::VectorField

/// \defgroup ModeOperator Mode Operator
/// \ingroup Operator
/// Operator operate on mode fields: \c Pimpact::ModeField

/// \defgroup MultiOperator Multi Operator
/// \ingroup Operator
/// Operator operate on multi fields: \c Pimpact::MultiField

/// \defgroup CompoundOperator Compound Operator
/// \ingroup Operator
/// Operator operate on compound fields: \c Pimpact::CompoundField

#include "Pimpact_CompoundOp.hpp"
#include "Pimpact_CompoundStokesOp.hpp"

#include "Pimpact_DivDtLinvGradOp.hpp"
#include "Pimpact_DivGradOp.hpp"
#include "Pimpact_DivHinvGradOp.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_DivOpGrad.hpp"

#include "Pimpact_DtHelmholtzOp.hpp"
#include "Pimpact_DtOp.hpp"

#include "Pimpact_EddyPrec.hpp"

#include "Pimpact_GradOp.hpp"

#include "Pimpact_HelmholtzInverseOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"

#include "Pimpact_InverseOperator.hpp"

#include "Pimpact_ModeOpWrap.hpp"
#include "Pimpact_MultiOpWrap.hpp"

#include "Pimpact_Nonlinear.hpp"





#endif // end of #ifndef PIMPACT_OPERATOR_HPP
