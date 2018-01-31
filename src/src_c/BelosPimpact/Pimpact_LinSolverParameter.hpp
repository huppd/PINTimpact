#pragma once
#ifndef PIMPACT_LINSOLVERPARAMETER_HPP
#define PIMPACT_LINSOLVERPARAMETER_HPP


#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "BelosTypes.hpp"




namespace Pimpact {


Teuchos::RCP<Teuchos::ParameterList > createLinSolverParameter(
  const std::string& solver_name="GMRES",
  double tol=1.e-6,
  int outfreq=100,
  const Teuchos::RCP<std::ostream>& outStream=Teuchos::rcp(&std::cout, false),
  int maxIters=100,
  const std::string& label="HelloWorld");



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_LINSOLVERPARAMETER_HPP
