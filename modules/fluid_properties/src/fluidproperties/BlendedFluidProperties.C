//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlendedFluidProperties.h"

registerMooseObject("FluidPropertiesApp", BlendedFluidProperties);

InputParameters
BlendedFluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();
  params.addRequiredParam<UserObjectName>("fp_liquid", "Liquid-phase properties");
  params.addRequiredParam<UserObjectName>("fp_gas", "Gas-phase properties");
  params.addRequiredParam<PostprocessorName>(
      "liquid_fraction", "Postprocessor giving the liquid volume fraction in [0,1]");
  return params;
}

BlendedFluidProperties::BlendedFluidProperties(const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters),
    _l(getUserObject<SinglePhaseFluidProperties>("fp_liquid")),
    _g(getUserObject<SinglePhaseFluidProperties>("fp_gas")),
    _phi(getPostprocessorValue("liquid_fraction"))
{
}

ADReal
BlendedFluidProperties::phi() const
{
  return std::clamp(_phi, 0.0, 1.0);
}

ADReal
BlendedFluidProperties::rho_from_p_T(ADReal p, ADReal T) const
{
  const ADReal f = phi();
  return f * _l.rho_from_p_T(p, T) + (1 - f) * _g.rho_from_p_T(p, T);
}

ADReal
BlendedFluidProperties::cp_from_p_T(ADReal p, ADReal T) const
{
  const ADReal f = phi();
  const auto rl = _l.rho_from_p_T(p, T), rg = _g.rho_from_p_T(p, T);
  return (f * rl * _l.cp_from_p_T(p, T) + (1 - f) * rg * _g.cp_from_p_T(p, T)) /
         (f * rl + (1 - f) * rg);
}

ADReal
BlendedFluidProperties::mu_from_p_T(ADReal p, ADReal T) const
{
  const ADReal f = phi();
  return exp(f * log(_l.mu_from_p_T(p, T)) + (1 - f) * log(_g.mu_from_p_T(p, T)));
}

ADReal
BlendedFluidProperties::k_from_p_T(ADReal p, ADReal T) const
{
  const ADReal f = phi();
  return f * _l.k_from_p_T(p, T) + (1 - f) * _g.k_from_p_T(p, T);
}
