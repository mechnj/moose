//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SinglePhaseFluidProperties.h"

/**
 * Fluid properties provided as a blend of a liquid and a gas fluid properties object, weighted by
 * a liquid volume fraction phi in [0, 1] that is read from a postprocessor.
 *
 * Blending rules:
 *  - rho: volume-fraction weighted,  rho = phi rho_l + (1 - phi) rho_g
 *  - cp : mass-weighted,             cp  = sum(phi_i rho_i cp_i) / sum(phi_i rho_i)
 *  - k  : volume-fraction weighted
 *  - mu : log-weighted (geometric),  mu  = mu_l^phi mu_g^(1 - phi)
 *
 * For each property the Real, Real-with-derivatives, ADReal and ADReal-with-derivatives
 * (p, T) interfaces are overridden together, so no overloads of the base class are hidden.
 * The liquid fraction is treated as a constant (non-AD) value for each evaluation, so no
 * derivatives with respect to it are propagated.
 */
class BlendedFluidProperties : public SinglePhaseFluidProperties
{
public:
  static InputParameters validParams();

  BlendedFluidProperties(const InputParameters & parameters);

  virtual std::string fluidName() const override;

  /// Liquid volume fraction, clamped to [0, 1]
  Real phi() const;

  // Declares Real, Real+derivatives, ADReal and ADReal+derivatives overrides plus the templated
  // implementations (see the definitions in BlendedFluidProperties.C)
  propfuncWithDefinitionOverride(rho, p, T);
  propfuncWithDefinitionOverride(cp, p, T);
  propfuncWithDefinitionOverride(k, p, T);
  propfuncWithDefinitionOverride(mu, p, T);

protected:
  /// Liquid-phase fluid properties
  const SinglePhaseFluidProperties & _fp_liquid;
  /// Gas-phase fluid properties
  const SinglePhaseFluidProperties & _fp_gas;
  /// Liquid volume fraction (raw, unclamped)
  const PostprocessorValue & _liquid_fraction;
};
