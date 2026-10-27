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
 * Fluid properties provided as a mixture of two other fluid properties with AD overrides
 */
class BlendedFluidProperties : public SinglePhaseFluidProperties
{
public:
  BlendedFluidProperties(const InputParameters & parameters);
  static InputParameters validParams();

  virtual ADReal phi() const;

  /**
   * Density from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @return density (kg/m$^3$)
   */
  virtual ADReal rho_from_p_T(ADReal p, ADReal T) const;

  /**
   * isobaric specific heat from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @return cp (J/kg/.K)
   */
  virtual ADReal cp_from_p_T(ADReal p, ADReal T) const;

  /**
   * thermal conductivity from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @return k (W/m./K)
   */
  virtual ADReal k_from_p_T(ADReal p, ADReal T) const;

  /**
   * dynamic viscosity from pressure and temperature
   *
   * @param[in] p   pressure (Pa)
   * @param[in] T   temperature (K)
   * @return mu (Pa.s)
   */
  virtual ADReal mu_from_p_T(ADReal p, ADReal T) const;

protected:
  const SinglePhaseFluidProperties &_l, &_g;
  const PostprocessorValue & _phi;
};
