//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlendedFluidProperties.h"

#include <algorithm>
#include <cmath>

registerMooseObject("FluidPropertiesApp", BlendedFluidProperties);

InputParameters
BlendedFluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();
  params.addRequiredParam<UserObjectName>("fp_liquid", "Liquid-phase fluid properties object");
  params.addRequiredParam<UserObjectName>("fp_gas", "Gas-phase fluid properties object");
  params.addRequiredParam<PostprocessorName>(
      "liquid_fraction", "Postprocessor giving the liquid volume fraction in [0, 1]");
  params.addClassDescription(
      "Fluid properties from blending a liquid and a gas fluid properties object using a liquid "
      "volume fraction.");
  return params;
}

BlendedFluidProperties::BlendedFluidProperties(const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters),
    _fp_liquid(getUserObject<SinglePhaseFluidProperties>("fp_liquid")),
    _fp_gas(getUserObject<SinglePhaseFluidProperties>("fp_gas")),
    _liquid_fraction(getPostprocessorValue("liquid_fraction"))
{
}

std::string
BlendedFluidProperties::fluidName() const
{
  return "blend(" + _fp_liquid.fluidName() + ", " + _fp_gas.fluidName() + ")";
}

Real
BlendedFluidProperties::phi() const
{
  return std::min(std::max(_liquid_fraction, Real(0)), Real(1));
}

// ---------------------------------------------------------------------------------------------
// rho(p, T)
// ---------------------------------------------------------------------------------------------

Real
BlendedFluidProperties::rho_from_p_T(Real p, Real T) const
{
  return rho_from_p_T_template(p, T);
}

void
BlendedFluidProperties::rho_from_p_T(
    Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho_from_p_T_template(p, T, rho, drho_dp, drho_dT);
}

ADReal
BlendedFluidProperties::rho_from_p_T(const ADReal & p, const ADReal & T) const
{
  return rho_from_p_T_template(p, T);
}

void
BlendedFluidProperties::rho_from_p_T(
    const ADReal & p, const ADReal & T, ADReal & rho, ADReal & drho_dp, ADReal & drho_dT) const
{
  rho_from_p_T_template(p, T, rho, drho_dp, drho_dT);
}

template <typename CppType>
CppType
BlendedFluidProperties::rho_from_p_T_template(const CppType & p, const CppType & T) const
{
  const Real f = phi();
  return f * _fp_liquid.rho_from_p_T(p, T) + (1 - f) * _fp_gas.rho_from_p_T(p, T);
}

template <typename CppType>
void
BlendedFluidProperties::rho_from_p_T_template(
    const CppType & p, const CppType & T, CppType & rho, CppType & drho_dp, CppType & drho_dT) const
{
  const Real f = phi();
  CppType rl, drl_dp, drl_dT, rg, drg_dp, drg_dT;
  _fp_liquid.rho_from_p_T(p, T, rl, drl_dp, drl_dT);
  _fp_gas.rho_from_p_T(p, T, rg, drg_dp, drg_dT);

  rho = f * rl + (1 - f) * rg;
  drho_dp = f * drl_dp + (1 - f) * drg_dp;
  drho_dT = f * drl_dT + (1 - f) * drg_dT;
}

// ---------------------------------------------------------------------------------------------
// cp(p, T): mass-weighted
// ---------------------------------------------------------------------------------------------

Real
BlendedFluidProperties::cp_from_p_T(Real p, Real T) const
{
  return cp_from_p_T_template(p, T);
}

void
BlendedFluidProperties::cp_from_p_T(Real p, Real T, Real & cp, Real & dcp_dp, Real & dcp_dT) const
{
  cp_from_p_T_template(p, T, cp, dcp_dp, dcp_dT);
}

ADReal
BlendedFluidProperties::cp_from_p_T(const ADReal & p, const ADReal & T) const
{
  return cp_from_p_T_template(p, T);
}

void
BlendedFluidProperties::cp_from_p_T(
    const ADReal & p, const ADReal & T, ADReal & cp, ADReal & dcp_dp, ADReal & dcp_dT) const
{
  cp_from_p_T_template(p, T, cp, dcp_dp, dcp_dT);
}

template <typename CppType>
CppType
BlendedFluidProperties::cp_from_p_T_template(const CppType & p, const CppType & T) const
{
  const Real f = phi();
  const CppType rl = _fp_liquid.rho_from_p_T(p, T);
  const CppType rg = _fp_gas.rho_from_p_T(p, T);
  const CppType cpl = _fp_liquid.cp_from_p_T(p, T);
  const CppType cpg = _fp_gas.cp_from_p_T(p, T);
  return (f * rl * cpl + (1 - f) * rg * cpg) / (f * rl + (1 - f) * rg);
}

template <typename CppType>
void
BlendedFluidProperties::cp_from_p_T_template(
    const CppType & p, const CppType & T, CppType & cp, CppType & dcp_dp, CppType & dcp_dT) const
{
  const Real f = phi();
  CppType rl, drl_dp, drl_dT, rg, drg_dp, drg_dT;
  CppType cpl, dcpl_dp, dcpl_dT, cpg, dcpg_dp, dcpg_dT;
  _fp_liquid.rho_from_p_T(p, T, rl, drl_dp, drl_dT);
  _fp_gas.rho_from_p_T(p, T, rg, drg_dp, drg_dT);
  _fp_liquid.cp_from_p_T(p, T, cpl, dcpl_dp, dcpl_dT);
  _fp_gas.cp_from_p_T(p, T, cpg, dcpg_dp, dcpg_dT);

  // cp = N / D, with N = sum(phi_i rho_i cp_i) and D = sum(phi_i rho_i)
  const CppType N = f * rl * cpl + (1 - f) * rg * cpg;
  const CppType D = f * rl + (1 - f) * rg;
  cp = N / D;

  const CppType dN_dp = f * (drl_dp * cpl + rl * dcpl_dp) + (1 - f) * (drg_dp * cpg + rg * dcpg_dp);
  const CppType dN_dT = f * (drl_dT * cpl + rl * dcpl_dT) + (1 - f) * (drg_dT * cpg + rg * dcpg_dT);
  const CppType dD_dp = f * drl_dp + (1 - f) * drg_dp;
  const CppType dD_dT = f * drl_dT + (1 - f) * drg_dT;

  dcp_dp = (dN_dp - cp * dD_dp) / D;
  dcp_dT = (dN_dT - cp * dD_dT) / D;
}

// ---------------------------------------------------------------------------------------------
// k(p, T): volume-fraction weighted
// ---------------------------------------------------------------------------------------------

Real
BlendedFluidProperties::k_from_p_T(Real p, Real T) const
{
  return k_from_p_T_template(p, T);
}

void
BlendedFluidProperties::k_from_p_T(Real p, Real T, Real & k, Real & dk_dp, Real & dk_dT) const
{
  k_from_p_T_template(p, T, k, dk_dp, dk_dT);
}

ADReal
BlendedFluidProperties::k_from_p_T(const ADReal & p, const ADReal & T) const
{
  return k_from_p_T_template(p, T);
}

void
BlendedFluidProperties::k_from_p_T(
    const ADReal & p, const ADReal & T, ADReal & k, ADReal & dk_dp, ADReal & dk_dT) const
{
  k_from_p_T_template(p, T, k, dk_dp, dk_dT);
}

template <typename CppType>
CppType
BlendedFluidProperties::k_from_p_T_template(const CppType & p, const CppType & T) const
{
  const Real f = phi();
  return f * _fp_liquid.k_from_p_T(p, T) + (1 - f) * _fp_gas.k_from_p_T(p, T);
}

template <typename CppType>
void
BlendedFluidProperties::k_from_p_T_template(
    const CppType & p, const CppType & T, CppType & k, CppType & dk_dp, CppType & dk_dT) const
{
  const Real f = phi();
  CppType kl, dkl_dp, dkl_dT, kg, dkg_dp, dkg_dT;
  _fp_liquid.k_from_p_T(p, T, kl, dkl_dp, dkl_dT);
  _fp_gas.k_from_p_T(p, T, kg, dkg_dp, dkg_dT);

  k = f * kl + (1 - f) * kg;
  dk_dp = f * dkl_dp + (1 - f) * dkg_dp;
  dk_dT = f * dkl_dT + (1 - f) * dkg_dT;
}

// ---------------------------------------------------------------------------------------------
// mu(p, T): log-weighted
// ---------------------------------------------------------------------------------------------

Real
BlendedFluidProperties::mu_from_p_T(Real p, Real T) const
{
  return mu_from_p_T_template(p, T);
}

void
BlendedFluidProperties::mu_from_p_T(Real p, Real T, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  mu_from_p_T_template(p, T, mu, dmu_dp, dmu_dT);
}

ADReal
BlendedFluidProperties::mu_from_p_T(const ADReal & p, const ADReal & T) const
{
  return mu_from_p_T_template(p, T);
}

void
BlendedFluidProperties::mu_from_p_T(
    const ADReal & p, const ADReal & T, ADReal & mu, ADReal & dmu_dp, ADReal & dmu_dT) const
{
  mu_from_p_T_template(p, T, mu, dmu_dp, dmu_dT);
}

template <typename CppType>
CppType
BlendedFluidProperties::mu_from_p_T_template(const CppType & p, const CppType & T) const
{
  using std::exp;
  using std::log;
  const Real f = phi();
  return exp(f * log(_fp_liquid.mu_from_p_T(p, T)) + (1 - f) * log(_fp_gas.mu_from_p_T(p, T)));
}

template <typename CppType>
void
BlendedFluidProperties::mu_from_p_T_template(
    const CppType & p, const CppType & T, CppType & mu, CppType & dmu_dp, CppType & dmu_dT) const
{
  using std::exp;
  using std::log;
  const Real f = phi();
  CppType mul, dmul_dp, dmul_dT, mug, dmug_dp, dmug_dT;
  _fp_liquid.mu_from_p_T(p, T, mul, dmul_dp, dmul_dT);
  _fp_gas.mu_from_p_T(p, T, mug, dmug_dp, dmug_dT);

  mu = exp(f * log(mul) + (1 - f) * log(mug));

  // d(mu)/dx = mu * (phi * dmu_l/dx / mu_l + (1 - phi) * dmu_g/dx / mu_g)
  dmu_dp = mu * (f * dmul_dp / mul + (1 - f) * dmug_dp / mug);
  dmu_dT = mu * (f * dmul_dT / mul + (1 - f) * dmug_dT / mug);
}
