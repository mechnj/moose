//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IncompressiblePathAction.h"
#include "IncompressibleMomentumSPBase.h"
#include "IncompressibleMomentumSPScalarKernel.h"
#include "CoupledPressureIncompressibleMomentumSPScalarKernel.h"
#include "IncompressibleEnergySPScalarKernel.h"
#include "FunctorInterface.h"
#include "Factory.h"
#include "ScalarCoupleable.h"

registerMooseAction("ThermalHydraulicsApp", IncompressiblePathAction, "add_kernel");
registerMooseAction("ThermalHydraulicsApp", IncompressiblePathAction, "add_variable");

InputParameters
IncompressiblePathAction::validParams()
{
  InputParameters params = Action::validParams();
  params += FunctorInterface::validParams();
  params.addClassDescription("Action implements a path-integrated incompressible flow solve.");
  // Boolean inputs
  params
      .addRequiredParam<const bool>(
          "operate_on_mass",
          true,
          "Whether to set up an ordinary (operating on mass) momentum kernel or a "
          "coupled pressure (operating on reference pressure drop) momentum kernel.") params
      .addParam<const bool>(
          "pipe_wall_CHT", true, "Whether to set up wall conjugate heat transfer.") params
      .addParam<const bool>(
          "is_implicit",
          false,
          "Whether an explicit (previous value calculation) or implicit (current value) is used");
  // Integer inputs
  params
      .addRequiredParam<const unsigned int>(
          "segment_count", 1, "Number of thermal/geometrical segments.") params
      .addParam<const unsigned int>(
          "segment_offset",
          0,
          "Number of previously created segments (to avoid renaming variables).")
      // User Object inputs
      params.addRequiredParam<UserObjectName>("fp",
                                              "The name of the user object for fluid properties");
  params.addParam<UserObjectName>("sp", "The name of the user object for solid properties");
  // Variable inputs
  params.addRequiredParam<NonlinearVariableName>("mass_flow_rate",
                                                 "The name of the mass flow rate variable.");
  params.addRequiredParam<NonlinearVariableName>(
      "pressure_drop", "The name of the reference pressure drop variable.");
  params.addRequiredParam<NonlinearVariableName>("inlet_temperature",
                                                 "The name of the inlet temperature variable.");
  params.addRequiredParam<NonlinearVariableName>("outlet_temperature",
                                                 "The name of the outlet temperature variable.");
  params.addRequiredParam<NonlinearVariableName>(
      "inlet_wall_temperature", "The name of the inlet wall temperature variable.");
  params.addRequiredParam<NonlinearVariableName>(
      "outlet_wall_temperature", "The name of the outlet wall temperature variable.");
  params
      .addParam<std::vector<NonlinearVariableName>>(
          "wall_temperatures",
          "The names of wall temperature variables to use if pipe_wall_cht is false.")
      // Functor inputs
      params.addRequiredParam<MooseFunctorName>("reference_pressure",
                                                "system reference pressure [Pa]");
  params.addParam<std::vector<MooseFunctorName>>("wall_area_inner",
                                                 std::vector<MooseFunctorName>({}),
                                                 "Cross-sectional area of inner wall layer [m^2]");
  params.addParam<std::vector<MooseFunctorName>>("wall_area_outer",
                                                 std::vector<MooseFunctorName>({}),
                                                 "Cross-sectional area of outer wall layer [m^2]");
  params.addParam<std::vector<MooseFunctorName>>(
      "interface_perimeter",
      std::vector<MooseFunctorName>({}),
      "Perimeter of the interface between the inner and outer wall nodes [m]");
  params.addParam<std::vector<MooseFunctorName>>(
      "interface_thickness",
      std::vector<MooseFunctorName>({}),
      "Radial conduction-path distance between the inner and outer wall nodes [m]");
  params.addRequiredParam<std::vector<MooseFunctorName>>(
      "flow_areas",
      std::vector<MooseFunctorName>({}),
      "Component flow areas per segment. Takes a vector of functors.");
  params.addRequiredParam<std::vector<MooseFunctorName>>(
      "wetted_perimeters",
      std::vector<MooseFunctorName>({}),
      "Component flow perimeters per segment. Takes a vector of functors.");
  params.addRequiredParam<std::vector<MooseFunctorName>>(
      "lengths",
      std::vector<MooseFunctorName>({}),
      "Component flow lengths per segment. Takes a vector of functors.");
  params.addParam<std::vector<MooseFunctorName>>(
      "dHs",
      std::vector<MooseFunctorName>({}),
      "Component height change from inlet to outlet. Takes a vector of functors.");
  params.addParam<std::vector<MooseFunctorName>>(
      "forms_losses",
      std::vector<MooseFunctorName>({}),
      "Forms loss coefficients per segment. Takes a vector of functors.");
  params.addParam<std::vector<MooseFunctorName>>(
      "pump_pressures",
      std::vector<MooseFunctorName>({}),
      "Pump pressure gains per segment [Pa]. Takes a vector of functors.");
  params.addParam<std::vector<MooseFunctorName>>(
      "roughnesses",
      std::vector<MooseFunctorName>({}),
      "Component wall roughnesses per segment [m]. Takes a vector of functors.");
  params.addParam<MooseFunctorName>(
      "g", PhysicalConstants::acceleration_of_gravity, "Gravitational acceleration [m/s]");

  return params;
}

IncompressiblePathAction::IncompressiblePathAction(const InputParameters & params)
  : Action(params), FunctorInterface(this)
{
}

void
IncompressiblePathAction::act()
{
  // Initialize parameters
  const bool regular_momentum = getParam<const bool>("operate_on_mass");
  const bool wall_cht = getParam<const bool>("pipe_wall_CHT");
  const unsigned int n_seg = getParam<const unsigned int>("segment_count");
  const unsigned int seg_off = getParam<const unsigned int>("segment_offset");
  std::vector<NonlinearVariableName> walltemps =
      getParam<std::vector<NonlinearVariableName>>("wall_temperatures");
  std::vector<MooseFunctorName> flow_areas =
      MooseBase::getParam<std::vector<MooseFunctorName>>("flow_areas");
  std::vector<MooseFunctorName> wetted_perimeters =
      MooseBase::getParam<std::vector<MooseFunctorName>>("wetted_perimeters");
  std::vector<MooseFunctorName> lengths =
      MooseBase::getParam<std::vector<MooseFunctorName>>("lengths");

  // Add temperature variables
  if (_current_task == "add_variable")
  {
    auto fe_type = AddVariableAction::feType(_pars);
    auto type = AddVariableAction::variableType(fe_type, false, array);
    auto var_params = _factory.getValidParams(type);
    var_params.set<MooseEnum>("family") = "SCALAR";
    // Add fluid temps
    for (unsigned int i = 0; i < n_seg; i++)
    {
      std::string var_name = "T" + Moose::stringify(i + seg_off);
      _problem->addVariable(type, var_name, var_params);
    }
    // Add wall temps
    if (wall_cht)
    {
      for (unsigned int i = 0; i < n_seg; i++)
      {
        std::string var_name = "Tw" + Moose::stringify(i + seg_off);
        _problem->addVariable(type, var_name, var_params);
      }
    }
  }
  else if (_current_task == "add_kernel")
  {
    // Add regular momentum kernel
    if (regular_momentum)
    {
      auto kernel_type = "ADIncompressibleMomentumSPScalarKernel";
      InputParameters params =
          _factory.getValidParams(kernel_type) params.set<NonlinearVariableName>("variable") =
              getParam<NonlinearVariableName>("mass_flow_rate");
      params.set<ScalarCoupleable::coupledScalarComponents>("reference_pressure_drop") =
          getParam<NonlinearVariableName>("reference_pressure_drop");
      std::vector<std::string> temps = {} for (unsigned int i = 0; i < n_seg; i++)
      {
        std::string var_name = "T" + Moose::stringify(i + seg_off);
        temps.push_back(var_name);
      }
      params.set<ScalarCoupleable::coupledScalarComponents>("temperatures") = temps;
      params.set<bool>("is_implicit") = getParam<bool>("is_implicit");
      params.set<UserObjectName>("fp") = getParam<UserObjectName>("fp");
      params.set<MooseFunctorName>("reference_pressure") =
          getParam<MooseFunctorName>("reference_pressure");
      params.set<std::vector<MooseFunctorName>>("areas") =
          getParam<std::vector<MooseFunctorName>>("flow_areas");
      params.set<std::vector<MooseFunctorName>>("perimeters") =
          getParam<std::vector<MooseFunctorName>>("wetted_perimeters");
      params.set<std::vector<MooseFunctorName>>("lengths") =
          getParam<std::vector<MooseFunctorName>>("lengths");
      params.set<std::vector<MooseFunctorName>>("dHs") =
          getParam<std::vector<MooseFunctorName>>("dHs");
      params.set<std::vector<MooseFunctorName>>("forms_losses") =
          getParam<std::vector<MooseFunctorName>>("forms_losses");
      params.set<std::vector<MooseFunctorName>>("pump_pressures") =
          getParam<std::vector<MooseFunctorName>>("pump_pressures");
      params.set<std::vector<MooseFunctorName>>("roughnesses") =
          getParam<std::vector<MooseFunctorName>>("roughnesses");
      params.set<MooseFunctorName>("g") = getParam<MooseFunctorName>("g");
      _problem->addKernel(kernel_type, getParam<NonlinearVariableName>("mass_flow_rate"), params);
    }
    // Add coupled pressure momentum kernel
    else
    {
      auto kernel_type = "ADCoupledPressureIncompressibleMomentumSPScalarKernel";
      InputParameters params =
          _factory.getValidParams(kernel_type) params.set<NonlinearVariableName>("variable") =
              getParam<NonlinearVariableName>("reference_pressure_drop");
      params.set<ScalarCoupleable::coupledScalarComponents>("coupled_mass_flow_rate") =
          getParam<NonlinearVariableName>("mass_flow_rate");
      std::vector<std::string> temps = {} for (unsigned int i = 0; i < n_seg; i++)
      {
        std::string var_name = "T" + Moose::stringify(i + seg_off);
        temps.push_back(var_name);
      }
      params.set<ScalarCoupleable::coupledScalarComponents>("temperatures") = temps;
      params.set<bool>("is_implicit") = getParam<bool>("is_implicit");
      params.set<UserObjectName>("fp") = getParam<UserObjectName>("fp");
      params.set<MooseFunctorName>("reference_pressure") =
          getParam<MooseFunctorName>("reference_pressure");
      params.set<std::vector<MooseFunctorName>>("areas") =
          getParam<std::vector<MooseFunctorName>>("flow_areas");
      params.set<std::vector<MooseFunctorName>>("perimeters") =
          getParam<std::vector<MooseFunctorName>>("wetted_perimeters");
      params.set<std::vector<MooseFunctorName>>("lengths") =
          getParam<std::vector<MooseFunctorName>>("lengths");
      params.set<std::vector<MooseFunctorName>>("dHs") =
          getParam<std::vector<MooseFunctorName>>("dHs");
      params.set<std::vector<MooseFunctorName>>("forms_losses") =
          getParam<std::vector<MooseFunctorName>>("forms_losses");
      params.set<std::vector<MooseFunctorName>>("pump_pressures") =
          getParam<std::vector<MooseFunctorName>>("pump_pressures");
      params.set<std::vector<MooseFunctorName>>("roughnesses") =
          getParam<std::vector<MooseFunctorName>>("roughnesses");
      params.set<MooseFunctorName>("g") = getParam<MooseFunctorName>("g");
      _problem->addKernel(
          kernel_type, getParam<NonlinearVariableName>("reference_pressure_drop"), params);
    }
    // Add inlet temperature kernel
    auto kernel_type = "ADIncompressibleEnergySPScalarKernel";
    InputParameters params =
        _factory.getValidParams(kernel_type) params.set<NonlinearVariableName>("variable") =
            "T" + Moose::stringify(seg_off);
    params.set<ScalarCoupleable::coupledScalarComponents>("mass_flow_rate") =
        getParam<NonlinearVariableName>("mass_flow_rate");
    params.set<ScalarCoupleable::coupledScalarComponents>("inlet_temperature") =
        getParam<NonlinearVariableName>("inlet_temperature");
    params.set<ScalarCoupleable::coupledScalarComponents>("outlet_temperature") =
        "T" + Moose::stringify(1 + seg_off);
    if (wall_cht)
    {
      params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") =
          "Tw" + Moose::stringify(seg_off);
    }
    else
    {
      params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") = walltemps[0];
    }
    params.set<bool>("is_implicit") = getParam<bool>("is_implicit");
    params.set<UserObjectName>("fp") = getParam<UserObjectName>("fp");
    params.set<MooseFunctorName>("reference_pressure") =
        getParam<MooseFunctorName>("reference_pressure");
    params.set<MooseFunctorName>("area") = flow_areas[0];
    params.set<MooseFunctorName>("perimeter") = wetted_perimeters[0];
    params.set<MooseFunctorName>("length") = lengths[0];
    _problem->addKernel(kernel_type, "T" + Moose::stringify(seg_off);, params);
    // Add outlet temperature kernel
    auto kernel_type = "ADIncompressibleEnergySPScalarKernel";
    InputParameters params =
        _factory.getValidParams(kernel_type) params.set<NonlinearVariableName>("variable") =
            "T" + Moose::stringify(seg_off + n_seg);
    params.set<ScalarCoupleable::coupledScalarComponents>("mass_flow_rate") =
        getParam<NonlinearVariableName>("mass_flow_rate");
    params.set<ScalarCoupleable::coupledScalarComponents>("inlet_temperature") =
        "T" + Moose::stringify(seg_off + n_seg - 1);
    params.set<ScalarCoupleable::coupledScalarComponents>("outlet_temperature") =
        getParam<NonlinearVariableName>("outlet_temperature");
    if (wall_cht)
    {
      params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") =
          "Tw" + Moose::stringify(seg_off + n_seg);
    }
    else
    {
      params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") =
          walltemps[walltemps.size()];
    }
    params.set<bool>("is_implicit") = getParam<bool>("is_implicit");
    params.set<UserObjectName>("fp") = getParam<UserObjectName>("fp");
    params.set<MooseFunctorName>("reference_pressure") =
        getParam<MooseFunctorName>("reference_pressure");
    params.set<MooseFunctorName>("area") = flow_areas[flow_areas.size()];
    params.set<MooseFunctorName>("perimeter") = wetted_perimeters[wetted_perimeters.size()];
    params.set<MooseFunctorName>("length") = lengths[lengths.size()];
    _problem->addKernel(kernel_type, "T" + Moose::stringify(seg_off + n_seg);, params);
    // Add all other temperature kernels
    for (unsigned int i = 1; i < n_seg - 1; i++)
    {
      auto kernel_type = "ADIncompressibleEnergySPScalarKernel";
      InputParameters params =
          _factory.getValidParams(kernel_type) params.set<NonlinearVariableName>("variable") =
              "T" + Moose::stringify(seg_off + i);
      params.set<ScalarCoupleable::coupledScalarComponents>("mass_flow_rate") =
          getParam<NonlinearVariableName>("mass_flow_rate");
      params.set<ScalarCoupleable::coupledScalarComponents>("inlet_temperature") =
          "T" + Moose::stringify(seg_off + i - 1);
      params.set<ScalarCoupleable::coupledScalarComponents>("outlet_temperature") =
          "T" + Moose::stringify(seg_off + i + 1);
      if (wall_cht)
      {
        params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") =
            "Tw" + Moose::stringify(seg_off + i);
      }
      else
      {
        params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") = walltemps[i];
      }
      params.set<bool>("is_implicit") = getParam<bool>("is_implicit");
      params.set<UserObjectName>("fp") = getParam<UserObjectName>("fp");
      params.set<MooseFunctorName>("reference_pressure") =
          getParam<MooseFunctorName>("reference_pressure");
      params.set<MooseFunctorName>("area") = flow_areas[i];
      params.set<MooseFunctorName>("perimeter") = wetted_perimeters[i];
      params.set<MooseFunctorName>("length") = lengths[i];
      _problem->addKernel(kernel_type, "T" + Moose::stringify(seg_off + i);, params);
    }
    if (wall_cht)
    {
      // Add inlet wall temperature kernel
      auto kernel_type = "ADInnerPipeWallTemperatureScalarKernel";
      InputParameters params =
          _factory.getValidParams(kernel_type) params.set<NonlinearVariableName>("variable") =
              "Tw" + Moose::stringify(seg_off);
      params.set<ScalarCoupleable::coupledScalarComponents>("mass_flow_rate") =
          getParam<NonlinearVariableName>("mass_flow_rate");
      params.set<ScalarCoupleable::coupledScalarComponents>("inlet_temperature") =
          getParam<NonlinearVariableName>("inlet_temperature");
      params.set<ScalarCoupleable::coupledScalarComponents>("outlet_temperature") =
          "T" + Moose::stringify(1 + seg_off);
      if (wall_cht)
      {
        params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") =
            "Tw" + Moose::stringify(seg_off);
      }
      else
      {
        params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") = walltemps[0];
      }
      params.set<bool>("is_implicit") = getParam<bool>("is_implicit");
      params.set<UserObjectName>("fp") = getParam<UserObjectName>("fp");
      params.set<MooseFunctorName>("reference_pressure") =
          getParam<MooseFunctorName>("reference_pressure");
      params.set<MooseFunctorName>("area") = flow_areas[0];
      params.set<MooseFunctorName>("perimeter") = wetted_perimeters[0];
      params.set<MooseFunctorName>("length") = lengths[0];
      _problem->addKernel(kernel_type, "T" + Moose::stringify(seg_off);, params);
      // Add outlet temperature kernel
      auto kernel_type = "ADIncompressibleEnergySPScalarKernel";
      InputParameters params =
          _factory.getValidParams(kernel_type) params.set<NonlinearVariableName>("variable") =
              "T" + Moose::stringify(seg_off + n_seg);
      params.set<ScalarCoupleable::coupledScalarComponents>("mass_flow_rate") =
          getParam<NonlinearVariableName>("mass_flow_rate");
      params.set<ScalarCoupleable::coupledScalarComponents>("inlet_temperature") =
          "T" + Moose::stringify(seg_off + n_seg - 1);
      params.set<ScalarCoupleable::coupledScalarComponents>("outlet_temperature") =
          getParam<NonlinearVariableName>("outlet_temperature");
      if (wall_cht)
      {
        params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") =
            "Tw" + Moose::stringify(seg_off + n_seg);
      }
      else
      {
        params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") =
            walltemps[walltemps.size()];
      }
      params.set<bool>("is_implicit") = getParam<bool>("is_implicit");
      params.set<UserObjectName>("fp") = getParam<UserObjectName>("fp");
      params.set<MooseFunctorName>("reference_pressure") =
          getParam<MooseFunctorName>("reference_pressure");
      params.set<MooseFunctorName>("area") = flow_areas[flow_areas.size()];
      params.set<MooseFunctorName>("perimeter") = wetted_perimeters[wetted_perimeters.size()];
      params.set<MooseFunctorName>("length") = lengths[lengths.size()];
      _problem->addKernel(kernel_type, "T" + Moose::stringify(seg_off + n_seg);, params);
      // Add all other temperature kernels
      for (unsigned int i = 1; i < n_seg - 1; i++)
      {
        auto kernel_type = "ADIncompressibleEnergySPScalarKernel";
        InputParameters params =
            _factory.getValidParams(kernel_type) params.set<NonlinearVariableName>("variable") =
                "T" + Moose::stringify(seg_off + i);
        params.set<ScalarCoupleable::coupledScalarComponents>("mass_flow_rate") =
            getParam<NonlinearVariableName>("mass_flow_rate");
        params.set<ScalarCoupleable::coupledScalarComponents>("inlet_temperature") =
            "T" + Moose::stringify(seg_off + i - 1);
        params.set<ScalarCoupleable::coupledScalarComponents>("outlet_temperature") =
            "T" + Moose::stringify(seg_off + i + 1);
        if (wall_cht)
        {
          params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") =
              "Tw" + Moose::stringify(seg_off + i);
        }
        else
        {
          params.set<ScalarCoupleable::coupledScalarComponents>("wall_temperature") = walltemps[i];
        }
        params.set<bool>("is_implicit") = getParam<bool>("is_implicit");
        params.set<UserObjectName>("fp") = getParam<UserObjectName>("fp");
        params.set<MooseFunctorName>("reference_pressure") =
            getParam<MooseFunctorName>("reference_pressure");
        params.set<MooseFunctorName>("area") = flow_areas[i];
        params.set<MooseFunctorName>("perimeter") = wetted_perimeters[i];
        params.set<MooseFunctorName>("length") = lengths[i];
        _problem->addKernel(kernel_type, "T" + Moose::stringify(seg_off + i);, params);
      }
    }
  }
}
