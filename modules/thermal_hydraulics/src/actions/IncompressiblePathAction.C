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

registerMooseAction("ThermalHydraulicsApp", IncompressiblePathAction, "add_kernel");
registerMooseAction("ThermalHydraulicsApp", IncompressiblePathAction, "add_variable");

InputParameters
IncompressiblePathAction::validParams()
{
  InputParameters params = Action::validParams();
  params += FunctorInterface::validParams();
  params += IncompressibleMomentumSPBase::validParams();
  params.addClassDescription("Action implements a path-integrated incompressible flow solve.");
  params
      .addParam<bool>("operate_on_mass",
                      true,
                      "Whether to set up an ordinary (operating on mass) momentum kernel or a "
                      "coupled pressure (operating on reference pressure drop) momentum kernel.")
          params.addRequiredParam<std::vector<NonlinearVariableName>>(
              "temperature_vars",
              "The names of the segment temperature variables in the simulation");
  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "wall_temperature_vars",
      "The names of the segment wall temperature variables in the simulation");

  return params;
}

IncompressiblePathAction::IncompressiblePathAction(const InputParameters & params)
  : Action(params),
    FunctorInterface(this),
    IncompressibleMomentumSPBase(params),
    _mc(getParam<std::vector<NonlinearVariableName>>("mass_flow_rate")),
    _dPc(getParam<std::vector<NonlinearVariableName>>("reference_pressure_drop")),
    _n_temps(getParam<std::vector<NonlinearVariableName>>("temperature_vars")),
    _n_wall_temps(getParam<std::vector<NonlinearVariableName>>("wall_temperature_vars")),
    _n_segments(this->template getParam<std::vector<MooseFunctorName>>("areas").size())
{
  const auto & area_names = MooseBase::getParam<std::vector<MooseFunctorName>>("areas");
  const auto & perimeter_names = MooseBase::getParam<std::vector<MooseFunctorName>>("perimeters");
  const auto & length_names = MooseBase::getParam<std::vector<MooseFunctorName>>("lengths");
  const auto & alpha_names = MooseBase::getParam<std::vector<MooseFunctorName>>("alphas");
  const auto & forms_loss_names =
      MooseBase::getParam<std::vector<MooseFunctorName>>("forms_losses");
  const auto & dPp_names = MooseBase::getParam<std::vector<MooseFunctorName>>("pump_pressures");
  const auto & roughness_names = MooseBase::getParam<std::vector<MooseFunctorName>>("roughnesses");
  if (_n_segments != area_names.size() || _n_segments != perimeter_names.size() ||
      _n_segments != length_names.size() || _n_segments != alpha_names.size() ||
      _n_segments != forms_loss_names.size() || _n_segments != dPp_names.size() ||
      _n_segments != roughness_names.size() || _n_segments != _n_temps - 2 ||
      _n_segments != _n_wall_temps)
  {
    mooseError(
        "Must provide consistent number of segments for each parameter! Including wall "
        "temperatures (need 2 additional segment temperatures for inlet + outlet temperatures)!");
  }
}

void
IncompressiblePathAction::act()
{
  std::vector<NonlinearVariableName> temperatures =
      getParam<std::vector<NonlinearVariableName>>("temperature_vars");
  std::vector<NonlinearVariableName> wall_temperatures =
      getParam<std::vector<NonlinearVariableName>>("wall_temperature_vars");

  // Do some error checking
  mooseAssert(variables.size() == 2, "Expected 2 variables, received " << variables.size());

  // Setup momentum kernel
  {
    InputParameters params = _factory.getValidParams("Diffusion");
    params.set<NonlinearVariableName>("variable") = variables[0];
    _problem->addKernel("Diffusion", "diff_u", params);
  }

  // Setup our Convection Kernel on the "u" variable coupled to the diffusion variable "v"
  {
    InputParameters params = _factory.getValidParams("ExampleConvection");
    params.set<NonlinearVariableName>("variable") = variables[0];
    //    params.addCoupledVar("some_variable", "The gradient of this var");
    vel_vec_variable.push_back(variables[1]);
    params.set<std::vector<VariableName>>("some_variable") = vel_vec_variable;
    _problem->addKernel("ExampleConvection", "conv", params);
  }

  // Setup out Diffusion Kernel on the "v" variable
  {
    InputParameters params = _factory.getValidParams("Diffusion");
    params.set<NonlinearVariableName>("variable") = variables[1];
    _problem->addKernel("Diffusion", "diff_v", params);
  }
}
