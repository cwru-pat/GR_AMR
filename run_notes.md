# Welcome to AMR version of CosmoGRaPH

Notes on run CosmoGRaPH with different matter fields

## General notes

A sample input file involving a static black hole is in input/static_blackhole.input.
It also contains comments on each parameter used.

Macros in file cosmo_macros.h need to be changed when running with different matter field
with enabling different features. See comments inside cosmo_macros.h


## Running vaccum simulation

Some notes on writing input files:

```
Main
{
...
  simulation_type = "vacuum"
...
}
..

VacuumSim{
  ic_type = "static_blackhole" or "kerr_blackhole" or "EdS" or awa_linear_wave
  
  boundary_type = "sommerfield" or "periodic"
}


```
User defined initial type can be made by adding new function in components/static/static_ic.cc
and add its entry in sims/vacuum.cc.

To stably evolve spacetime containing blackhole, gamma-driver gauge is recommended to use.
In order to employ this gauge, macros USE_BSSN_SHIFT and USE_GAMMA_DRIVER should be set to TRUE.
Proper gauge choices need to be set in the input file.

## Running scalar simulation

Some notes on writing input files:

```
Main
{
...
  simulation_type = "scalar"
...
}
...

ScalarSim{
  ic_type = "scalar_collapse", or "scalar_gaussian_random"
  boundary_type = "periodic"
  stop_after_setting_init = FALSE
}

...

Scalar{
...
// Could work with the name of parameters for different choices
// from scalarPotentialHandler.cc
  potential_type = "Constant", "Quadratic", "Exp_p"
...
}

```
User defined initial type can be made by adding new function in components/scalar/scalar.cc
and add its entry in sims/scalar.cc
