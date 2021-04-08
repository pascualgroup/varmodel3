# varmodel3

Ed Baskerville<br>
Frédéric Labbé

Var gene evolution model(s), implemented in Julia.

Based on previous C++ implementations (`varmodel` and `varmodel2`) by Ed Baskerville and Qixin He.

Model details are described inline in comments in the code; see *Code Organization* below to get oriented.

## Quickstart

### Setup

Before doing anything, run

```julia
./install-packages.jl
```

to install packages required by this code.

### Single runs

To do a single ad-hoc run of the model directly in Julia, copy the script `examples/julia/run.jl` into a new experiment directory, and modify and run as described in the comment string.

To perform a run with an existing parameters file in JSON format, copy the parameters file into a new experiment directory, and use the script `varmodel3/run.jl` as described in the comment string. To see how to generate parameters from JSON, see `examples/json/generate-json.jl`.

### Parameter sweeps

To do a parameter sweep, copy the `examples/sweep` directory, and modify/run as described in the comment string in `generate-sweep.jl`.

## Model Overview

TODO

## Code Organization

TODO
