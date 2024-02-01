#!/usr/bin/env julia

"""
This script installs packages required by the model code.

This process is not automatic. If the code is modified to require more
packages, they need to be added manually here.
"""

import Pkg

Pkg.add("StaticArrays")
Pkg.add("StatsBase")
Pkg.add("Distributions")
Pkg.add("Parameters")
Pkg.add("JSON")
Pkg.add("JSON3")
Pkg.add("StructTypes")
Pkg.add("SQLite")
Pkg.add("StatsAPI")
Pkg.add("HypothesisTests")
