"""
This file is loaded by run scripts before constructing parameters.

Prioritizing simplicity over precision, package requirements are specified
here and shared across all other code.
"""

using StaticArrays

using Dates

using Random
using StatsBase
using Distributions

using Parameters
using JSON
using DelimitedFiles
using SQLite

include("src/parameters.jl")
