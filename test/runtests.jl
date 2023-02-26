include("../src/MUonEPede.jl")

using .MUonEPede

using Test
using StaticArrays
using Optim, LineSearches

modules = MUonEPede.MUonEStation{Float32}("../MUonEStructure_December2022.xml")

include("transforms_test.jl")
#include("mille_test.jl")
include("intersection_test.jl")
#include("types_test.jl")
include("fit_test.jl")
