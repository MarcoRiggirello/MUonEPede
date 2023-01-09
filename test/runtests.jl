include("../src/MUonEPede.jl")

using .MUonEPede

using Test
using StaticArrays

modules = getmodules("../MUonEStructure_December2022.xml")

include("intersection_test.jl")
#include("mille_test.jl")
include("transforms_test.jl")
#include("types_test.jl")
include("fit_test.jl")
