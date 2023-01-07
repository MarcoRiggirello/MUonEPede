# Marco Riggirello

"""
CMS/MUonE BT data preprocessing for alignment with Millepede II
"""
module MUonEPede

using LinearAlgebra, Rotations, StaticArrays
using UnROOT, FortranFiles, EzXML
using Optim, LineSearches
using PythonCall
using ProgressBars

export MUonEModule, MUonEStation, Stub, Track
export strip_to_local, local_to_global, global_to_local
export intersection, interpolate
export trackfit
export getmodules, generatebin, generatebinmc
export residuals

include("types.jl")
include("io.jl")
include("intersection.jl")
include("fit.jl")
include("mille.jl")
include("transforms.jl")

include("montecarlo.jl")
include("residuals.jl")

end
