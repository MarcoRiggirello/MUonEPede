# Marco Riggirello

"""
CMS/MUonE BT data preprocessing for alignment with Millepede II
"""
module MUonEPede

using LinearAlgebra, Rotations, StaticArrays
using UnROOT, FortranFiles, EzXML
using ProgressBars

export MUonEModule, Track
export strip_to_local, local_to_global, global_to_local
export intersection, interpolate
export getmodules, generatebin, generatebinmc

include("types.jl")
include("io.jl")
include("intersection.jl")
include("mille.jl")
include("transforms.jl")

include("montecarlo.jl")

end
