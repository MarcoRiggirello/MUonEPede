# Marco Riggirello

"""
CMS/MUonE BT data preprocessing for alignment with Millepede II
"""
module MUonEPede

using LinearAlgebra, Rotations, StaticArrays
using UnROOT, FortranFiles, EzXML
using ProgressBars

export MUonEModule, Track
export local_to_global, global_to_local
export intersection, interpolate
export getmodules, generatebin

include("io.jl")
include("intersection.jl")
include("mille.jl")
include("transformations.jl")
include("types.jl")

end
