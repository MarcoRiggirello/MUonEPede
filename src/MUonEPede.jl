# Marco Riggirello

"""
CMS/MUonE BT data preprocessing for alignment with Millepede II
"""
module MUonEPede

import LinearAlgebra: ⋅, Diagonal
using Rotations, StaticArrays
using UnROOT, FortranFiles, EzXML, DelimitedFiles
using Optim, LineSearches
using PythonCall
using ProgressBars

export MUonEModule, MUonEStation, Stub, StubSet, Track
export generatebin, residuals

include("types.jl")

include("io.jl")
include("intersection.jl")
include("fit.jl")
include("mille.jl")
include("transforms.jl")

include("montecarlo.jl")

end
