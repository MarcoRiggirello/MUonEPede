# Marco Riggirello

"""
CMS/MUonE BT data preprocessing for alignment with Millepede II
"""
module MUonEPede

#import Statistics: median
using LinearAlgebra
using Rotations, StaticArrays
using UnROOT, FortranFiles, EzXML, DelimitedFiles
using Optim, NLSolversBase, LineSearches
using PythonCall
using ProgressBars

const ROOT = pyimport("ROOT")

export MUonEModule, MUonEStation, Stub, StubSet, Track
export toymontecarlo, align

include("types.jl")

include("io.jl")
include("intersection.jl")
include("fit.jl")
include("mille.jl")
include("residuals.jl")
include("transforms.jl")

include("montecarlo.jl")

end
