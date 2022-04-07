module Condor

#using Falcons
using Healpix
using PyCall
using Base.Threads
using StaticArrays
using ReferenceFrameRotations
using LinearAlgebra
using ProgressMeter
using HDF5
using Interpolations
const hp = pyimport("healpy")
np = pyimport("numpy")
pd = pyimport("pandas")
using WignerD

function __init__()
    copy!(hp, pyimport_conda("healpy", "healpy"))
end

include("./function/beam_func.jl")
include("./function/convolution_func.jl")
include("./function/adding_YN_func.jl")

export get_2Ddata, symmetrizer, dBi, gen_beammap, beam_pointor, AlmPair, gen_Blm, gen_GaussBeammap
export EffectiveBeamConvolution, GaussBeamConvolution, convolvor

end
