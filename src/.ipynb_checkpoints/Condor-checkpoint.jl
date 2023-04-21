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
using Falcons
#const hp = pyimport("healpy")
#np = pyimport("numpy")
#pd = pyimport("pandas")
using WignerD
using FFTW
using NPZ

#=
function __init__()
    copy!(hp, pyimport_conda("healpy", "healpy"))
end
=#

#include("./function/adding_YN_func.jl")
include("./function/tools.jl")
include("./function/convolution.jl")
include("./function/scanning.jl")
include("./function/structure_define.jl")


export get_2Ddata, symmetrizer, dBi, gen_beammap, beam_pointor, AlmPair, gen_Blm, gen_GaussBeammap
export EffectiveBeamConvolution, GaussBeamConvolution, convolvor
export gauss_3d_xyz, make_beam_TQU, unique_theta, FFTConv_demo, get_psi_make_TOD, solver_matrix, for_healpy_order, pix_calcmax, FFTConv_demo_onlyT, get_psi_make_TOD_T, unique_theta_detect
export FFTConvolution_QU, FFTConvolution_T, get_psi_make_TOD_TQU_HWP, solver_matrix_TQU
export get_pointings_theta_phi_psi_alpha_pix_tod, tod_convolution_idalhwp
export ConvolutionParams, gen_ConvolutionParams
end
