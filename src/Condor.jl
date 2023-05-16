module Condor

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
export get_pointings_theta_phi_psi_alpha_pix_tod, tod_convolution_idalhwp, tod_convolution_like_mc, R_complex, 
export ConvolutionParams, gen_ConvolutionParams
using Statistics

#include("./function/beam_func.jl")
#include("./function/convolution_func.jl")
#include("./function/adding_YN_func.jl")

include("./function/beams/beam_alm.jl")
include("./function/beams/beam_grid.jl")
include("./function/beams/beam_cut.jl")
include("./function/beams/beam_polar.jl")
include("./function/beams/beam_convert.jl")
include("./function/beams/beam2alm.jl")
include("./function/beams/beam_systematics.jl")

#export get_2Ddata, symmetrizer, dbi, gen_beammap, beam_pointor, AlmPair, gen_Blm, gen_GaussBeammap
#export EffectiveBeamConvolution, GaussBeamConvolution, convolvor
export gauss_3d_xyz, make_beam_TQU, unique_theta, FFTConv_demo, get_psi_make_TOD, solver_matrix, for_healpy_order, pix_calcmax
export bmpolar, bm_grid2polar, truncate_alm, bm_alm_init, bm_grid_read, bm_cut_read, bm_polar_init, bm_polar_normalise!, beam2alm, symmetrize, fwhm2sigma, angular_gaussbeam, fwhm2std, bmcut, bm_cut2polar, absdbi, dbi

end
