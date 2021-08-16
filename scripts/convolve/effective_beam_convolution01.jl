using Falcons
using Healpix
#using Plots
using NPZ
#using StatsBase
using BenchmarkTools
using PyCall
using TickTock
using Base.Threads
using StaticArrays
using PyPlot
using Formatting
using ReferenceFrameRotations
using LinearAlgebra
using ProgressMeter
using HDF5
using Interpolations

hp = pyimport("healpy")
np = pyimport("numpy")
pd = pyimport("pandas")

include("../convolvor_function.jl")

data = npzread("/home/cmb/yusuket/program/beam_study/npz_beamdata/MFT_88.5GHz_000.0_166.7_xpol_v2.npz")
Nx = 360
Ny = 721
beam2d, x, y = get_2Ddata(data, Nx, Ny)

res = Resolution(512)
cmb_path = "/home/cmb/yusuket/program/MapData/Nside512/lensed_r0_512_non_smooth.fits"
skymap = hp.read_map(cmb_path, field=(0,1,2));
lmax = 3res.nside - 1 
alm = hp.map2alm(skymap, lmax = lmax, datapath = "/home/cmb/yusuket/program/MapData/healpy-data/.", use_pixel_weights=true)
skyalm = AlmPair(alm, conj.(alm))

idx = parse(Int64, ARGS[1])
convolved_map = EffectiveBeamConvolution(res, skyalm, beam2d, idx)

println("Saving...")
np.savez_compressed("/group/cmb/litebird/usr/ytakase/beam_convolution/devided_maps/outputs/mapdata/idx=$idx", convolved_map = convolved_map)
println("=====Finish=====")

