using Condor
using NPZ
using Healpix

idx=parse(Int64, ARGS[1])
nside = parse(Int64, ARGS[2])
dir_1 = ARGS[3]
dir=dir_1*"test_$idx.hdf5"
alm_path=""
blm_path=""
alm_sky = npzread("test_alm.npy")
blm = npzread("test_blm.npy")

npix = nside2npix(nside)
res = Resolution(nside)
lmax = 3nside-1

unique_θ = unique_theta(npix, res);

FFTConv_demo_onlyT(alm_sky, blm, unique_θ, lmax, nside, idx, dir)