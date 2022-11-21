using Condor
using DataFrames
using NPZ
using WignerD
using FFTW

include("/home/cmb/yusuket/program/scan_strategy/falcons/FP_config/fp_config03.jl")
mft_path = "/home/cmb/yusuket/program/scan_strategy/falcons/FP_config/MFTFP4Falcons210625.csv"

target_det = get_MHFT_beam_pointing("MFT", wafer=5, pixel=56)

ss = gen_ScanningStrategy()

day = 60 * 60 * 24
year = day * 365

idx=parse(Int64, ARGS[1])
nside = parse(Int64, ARGS[2])
dir_1=ARGS[3]
dir=dir_1*"test_$idx.hdf5"
Hz = parse(Int64, ARGS[4])
dir_save=dir_1*"map_$nside=$idx=$Hz"

ss.nside = nside
ss.duration = year #[sec]
ss.sampling_rate = Hz #[Hz]
ss.alpha = 45 #[degree]
ss.beta = 50 #[degree]
ss.prec_rpm = period2rpm(192.348)
ss.spin_rpm = 0.05 #[rpm]
ss.hwp_rpm = 0.0 #[rpm]
ss.start_point = "pole" #You can choose "pole" or "equator"
ss.coord="G"

ss.FP_theta = [0] #[target_det.theta[1]]
ss.FP_phi = [0] #[target_det.phi[1]] .+ 30

npix = nside2npix(ss.nside)
res = Resolution(ss.nside)
lmax = 3ss.nside-1

d_var , h= get_psi_make_TOD_T(ss, division = 1600, idx =idx, map_div=4, dir=dir)

#@show d_var , h
np.save(dir_save, d_var)
#conved = solver_matrix(d_var,h)

#np.save(dir, conved)