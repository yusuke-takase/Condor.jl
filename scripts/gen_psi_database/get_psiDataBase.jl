using Falcons
using Healpix
using NPZ
using PyCall
using Base.Threads
using StaticArrays
using Formatting
using LinearAlgebra
using HDF5

function save_psiDB(save_dir, idx, db)
    h5open(save_dir * "idx=$idx.h5", "w") do file
        write(file,"map_div", db[1])
        write(file,"idx", db[2])
        for i in eachindex(db[3])
            #println(i)
            write(file,"/psi/$i", db[3][i])
        end
    end
end

hp = pyimport("healpy")
np = pyimport("numpy")
pd = pyimport("pandas")
glob = pyimport("glob")
re = pyimport("re")
os = pyimport("os")

day = 60 * 60 * 24
year = day * 365

ss = gen_ScanningStrategy()
ss.nside = 512
ss.duration = 3year #[sec]
ss.sampling_rate = 19 #[Hz]
ss.FP_theta = [0]#mft_focalplane_pointings[:,3]
ss.FP_phi = [0]#mft_focalplane_pointings[:,4].+90
ss.alpha = 45 #[degree]
ss.beta = 50 #[degree]
ss.prec_rpm = period2rpm(192.348) #[sec]
ss.spin_rpm = 0.05 #[rpm]
ss.hwp_rpm = 0 #[rpm]
ss.start_point = "equator"

idx = parse(Int64, ARGS[1])

dir = "/group/cmb/litebird/usr/ytakase/psi_db/"
det_pixel = "3_30/"#"$ss.FP_theta[1]_$ss.FP_phi[1]"#将来はピクセルに置き換えるべき
save_dir = dir * det_pixel

if os.path.exists(save_dir) == false
    os.makedirs(save_dir)
end

db = get_psiDataBase(ss, division=32*10, idx=idx, map_div=3);

@show sizeof(db)

save_psiDB(save_dir, idx, db)
#np.savez_compressed(save_dir*det_pixel*"idx=$idx", psi=psi_db)

