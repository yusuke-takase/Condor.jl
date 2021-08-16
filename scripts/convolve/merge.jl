using NPZ
using PyCall
using Healpix

np = pyimport("numpy")

spath = "/group/cmb/litebird/usr/ytakase/beam_convolution/devided_maps/outputs/mapdata"

nside = 512
npix = nside2npix(nside)

RunRange = 1536
outmap = zeros(2, npix)
for i in 1:RunRange
    outmap .+= npzread(spath*"/idx=$i.npz")["convolved_map"]
end

np.savez_compressed("/group/cmb/litebird/usr/ytakase/beam_convolution/devided_maps/result/hiroaki_conv", outmap = outmap)