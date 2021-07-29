function convolvor(res::Resolution, sky::AlmPair, Beam::AlmPair)
    lmax = 3res.nside - 1 
    E = 2 
    B = 3
    """ First term: m=>0, Second term: m<0 """
    AlmᴱBlmᴱ = sum(sky.alm[E,:] .* Beam.almconj[E,:]) + sum(sky.almconj[E,lmax+2:end] .* Beam.alm[E,lmax+2:end])
    AlmᴮBlmᴮ = sum(sky.alm[B,:] .* Beam.almconj[B,:]) + sum(sky.almconj[B,lmax+2:end] .* Beam.alm[B,lmax+2:end])
    AlmᴱBlmᴮ = sum(sky.alm[E,:] .* Beam.almconj[B,:]) + sum(sky.almconj[E,lmax+2:end] .* Beam.alm[B,lmax+2:end])
    AlmᴮBlmᴱ = sum(sky.alm[B,:] .* Beam.almconj[E,:]) + sum(sky.almconj[B,lmax+2:end] .* Beam.alm[E,lmax+2:end])
    return [real(AlmᴱBlmᴱ + AlmᴮBlmᴮ), real(-AlmᴱBlmᴮ + AlmᴮBlmᴱ)]
end

function EffectiveBeamConvolution(res::Resolution, sky::AlmPair, beam2d, idx, dir)
    #dir = "/group/cmb/litebird/usr/ytakase/psi_db/3_30/"
    file = h5open(dir * "idx=$idx.h5", "r")
    ariaOfPixels = length(read(file,"psi"))
    convolved_map = zeros(2, res.numOfPixels)
    run_idx_range = read(file,"map_div") * 512
    
    println("ariaOfPixels : $ariaOfPixels")
    println("Run range of idx : 1 to $run_idx_range")
    p = Progress(ariaOfPixels)
    @time for i in 1:ariaOfPixels
        pixid = ariaOfPixels * (idx - 1) + i
        psi_in_ipix = read(file,"psi/$i")
        effective_beam = symmetrizer(psi_in_ipix, beam2d)
        beammap = gen_beammap(res, effective_beam, pixid)       
        Blm = gen_Blm(res, beammap)
        pol_signal = convolvor(res, sky, Blm)
        convolved_map[:,pixid] = pol_signal
        next!(p)
    end
    return convolved_map
end