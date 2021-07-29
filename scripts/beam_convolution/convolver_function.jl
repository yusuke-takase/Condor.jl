function get_2Ddata(data, Nx, Ny)
    beam2d = reshape(data["copol"], (Ny, Nx))'
    theta = deg2rad.(reshape(data["theta"], (Ny, Nx)))'
    phi = deg2rad.(reshape(data["phi"], (Ny, Nx)))'
    x = rad2deg.(theta .* cos.( phi ))
    y = rad2deg.(theta .* sin.( phi ))
    return (beam2d, x, y)
end

function symmetrizer(psi_in_ipix, beam2d, Nx, Ny)
    beam_reshape = [beam2d'[1:Nx+1,:] rot180(beam2d'[Nx+1:end,:])]
    
    SIZE = size(beam_reshape)
    beam_sym = zeros(SIZE[1], SIZE[2])
    n = length(psi_in_ipix)
    
    for i in 1:n
        rot_i = Int(round(rad2deg(2psi_in_ipix[i])))
        #println(rot_i)
        beam_sym .+= circshift(beam_reshape, (0, rot_i))
    end
    
    beam_sym = beam_sym ./ n
    beam_reconst_sym = [beam_sym[:,1:Nx] ;reverse(rot180(beam_sym[1:end-1,Nx+1:Ny-1]), dims=2)]'
    
    return beam_reconst_sym
end

@. dBi(x) =10log10(x)

function gen_beammap(res::Resolution, beam2d, target_pixel)
    #beam_reshape = rotl90([beam2d'[1:361,:] rotl90(beam2d'[361:end,:],2)],2)
    beam_reshape = rot180([ beam2d'[1:Nx+1,:] reverse(rot180(beam2d'[Nx+1:end,:]), dims=2) ])
    itp = interpolate(beam_reshape, BSpline(Linear()))
    
    itp = extrapolate(itp, Flat())
    
    target_ang = pix2angRing(res, target_pixel)
    ipix = beam_pointor(res, target_ang[1], target_ang[2])
    
    beam_map = zeros(res.numOfPixels)
    @inbounds @threads for i in 1:res.numOfPixels
        ang = @views rad2deg.(pix2angRing(res, ipix[i]))
        beam_map[i] = itp(ang[1]*2, ang[2]*2)
    end
    return beam_map
end

function quaternion(theta, rotate_axis)
    V = @views @SVector [cos(theta/2), rotate_axis[1]*sin(theta/2), rotate_axis[2]*sin(theta/2), rotate_axis[3]*sin(theta/2)]
    return Quaternion(V)
end

function beam_pointor(res::Resolution, theta, phi)
    transposed_ipix = zeros(Int64, res.numOfPixels)
    qy = @views quaternion(theta, [0,-1,0])
    qz = @views quaternion(phi, [0,0,-1])
    @inbounds @threads for i in 1:res.numOfPixels
        exang = pix2angRing(res, i)
        vec = @views ang2vec(exang[1], exang[2])
        q_vec = Quaternion(0, vec)
        #Q = qy*qz
        P = qy*qz * q_vec / qy/qz
        ang = @views vec2ang(P.q1, P.q2, P.q3)
        transposed_ipix[i] = ang2pixRing(res, ang[1], ang[2])
    end
    return transposed_ipix
end


mutable struct AlmPair{T<:Matrix{ComplexF64}}
    alm::T
    almconj::T
end

function gen_Blm(res::Resolution, beammap)
    lmax = 3res.nside - 1
    beammap_pol = zeros(3, res.numOfPixels)
    beammap_pol[2,:] .= beammap .* (res.numOfPixels./ (4π*sum(beammap)))
    Blm = hp.map2alm(beammap_pol, lmax = lmax, datapath = "/home/cmb/yusuket/program/MapData/healpy-data/.", use_pixel_weights=true)
    return AlmPair(Blm, conj.(Blm))
end

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

function EffectiveBeamConvolution(res::Resolution, sky::AlmPair, beam2d, idx)
    Nx = 360
    Ny = 721
    dir = "/group/cmb/litebird/usr/ytakase/psi_db/3_30/"
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
        effective_beam = symmetrizer(psi_in_ipix, beam2d, Nx, Ny)
        beammap = gen_beammap(res, effective_beam, pixid)       
        Blm = gen_Blm(res, beammap)
        pol_signal = convolvor(res, sky, Blm)
        convolved_map[:,pixid] = pol_signal
        next!(p)
    end
    return convolved_map
end