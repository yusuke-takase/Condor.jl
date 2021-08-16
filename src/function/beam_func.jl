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
    Nx, Ny = size(beam2d)
    beam_reshape = rot180([ beam2d'[1:Nx+1,:] reverse(rot180(beam2d'[Nx+1:end,:]), dims=2) ])
    itp = Interpolations.interpolate(beam_reshape, BSpline(Linear()))
    
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
        q_vec = Quaternion(0, vec[1], vec[2], vec[3])
        #Q = qy*qz
        P = qy*qz * q_vec / qz/qy
        ang = @views vec2ang(P.q1, P.q2, P.q3)
        transposed_ipix[i] = ang2pixRing(res, ang[1], ang[2])
    end
    return transposed_ipix
end


mutable struct AlmPair{T<:Matrix{ComplexF64}}
    alm::T
    almconj::T
end

function gen_Blm(res::Resolution, beammap, weight_path)
    lmax = 3res.nside - 1
    beammap_pol = zeros(3, res.numOfPixels)
    beammap_pol[2,:] .= beammap .* (res.numOfPixels ./ (4π*sum(beammap)))
    Blm = hp.map2alm(beammap_pol, lmax = lmax, datapath = weight_path, use_pixel_weights=true)
    return AlmPair(Blm, conj.(Blm))
end

function gen_GaussBeammap(base, fwhm, peak_vec)
    sigma = fwhm/(2.0*sqrt(2.0*log(2.0))) 
    r = @. sqrt((base[1] - peak_vec[1])^2 + (base[2] - peak_vec[2])^2 + (base[3] - peak_vec[3])^2)
    theta = zeros(length(r))
    @inbounds @threads for i in eachindex(r)
        r[i] = ifelse(r[i]>2.0, 2.0, r[i])
        theta[i] = 2.0*asin(r[i]/2.0)
    end
    return @. exp(-theta^2/(2.0*sigma^2))/((sigma*sqrt(2*π)))/((sigma*sqrt(2*π)))
end
