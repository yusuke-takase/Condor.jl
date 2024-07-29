mutable struct ConvolutionParams{I<:Int, MC<:Matrix{Complex{Float64}}, VI<:Vector{Int64}}
    npix::I
    lmax::I
    alm::MC
    blm::MC
    l_range::VI
end

function gen_ConvolutionParams(;
        npix=12*128^2,
        lmax=3*128-1,
        alm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        blm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        l_range = [0,0]
    )
    return ConvolutionParams(
        npix,
        lmax,
        alm,
        blm,
        l_range
    )
end

mutable struct ConvolutionParams_test_{I<:Int, MC<:Matrix{Complex{Float64}}, VI<:Vector{Int64}, VF<:Vector{Float64}}
    nside::I
    npix::I
    lmax::I
    sky_lmax::I
    beam_lmax::I
    beam_mmax::I
    alm::MC
    blm::MC
    l_range::VI
    phi::VF
    theta::VF
    psi::VF
    alpha::VF
end

function gen_test2(;
        nside = 128,
        npix=12*nside^2,
        lmax = 128*3-1,
        sky_lmax=3*nside-1,
        beam_lmax = 3*nside-1, 
        beam_mmax = 3*nside-1, 
        alm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        blm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        l_range = [0,0],
        phi = [0.0,0.0],
        theta = [0.0,0.0],
        psi = [0.0,0.0],
        alpha = [0.0,0.0]
    )
    return ConvolutionParams_test_(
        nside,
        npix,
        lmax,
        sky_lmax,
        beam_lmax,
        beam_mmax,
        alm,
        blm,
        l_range,
        phi,
        theta,
        psi,
        alpha
    )
end