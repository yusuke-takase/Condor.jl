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