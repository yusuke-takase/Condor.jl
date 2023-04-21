mutable struct ConvolutionParams{I<:Int, MC<:Matrix{Complex{Float64}}, VI<:Vector{Int64}}
    nside::I
    lmax::I
    alm::MC
    blm::MC
    l_range::VI
end

function gen_ConvolutionParams(;
        nside=128,
        lmax=3nside-1,
        alm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        blm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        l_range = [0,0]
    )
    return ConvolutionParams(
        nside,
        lmax,
        alm,
        blm,
        l_range
    )
end