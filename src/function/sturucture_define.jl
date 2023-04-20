mutable struct ConvolutionParams{I<:Int, MC<:Matrix{Complex{Float64}}}
    nside::I
    lmax::I
    alm::MC
    blm::MC
end

function gen_ConvolutionParams(;
        nside=128,
        lmax=3nside-1,
        alm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im],
        blm = [1.0+1im 1.0+1im;1.0+1im 1.0+1im]
    )
    return ConvolutionParams(
        nside,
        lmax,
        alm,
        blm
    )
end