mutable struct ConvolutionParams{I<:Int, VC<:Vector{Complex{Float64}}}
    nside::I
    lmax::I
    alm::VC
    blm::VC
end

function gen_convolutionParams(;
        nside=128,
        lmax=3nside-1,
        alm = [1.0+1im],
        blm = [1.0+1im]
    )
    return ConvolutionParams(
        nside,
        lmax,
        alm,
        blm
    )
end