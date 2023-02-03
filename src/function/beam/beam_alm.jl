mutable struct bmalm{I<:Int, F<:Float64,AA<:Array}
    lmax::I
    mmax::I
    nstokes::I
    a::AA
end

# Initialise a bmalm type with lmax, mmax, and nstokes Stokes
# parameters.

function bm_alm_init(lmax, mmax, nstokes)
    mmax_temp = mmax
    if lmax < 0
        error("bm_alm_init: lmax < 0")
    else if mmax < 0
        error("bm_alm_init: mmax < 0")
    else if mmax > lmax
    println("Warning from bm_alm_init: mmax > lmax: setting mmax = lmax")
    mmax_temp = lmax
    end if

    if !(nstokes == 1 || nstokes == 3 || nstokes == 4 )
    error("bm_alm_init: nstokes must be  1, 3 or 4")
    end if

    a = zeros(ComplexF32, (nstokes, lmax, mmax_temp))
    return bmalm(lmax, mmax_temp, nstokes, a)
end


# Free bmalm type.

function bm_alm_free(alm::Alm)
    alm = nothing
end

function truncate_alm(alm::Alm, lmax, mmax)
    alm_new = Alm(lmax,mmax)
    size = numberOfAlms(lmax, min(mmax,alm.mmax))
    for idx in 1:size
        m = convert(Int64, ceil(((2 * lmax + 1) - sqrt((2 * lmax + 1) ^ 2 - 8 * (idx - 1 - lmax))) / 2))
        l = idx - 1 - m * div(2 * lmax + 1 - m,2)
        alm_new.alm[idx] = alm.alm[almIndex(alm, l, m)]
    end
    return alm_new
end

function truncate_alm(alm::Vector, lmax, mmax)
    alm_t = truncate_alm(alm[1], lmax, mmax)
    alm_e = truncate_alm(alm[2], lmax, mmax)
    alm_b = truncate_alm(alm[3], lmax, mmax)
    return [alm_t, alm_e, alm_b]
end
                
@. gaussbeam(θ, σ) = exp(-θ^2/(2*σ^2))/(σ*sqrt(2*pi))/(σ*sqrt(2*pi))
sigma(fwhm) = fwhm/(2*sqrt(2*log(2))) 
