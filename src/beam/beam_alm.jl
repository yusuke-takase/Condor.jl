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

function bm_alm_free(alm)
    alm = nothing
end bm_alm_free(alm)