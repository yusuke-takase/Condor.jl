include("beam_grid.jl")
include("beam_polar.jl")

function bm_grid2polar(beam::bmgrid; copol::String = "x")

    # Check metadata of input beam.
    @assert beam.ncomp == 2 "Error in bm_grid2polar: beam has wrong number of components"

    @assert beam.kgrid == 7 "Error in bm_grid2polar: beam is not on theta-phi grid"

    @assert abs(beam.xs) <= 1e-5 "Error in bm_grid2polar: phi coordinates does not start at zero"
    # Amplitudes at phi = 0 degrees are repeated at phi = 360 degrees,
    # so nphi is set to nx - 1.

    @assert abs(beam.xe - beam.xs -360.0) <= 1e-5 "Error in bm_grid2polar: phi range is not 360 degrees"

    nphi = beam.nx - 1
    ntheta = beam.ny

    # Note theta min and max angles are in radians.

    theta_min = beam.ys * pi / 180.0
    theta_max = beam.ye * pi / 180.0

    swaptheta = theta_min > theta_max
    if swaptheta
        println("Warning: swapping theta direction")
        theta_min = beam.ye * pi / 180.0
        theta_max = beam.ys * pi / 180.0
    end
    beam_p = bm_polar_init(nphi, ntheta, theta_min, theta_max)
    # dphi = 2. * pi / nphi

    if beam.kcomp == 3
        # Beam is expressed in linear co and cx components.
        if copol == "x"
            sign = -1
        elseif copol == "y"
            sign = 1
        else
            error("Error in bm_grid2polar: unknown value for copol")
        end
        for itheta in 1:ntheta
            itheta2 = itheta
            if swaptheta
                itheta2 = ntheta + 1 -itheta
            end
            for iphi in 1:nphi
                c = beam.amp[1, iphi, itheta]
                x = beam.amp[2, iphi, itheta]

                modc2 = abs(c)^2.
                modx2 = abs(x)^2.
                acaxs = c * conj(x)

                beam_p.stokes[1, iphi, itheta2] = modc2 + modx2
                beam_p.stokes[2, iphi, itheta2] = sign * (modc2 - modx2)
                beam_p.stokes[3, iphi, itheta2] = sign * 2.0 * real(acaxs)
                beam_p.stokes[4, iphi, itheta2] = 2.0 * imag(acaxs)
            end
        end
    
    elseif beam_g.kcomp == 9
        # Beam is expressed in |E| and sqrt(rhc/lhc) format.  Ignoring polarisation!

        for itheta in 1:ntheta
            itheta2 = itheta
            if swaptheta
                itheta2 = ntheta + 1 -itheta
            end
            for iphi in 1:nphi

                c = beam.amp(1, iphi, itheta)
                modc2 = abs(c)^2

                beam_p.stokes[1, iphi, itheta2] = modc2
                beam_p.stokes[2, iphi, itheta2] = 0.0
                beam_p.stokes[3, iphi, itheta2] = 0.0
                beam_p.stokes[4, iphi, itheta2] = 0.0
            end
        end
    else
        error("Error in bm_grid2square: beam is not in supported grid sub-format")
    end
    return beam_p
    end

    
