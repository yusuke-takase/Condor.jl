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
                itheta2 = ntheta + 1 - itheta
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
                itheta2 = ntheta + 1 - itheta
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


function bm_cut2polar(beam_c::bmcut,; copol="x")
    ncut, np, nphi, ntheta, icut, iphi, itheta, sign = 0, 0, 0, 0, 0, 0, 0, 0
    theta_min, theta_max, dphi = 0.0, 0.0, 0.0
    modc2, modx2 = 0.0, 0.0
    c, x, acaxs = 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im
    amp_tmp = Array{Complex{Float64}}[]

    # Check metadata of input beam.
    if beam_c.icomp != 3
        error("Error in bm_cut2polar: beam is not in linear co and cx components")
    end
    if beam_c.icon != 1
        error("Error in bm_cut2polar: beam is not in phi cuts")
    end
    if beam_c.ncomp != 2
        error("Error in bm_cut2polar: beam has the wrong number of components")
    end

    # Work out the number of theta and phi values in the cuts.
    ncut = beam_c.ncut
    np = beam_c.np

    nphi = 2 * ncut
    ntheta = div(np + 1, 2)

    # Note theta min and max angles are in radians.
    theta_min = 0.0
    theta_max = abs(beam_c.sa) * pi / 180.0

    beam_p = bm_polar_init(nphi, ntheta, theta_min, theta_max)

    # Reshape amplitude array into a temporary array.
    amp_tmp = zeros(Complex{Float64}, 2, nphi, ntheta)

    for icut = 1:ncut
        amp_tmp[:, icut, :] = beam_c.amp[:, ntheta:np, icut]
        amp_tmp[:, ncut + icut, :] = beam_c.amp[:, ntheta:-1:1, icut]
    end

    # Convert amplitudes into Stokes parameters.
    dphi = 2π / nphi

    if copol == "x"
        sign = -1
    elseif copol == "y"
        sign = 1
    else
        error("Error in bm_cut2polar: unknown value for copol")
    end

    for itheta = 1:ntheta
        for iphi = 1:nphi
            c = amp_tmp[1, iphi, itheta]
            x = amp_tmp[2, iphi, itheta]

            modc2 = abs(c)^2
            modx2 = abs(x)^2
            acaxs = c * conj(x)

            beam_p.stokes[1, iphi, itheta] = modc2 + modx2
            beam_p.stokes[2, iphi, itheta] = sign * (modc2 - modx2)
            beam_p.stokes[3, iphi, itheta] = sign * 2.0 * real(acaxs)
            beam_p.stokes[4, iphi, itheta] = 2.0 * imag(acaxs)
        end
    end
    return beam_p
end