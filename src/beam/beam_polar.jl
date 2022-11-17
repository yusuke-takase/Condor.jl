mutable struct bmpolar{I<:Int, F<:Float64,AA<:Array}
    nphi::I
    ntheta::I
    theta_min::F # in radians
    theta_max::F # in radians
    stokes::AA
end

# Initialise bmpolar type
function bm_polar_init(nphi, ntheta, theta_min, theta_max)
    stocks = zeros(Float64, (4, nphi, ntheta))
    return bmpolar(nphi, ntheta, theta_min, theta_max, stocks)
end

function bm_polar_free(beam::bmpolar)
    beam = nothing
end


function bm_polar_stokesrotate(beam::bmpolar)
    dphi = 2 * pi / beam.nphi
    dtheta = (beam.theta_max - beam.theta_min) / (beam.ntheta-1)

    for itheta in 1:beam.ntheta
        theta = beam.theta_min + (itheta-1) * dtheta
        if theta == 0.0
            continue
        end
        for iphi in 1:beam.nphi

        phi = (iphi-1) * dphi
        cos2phi = cos(2.0 * phi)
        sin2phi = sin(2.0 * phi)

        # Change co-polar basis from y to x.

        q = -beam.stokes[2, iphi, itheta]
        u = -beam.stokes[3, iphi, itheta]

        # Rotate to polar basis.

        beam.stokes[2, iphi, itheta] =  q * cos2phi + u * sin2phi
        beam.stokes[3, iphi, itheta] = -q * sin2phi + u * cos2phi

        end
    end
    return beam
end

function bm_polar_normalise(beam::bmpolar, norm)
    if norm == "unity"
        beam.stokes = beam.stokes * 0.5
    elseif norm == "four_pi"
        beam.stokes = beam.stokes * 0.5 / (4*pi)
    elseif norm == "eight_pi"
        beam.stokes = beam.stokes * 0.25 / (4*pi)
    else
        error("Unknown normalisation convention")
    end
    return beam
end
