using Healpix

function beam2alm(args::Dict)
    main_format = " "
    if haskey(args, "beam_main_file_polar")
        main_format = "polar"
    elseif haskey(args, "beam_main_file_square")
        main_format = "square"
    end
    
    full_format = " "
    if haskey(args, "beam_full_file")
        full_format = "polar"
    end
    
    main_present = false
    if main_format == "square" || main_format =="polar"
        main_present = true
    elseif main_format == " "
        error("Unknown format for main beam") 
    end
    
    full_present = false
    if full_format == "polar"
        full_present = true
    elseif full_format == " "
    else
        error("Unknown format for full-sky beam") 
    end
    
    if !main_present && !full_present
        error("There must be at least one input beam data file")
    end
    
    if main_present
        main = args["beam_main_file_$main_format"]
    end
    
    if full_present
        full = args["beam_full_file"]
    end
    
    lmax = args["beam_lmax"]
    mmax = args["beam_mmax"]
    
    if main_format == "square"
        ntheta = args["beam_ntheta"]
        nphi = args["beam_nphi"]
    end
    
    healpix_output = args["beam_healpix_output"]
    if healpix_output
        nside = args["beam_nside"]
    else
        nside = 3 * lmax -1 
    end
    
    
    # alm = bm_alm_init(lmax, mmax, 3)
    println("Calculating multipoles up to lmax = $lmax, mmax = $mmax")
    main = bm_polar_stokesrotate(main)
    
    beammap = PolarizedHealpixMap{Float64, RingOrder}(nside)
    beammap.i.pixels[:] .= 0
    beammap.q.pixels[:] .= 0
    beammap.u.pixels[:] .= 0
    resol = Resolution(256)
    
    for ipix in 1:nside2npix(nside)
        theta, phi = pix2angRing(resol, ipix)
        if main_present && theta<main.theta_max
            beammap.i.pixels[ipix] = bm_polar_get_value(main, theta, phi, 1)
            beammap.q.pixels[ipix] = bm_polar_get_value(main, theta, phi, 2)
            beammap.u.pixels[ipix] = bm_polar_get_value(main, theta, phi, 3)
        elseif full_present
            beammap.i.pixels[ipix] = bm_polar_get_value(full, theta, phi, 1)
            beammap.q.pixels[ipix] = bm_polar_get_value(full, theta, phi, 2)
            beammap.u.pixels[ipix] = bm_polar_get_value(full, theta, phi, 3)
        else
            continue
        end
    end
    
    alm = map2alm(beammap, lmax = lmax, mmax = mmax)
    if healpix_output
        return beammap, alm
    else
        return alm 
    end
end

function bm_polar_get_value(beam::bmpolar, theta, phi, x)
    dth = (beam.theta_max-beam.theta_min)/(beam.ntheta-1)
    dph = 2. * pi / beam.nphi

    ith1 = 1 + trunc(Int, theta/dth)
    ith1 = max(1, min(beam.ntheta-1,ith1))
    ith2 = ith1+1
    iph1 = 1 + trunc(Int, phi/dph)
    iph1 = max(1, min(beam.nphi,iph1))
    iph2 = iph1+1
    if (iph2>beam.nphi) 
        iph2=1
    end

    th1 = beam.theta_min + (ith1 - 1)*dth
    wth = 1. - (theta-th1)/dth

    ph1 = (iph1 - 1)*dph
    wph = 1. - (phi-ph1)/dph
    value = wth * (wph * beam.stokes[x, iph1, ith1] + (1. - wph) * beam.stokes[x, iph2, ith1]) + 
            (1. - wth)* (wph * beam.stokes[x, iph1, ith2] + (1. - wph) * beam.stokes[x, iph2, ith2])
    return value
end