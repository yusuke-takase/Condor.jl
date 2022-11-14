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
    
    
    map = PolarizedHealpixMap{Float64, RingOrder}(nside)
    map.i.pixels[:] .= UNSEEN
    map.q.pixels[:] .= UNSEEN
    map.u.pixels[:] .= UNSEEN
    resol = Resolution(256)
    
    for ipix in 1:nside2npix(nside)
        theta, phi = pix2angRing(resol, ipix)
        if main_present && theta<main.theta_max
            map.i.pixels[ipix] = bm_polar_get_value(main, theta, phi, 1)
            map.q.pixels[ipix] = bm_polar_get_value(main, theta, phi, 2)
            map.u.pixels[ipix] = bm_polar_get_value(main, theta, phi, 3)
        elseif full_present
            map.i.pixels[ipix] = bm_polar_get_value(full, theta, phi, 1)
            map.q.pixels[ipix] = bm_polar_get_value(full, theta, phi, 2)
            map.u.pixels[ipix] = bm_polar_get_value(full, theta, phi, 3)
        else
            continue
        end
    end
    
    alm = map2alm(map, lmax = lmax, mmax = mmax)
    if healpix_output
        return map, alm
    else
        return alm 
    end
end

function bm_polar_get_value(beam::bmpolar, theta, phi, x)
    dth = (pi/2.) / 1000.
    dph = 2. * pi / 1000.

    ith1 = trunc(Int, theta/dth)
    ith1 = max(1, min(1000,ith1))
    ith2 = ith1+1
    iph1 = trunc(Int, phi/dph)
    iph1 = max(1, min(1000,iph1))
    iph2 = iph1+1
    if (iph2>1000) 
        iph2=1
    end

    th1 = 0. + (ith1 - 1)*dth
    wth = 1. - (theta-th1)/dth

    ph1 = (iph1 - 1)*dph
    wph = 1. - (phi-ph1)/dph
    value = wth * (wph * beam.stokes[x, iph1, ith1] + (1. - wph) * beam.stokes[x, iph2, ith1]) + 
            (1. - wth)* (wph * beam.stokes[x, iph1, ith2] + (1. - wph) * beam.stokes[x, iph2, ith2])
    return value
end