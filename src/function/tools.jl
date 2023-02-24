function gauss_3d_xyz(NPIX, fwhm, xpeak, ypeak, zpeak, res)
    sigma = fwhm/(2*sqrt(2*log(2))) 
    G = zeros(NPIX)
    for i in 1:NPIX
        ang = pix2angRing(res, i)
        x,y,z = ang2vec(ang[1], ang[2])
        r = sqrt((x - xpeak)^2 + (y - ypeak)^2 + (z - zpeak)^2)
        if r > 2.0
            r = 2.0
        end
        theta = 2*asin(Float64(r/2))
        G[i] = exp(-theta^2/(2*sigma^2))/((sigma*sqrt(2*pi)))/((sigma*sqrt(2*pi)))
    end
    return G
end

function make_beam_TQU(beam_map, NPIX, res)
    map_TQU = zeros(3, NPIX)
    map_TQU[1,:] .= beam_map
    for i in 1:NPIX
        ang = pix2angRing(res, i)
        map_TQU[2,i] = beam_map[i]*cos(2*ang[2])
        map_TQU[3,i] = beam_map[i]*sin(2*ang[2])*-1
    end
    return map_TQU
end

function make_order_alm(alm, lmax, l) #almをl固定 m: 0~lまでの順番にする
    num_l = l + 1 
    new_alm = zeros(ComplexF64, num_l)
    for m in 0:l
        idx = for_healpy_order(l, m, lmax)
        new_alm[m+1] = alm[idx]
    end
    return new_alm
end

function make_order_alm_2(alm, lmax::Integer, l::Integer, mmax::Integer) #almをl固定 m: -l~l までの順番にする
    new_alm = zeros(ComplexF64, 2mmax+1)
    for m in -mmax:mmax
        if m < 0
            idx = for_healpy_order(l, -m, lmax)
            new_alm[m + mmax + 1] = conj(alm[idx])*(-1)^m
        else
            idx = for_healpy_order(l, m, lmax)
            new_alm[m+1 + mmax] = alm[idx]
        end
    end
    return new_alm
end

function for_healpy_order(l, m::Integer, lmax::Integer)
    a = 1
    if m != 0
        for i in 1:m
            a += lmax + 2 - i
        end
    end
    a += l - m
    return a
end

function unique_theta_detect(num, NSIDE, NPIX)
    n = num-1
    a1 = 4
    d = 4
    if num < NSIDE + 1
        start = 1/2*n*(2*a1+(n-1)*d)+1
        stop  = 1/2*num*(2*a1+(num-1)*d)
    elseif num < 3NSIDE+1
        n_2 = num - NSIDE
        n = NSIDE - 1
        start = 1/2*NSIDE*(2*a1+(NSIDE-1)*d) + 4*(n_2-1)*NSIDE + 1
        stop = 1/2*NSIDE*(2*a1+(NSIDE-1)*d) + 4*(n_2)*NSIDE 
    else
        n_2 = 4NSIDE-1 - num
        start = NPIX - 1/2*(n_2+1)*(2*a1+(n_2)*d) +1
        stop = NPIX - 1/2*(n_2)*(2*a1+(n_2-1)*d)
    end
            
    return Int(start), Int(stop)
end

function pix_calcmax(idx,nside)
    calcmax=0
    calcmax_temp = 8nside
    initial_max = 6nside
    if idx < nside + 1
        a = divrem(initial_max,idx*8)
        calcmax = (a[1]+1)*idx*8
    end
    if idx > nside
        if idx < 3nside+1
            calcmax = 8nside
            #println("here")
        end
        if idx > 3nside
            a = divrem(initial_max,(4nside-idx)*8)
            calcmax = (a[1]+1)*(4nside-idx)*4*2
        end
    end
    return calcmax
end  

function　unique_theta(NPIX, res)
    θ= zeros(NPIX)
    #map_TQU[1,:] .= beam_map
    for i in 1:NPIX
        ang = pix2angRing(res, i)
        θ[i] = ang[1]
    end
    return unique(θ)
end
