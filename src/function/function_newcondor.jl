function gauss_3d_xyz(NPIX, fwhm, xpeak, ypeak, zpeak)
    sigma = fwhm/(2*sqrt(2*np.log(2))) 
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

function rotater(lmax, θ, φ, ψ, Blm)
    blm_temp = zeros(ComplexF64, length(Blm))
    for l in 0:lmax
        Blm_l = make_order_alm_2(Blm, lmax, l, l)
        for m in 0:l
            index_result =  for_healpy_order(l, m, lmax)
            for n in -l:l
                W = WignerD.wignerDjmn(l, m, n, φ, θ, ψ)
                blm_temp[index_result] += Blm_l[n.+l.+1]*W
            end
        end
    end
    return blm_temp
end