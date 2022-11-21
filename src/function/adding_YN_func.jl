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

#=
function for_healpy_order(l, m::Integer, lmax::Integer)
    return Int(m.*(2 .*lmax .+ 1 .- m)/2 .+ l .+ 1)
end
=#

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

function l_convolver(alm, blm)
    return dot(blm, alm)
end

function l_rotater_pol(lmax, θ, φ, ψ, Blm_l, l)
    blm_temp = zeros(ComplexF64, 2l+1)
    #Blm_l = make_order_alm_2(Blm, lmax, l, l)
    for m in -l:l
        for n in -l:l
            W = @views WignerD.wignerDjmn(l, m, n, φ, θ, ψ)
            blm_temp[m+l+1] += Blm_l[n.+l.+1]*W
        end
    end
    return blm_temp
end

function test_l_calculation_pol(alm_E,alm_B, Blm_E, Blm_B, lmax, npix, l_calc, ψ, res)
    #ψ = 0.0
    conv = zeros(ComplexF64, npix)
    for l in 0:l_calc
        @show l
        alm_E_proccesed = @views make_order_alm_2(alm_E, lmax, l, l)
        alm_B_proccesed = @views make_order_alm_2(alm_B, lmax, l, l)
        Blm_E_proccesed = @views make_order_alm_2(Blm_E, lmax, l, l)
        Blm_B_proccesed = @views make_order_alm_2(Blm_B, lmax, l, l)
        _2alm = @views alm_E_proccesed .+ 1im*alm_B_proccesed
        _2Blm = @views Blm_E_proccesed .+ 1im*Blm_B_proccesed
        for i in 1:npix
            θ, φ = @views pix2angRing(res, i)
            rotated_blm = @views l_rotater_pol(lmax, θ, φ, ψ, _2Blm, l)
            conv[i] += @views l_convolver(_2alm, rotated_blm)
        end
    end
    return conv
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

function FFTConv_demo(alm, blm, unique_θ, lmax, nside, idx, dir)
    ini_calcmax = 12nside
    Tlmn_mtr = zeros(ComplexF64, ini_calcmax, ini_calcmax)
    QiU_φψ = zeros(ComplexF64, ini_calcmax, ini_calcmax)
    calcmax_psi = 8nside
    npix = nside2npix(nside)
    #@show θ
    for θ in idx:idx
        calcmax_phi = pix_calcmax(θ,nside)
        Tlmn_mtr = zeros(ComplexF64, ini_calcmax, ini_calcmax)
        l = 0
        m = 0
        n = 0
        _2alm = alm[2,:] .+ 1im*alm[3,:]
        _2blm = blm[2,:] .+ 1im*blm[3,:]
        Wigner_l = WignerD.wignerd(l, unique_θ[θ])
        Tlmn_mtr[1+m,1+n] += Wigner_l[1,1]*(alm[2,1] + 1im*alm[3,1])*conj(blm[2,1] + 1im*blm[3,1])
    
        @time for l in lmax:-1:1
            @show l
            m = 0
            p_alm = for_healpy_order(l, m, lmax)
            n = 0
            p_blm = for_healpy_order(l, n, lmax)
            Wigner_l = WignerD.wignerd(l, unique_θ[θ])
            Tlmn_mtr[1+m,1+n] += Wigner_l[l+1,l+1]*(alm[2,p_alm]+1im*alm[3,p_alm])*conj(blm[2,p_blm]+1im*blm[3,p_blm])
            m=0
            for n in 1:l
                p_blm = for_healpy_order(l, n, lmax)
                Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1+n]*(alm[2,p_alm]+1im*alm[3,p_alm])*conj(blm[2,p_blm]+1im*blm[3,p_blm])
                Tlmn_mtr[1+m, calcmax_psi+1-n] += Wigner_l[l+1,l+1-n]*(alm[2,p_alm]+1im*alm[3,p_alm])*(blm[2,p_blm]-1im*blm[3,p_blm])*(-1)^n
            end
            n = 0
            p_blm = for_healpy_order(l, n, lmax)
            for m in 1:l
                p_alm = for_healpy_order(l, m, lmax)
                Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1]*(alm[2,p_alm]+1im*alm[3,p_alm])*conj(blm[2,p_blm]+1im*blm[3,p_blm])
                Tlmn_mtr[calcmax_phi+1-m,1+n] += Wigner_l[l+1-m,l+1]*conj(alm[2,p_alm]-1im*alm[3,p_alm])*(-1)^(m)*conj(blm[2,p_blm]+1im*blm[3,p_blm])
            end
            for m in 1:l
                for n in 1:l
                    p_alm = for_healpy_order(l, m, lmax)
                    p_blm = for_healpy_order(l, n, lmax)
                    Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1+n]*(alm[2,p_alm]+1im*alm[3,p_alm])*conj(blm[2,p_blm]+1im*blm[3,p_blm])
                    Tlmn_mtr[calcmax_phi+1-m,1+n] += Wigner_l[l+1-m,l+1+n]*conj(alm[2,p_alm]-1im*alm[3,p_alm])*(-1)^(m)*conj(blm[2,p_blm]+1im*blm[3,p_blm])
                    Tlmn_mtr[1+m, calcmax_psi+1-n] += Wigner_l[l+1+m,l+1-n]*(alm[2,p_alm]+1im*alm[3,p_alm])*(blm[2,p_blm]-1im*blm[3,p_blm])*(-1)^n
                    Tlmn_mtr[calcmax_phi+1-m, calcmax_psi+1-n] += Wigner_l[l+1-m,l+1-n]*conj(alm[2,p_alm].-1im*alm[3,p_alm])*(-1)^(m)*(blm[2,p_blm]-1im*blm[3,p_blm])*(-1)^n
                end
            end
            Wigner_l=0
            #GC.gc()
        end
        QiU_φψ[1:calcmax_phi,1:calcmax_psi] = ifft(Tlmn_mtr[1:calcmax_phi,1:calcmax_psi].*calcmax_phi.*calcmax_psi)
    end
    
    pixmin,pixmax=unique_theta_detect(idx, nside,npix)
    fid = h5open(dir, "w")
    resol = Resolution(nside)
    calcmax = pix_calcmax(idx,nside)
    for i in pixmin:pixmax
        ang=pix2angRing(resol,i)
        position =Int(round(ang[2]/(2pi)*(calcmax)))+1
        if position <= 0
            position = 1
        end
        write(fid,"i=$i",QiU_φψ[position, 1:calcmax_psi])
    end
    close(fid)
    #return QiU_φψ[:,:]
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


function get_psi_make_TOD(ss::ScanningStrategy,; division::Int, idx, map_div, dir)
    nside = ss.nside
    resol = Resolution(ss.nside)
    npix = nside2npix(ss.nside)
    map_division = map_div*ss.nside
    pixmin,pixmax=unique_theta_detect(idx, ss.nside,npix)
    
    month = Int(ss.duration / division)
    ω_hwp = rpm2angfreq(ss.hwp_rpm)
    
    psi_db = [Float64[] for i in 1:(pixmax-pixmin+1)]
    d = zeros(ComplexF64, 3, (pixmax-pixmin+1))
    result_h = zeros(ComplexF64, 2, (pixmax-pixmin+1))

    BEGIN = 0
    p = Progress(division)
    calcmax = pix_calcmax(idx,ss.nside)
    calcmax_psi = 8ss.nside
    φposition_array = zeros(Int32, (pixmax-pixmin+1))
    fid = h5open(dir, "r")
    @views @inbounds for i = 1:division
        END = i .* month
        pix_tod, psi_tod, time_array = get_pointing_pixels(ss, BEGIN, END)
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds for k = eachindex(psi_tod[:,1])
                t = time_array[k]
                ipix = pix_tod_jth_det[k]
                #psi = 2ω_hwp*t - psi_tod_jth_det[k]
                #position_psi = Int(round(psi/(2pi)*(calcmax_phi)))+1
                if pixmin <= ipix <= pixmax
                    #push!(psi_db[pixmax+1-ipix], psi)
                    psi = - psi_tod_jth_det[k] .+pi
                    #psi = 2.0 .*ω_hwp.*t .- psi_tod_jth_det[k] .+pi
                    #=
                    if psi < 0
                        psi = 2pi .+ psi
                    end
                    =#
                    @views position_psi = Int(round(psi.*(calcmax_psi./(2 .*pi))))
                    
                    if position_psi == 0
                        position_psi = 1
                    end
                    #position_psi =1
                    @views re_psi = (position_psi.-1).*(2pi./calcmax_psi) .+ 4.0 .*mod2pi(ω_hwp).*t
                    #@show psi ,position_psi, re_psi 
                    #@show re_psi
                    position = 1 .+ipix .- pixmin
                    @views QiU_φψ=read(fid,"i=$ipix")
                    @views d[1,position] += real(QiU_φψ[position_psi].*exp(-2im.*4.0 .*mod2pi(ω_hwp).*t)).*exp(2im*re_psi)
                    @views d[2,position] += real(QiU_φψ[position_psi].*exp(-2im.*4.0 .*mod2pi(ω_hwp).*t)).*exp(-2im*re_psi)
                    @views d[3,position] += 1
                    @views result_h[1,position] += exp(2 .*1im.*re_psi)
                    @views result_h[2,position] += exp(4 .*1im.*re_psi)
                end
            end
        end
        BEGIN = END
        next!(p)
    end
    close(fid)
    #return  psi_db, d
    return  [d[1,:]./Int.(d[3,:]) d[2,:]./Int.(d[3,:])],  [result_h[1,:]./Int.(d[3,:]) result_h[2,:]./Int.(d[3,:])]
end



function get_psi_make_TOD_fullsky(ss::ScanningStrategy,; division::Int, map_div, QiU_φψ)
 
    d = [ComplexF64[] for i in 1:4*ss.nside-1]
    for i in 1:4*ss.nside-1
         @time d[i] = get_psi_make_TOD(ss, division = 32, idx =i, map_div=4, QiU_φψ=QiU_φψ)
    end
    return  d
end

function solver_matrix(dj,h)
    result_calc = zeros(ComplexF64, 2,length(h[:,1]))
    d_vector =  zeros(ComplexF64,3)
    h_matrix = zeros(ComplexF64, 3, 3)
    @inbounds for i in 1:length(h[:,1])
        h_matrix = @SMatrix [
        1.0/4.0*(h[i,2]) 1.0/4.0;
        1.0/4.0 1.0/4.0*conj(h[i,2])
        ]
        #d_vector = @SVector [dj[1,i]/2.0, dj[2,i]/2.0]
        d_vector = @SVector [dj[i,1]/2.0, dj[i,2]/2.0]
        @show h_matrix
        @show inv(h_matrix)
        result_calc[:,i] .= h_matrix \ d_vector
    end
    return result_calc
end
#=
function solver_matrix_full(h, dj)
    result_calc = zeros(ComplexF64, 3, npix)
    d_vector =  zeros(ComplexF64,3)
    h_matrix = zeros(ComplexF64, 3, 3)
    @inbounds for i in 1:npix
        h_matrix = @SMatrix [
        1.0 1.0/2*(h[1,i]) 1.0/2.0*conj(h[1,i]);
        1.0/2*(h[1,i]) 1.0/4.0*(h[2,i]) 1.0/4.0;
        1.0/2*conj(h[1,i]) 1.0/4.0 1.0/4.0*conj(h[2,i])
        ]
        d_vector = @SVector [dj[1,i], dj[2,i]/2.0, dj[3,i]/2.0]
        result_calc[:,i] .= h_matrix \ d_vector
    end
    return result_calc
end
=#
function get_psi_make_TOD_T(ss::ScanningStrategy,; division::Int, idx, map_div, dir)
    nside = ss.nside
    resol = Resolution(ss.nside)
    npix = nside2npix(ss.nside)
    map_division = map_div*ss.nside
    pixmin,pixmax=unique_theta_detect(idx, ss.nside,npix)
    
    month = Int(ss.duration / division)
    ω_hwp = rpm2angfreq(ss.hwp_rpm)
    
    psi_db = [Float64[] for i in 1:(pixmax-pixmin+1)]
    d = zeros(ComplexF64, 4, (pixmax-pixmin+1))
    result_h = zeros(ComplexF64, 2, (pixmax-pixmin+1))

    BEGIN = 0
    p = Progress(division)
    calcmax = pix_calcmax(idx,ss.nside)
    calcmax_psi = 8ss.nside
    φposition_array = zeros(Int32, (pixmax-pixmin+1))
    fid = h5open(dir, "r")
    @views @inbounds for i = 1:division
        END = i .* month
        pix_tod, psi_tod, time_array = get_pointing_pixels(ss, BEGIN, END)
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds for k = eachindex(psi_tod[:,1])
                t = time_array[k]
                ipix = pix_tod_jth_det[k]
                #psi = 2ω_hwp*t - psi_tod_jth_det[k]
                #position_psi = Int(round(psi/(2pi)*(calcmax_phi)))+1
                if pixmin <= ipix <= pixmax
                    #push!(psi_db[pixmax+1-ipix], psi)
                    #psi = - psi_tod_jth_det[k] .+pi
                    psi = 2.0 .*ω_hwp.*t .- psi_tod_jth_det[k] .+pi
                    #=
                    if psi < 0
                        psi = 2pi .+ psi
                    end
                    =#
                    @views position_psi = Int(round(psi.*(calcmax_psi./(2 .*pi))))
                    
                    if position_psi == 0
                        position_psi = 1
                    end
                    #position_psi =1
                    @views re_psi = (position_psi.-1).*(2pi./calcmax_psi) .+ 4.0 .*mod2pi(ω_hwp).*t
                    #@show psi ,position_psi, re_psi 
                    #@show re_psi
                    position = 1 .+ipix .- pixmin
                    @views T_φψ=read(fid,"i=$ipix")
                    @views d[1,position] += real(T_φψ[position_psi])
                    #@views d[2,position] += T_φψ[position_psi].*exp(-2im.*4.0 .*mod2pi(ω_hwp).*t).*exp(2im*re_psi)
                    #@views d[3,position] += T_φψ[position_psi].*exp(-2im.*4.0 .*mod2pi(ω_hwp).*t).*exp(-2im*re_psi)
                    @views d[4,position] += 1
                    @views result_h[1,position] += exp(2 .*1im.*re_psi)
                    @views result_h[2,position] += exp(4 .*1im.*re_psi)
                end
            end
        end
        BEGIN = END
        next!(p)
    end
    close(fid)
    #return  psi_db, d
    return  [d[1,:]./Int.(d[4,:]) d[2,:]./Int.(d[4,:]) d[3,:]./Int.(d[4,:]) ],  [result_h[1,:]./Int.(d[4,:]) result_h[2,:]./Int.(d[4,:])]
end

function solver_matrix_full(h, dj)
    result_calc = zeros(ComplexF64, 3, npix)
    d_vector =  zeros(ComplexF64,3)
    h_matrix = zeros(ComplexF64, 3, 3)
    @inbounds for i in 1:npix
        h_matrix = @SMatrix [
        1.0 1.0/2*(h[1,i]) 1.0/2.0*conj(h[1,i]);
        1.0/2*(h[1,i]) 1.0/4.0*(h[2,i]) 1.0/4.0;
        1.0/2*conj(h[1,i]) 1.0/4.0 1.0/4.0*conj(h[2,i])
        ]
        d_vector = @SVector [dj[i,1], dj[i,2]/2.0, dj[i,3]/2.0]
        result_calc[:,i] .= h_matrix \ d_vector
    end
    return result_calc
end

function FFTConv_demo_onlyT(alm, blm, unique_θ, lmax, nside, idx, dir)
    ini_calcmax = 12nside
    Tlmn_mtr = zeros(ComplexF64, ini_calcmax, ini_calcmax)
    T_φψ = zeros(ComplexF64, ini_calcmax, ini_calcmax)
    calcmax_psi = 8nside
    npix = nside2npix(nside)
    #@show θ
    for θ in idx:idx
        calcmax_phi = pix_calcmax(θ,nside)
        Tlmn_mtr = zeros(ComplexF64, ini_calcmax, ini_calcmax)
        l = 0
        m = 0
        n = 0
        alm0 = alm[1,:]
        blm0 = blm[1,:]
        #_2alm = alm[2,:] .+ 1im*alm[3,:]
        #_2blm = blm[2,:] .+ 1im*blm[3,:]
        Wigner_l = WignerD.wignerd(l, unique_θ[θ])
        #Tlmn_mtr[1+m,1+n] += Wigner_l[1,1]*(alm[2,1] + 1im*alm[3,1])*conj(blm[2,1] + 1im*blm[3,1])
        Tlmn_mtr[1+m,1+n] += Wigner_l[1,1]*(alm[1,1])*conj(blm[1,1] )
    
        @time for l in lmax:-1:1
            @show l
            m = 0
            p_alm = for_healpy_order(l, m, lmax)
            n = 0
            p_blm = for_healpy_order(l, n, lmax)
            Wigner_l = WignerD.wignerd(l, unique_θ[θ])
            #Tlmn_mtr[1+m,1+n] += Wigner_l[l+1,l+1]*(alm[2,p_alm]+1im*alm[3,p_alm])*conj(blm[2,p_blm]+1im*blm[3,p_blm])
            Tlmn_mtr[1+m,1+n] += Wigner_l[l+1,l+1]*(alm[1,p_alm])*conj(blm[1,p_blm])
            m=0
            for n in 1:l
                p_blm = for_healpy_order(l, n, lmax)
                #Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1+n]*(alm[2,p_alm]+1im*alm[3,p_alm])*conj(blm[2,p_blm]+1im*blm[3,p_blm])
                #Tlmn_mtr[1+m, calcmax_psi+1-n] += Wigner_l[l+1,l+1-n]*(alm[2,p_alm]+1im*alm[3,p_alm])*(blm[2,p_blm]-1im*blm[3,p_blm])*(-1)^n
                Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1+n]*(alm[1,p_alm])*conj(blm[1,p_blm])
                Tlmn_mtr[1+m, calcmax_psi+1-n] += Wigner_l[l+1,l+1-n]*(alm[1,p_alm])*(blm[1,p_blm])*(-1)^n
            end
            n = 0
            p_blm = for_healpy_order(l, n, lmax)
            for m in 1:l
                p_alm = for_healpy_order(l, m, lmax)
                #Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1]*(alm[2,p_alm]+1im*alm[3,p_alm])*conj(blm[2,p_blm]+1im*blm[3,p_blm])
                #Tlmn_mtr[calcmax_phi+1-m,1+n] += Wigner_l[l+1-m,l+1]*conj(alm[2,p_alm]-1im*alm[3,p_alm])*(-1)^(m)*conj(blm[2,p_blm]+1im*blm[3,p_blm])
                Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1]*(alm[1,p_alm])*conj(blm[1,p_blm])
                Tlmn_mtr[calcmax_phi+1-m,1+n] += Wigner_l[l+1-m,l+1]*conj(alm[1,p_alm])*(-1)^(m)*conj(blm[1,p_blm])
            end
            for m in 1:l
                for n in 1:l
                    p_alm = for_healpy_order(l, m, lmax)
                    p_blm = for_healpy_order(l, n, lmax)
                    #Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1+n]*(alm[2,p_alm]+1im*alm[3,p_alm])*conj(blm[2,p_blm]+1im*blm[3,p_blm])
                    #Tlmn_mtr[calcmax_phi+1-m,1+n] += Wigner_l[l+1-m,l+1+n]*conj(alm[2,p_alm]-1im*alm[3,p_alm])*(-1)^(m)*conj(blm[2,p_blm]+1im*blm[3,p_blm])
                    #Tlmn_mtr[1+m, calcmax_psi+1-n] += Wigner_l[l+1+m,l+1-n]*(alm[2,p_alm]+1im*alm[3,p_alm])*(blm[2,p_blm]-1im*blm[3,p_blm])*(-1)^n
                    #Tlmn_mtr[calcmax_phi+1-m, calcmax_psi+1-n] += Wigner_l[l+1-m,l+1-n]*conj(alm[2,p_alm].-1im*alm[3,p_alm])*(-1)^(m)*(blm[2,p_blm]-1im*blm[3,p_blm])*(-1)^n
                    Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1+n]*(alm[1,p_alm])*conj(blm[1,p_blm])
                    Tlmn_mtr[calcmax_phi+1-m,1+n] += Wigner_l[l+1-m,l+1+n]*conj(alm[1,p_alm])*(-1)^(m)*conj(blm[1,p_blm])
                    Tlmn_mtr[1+m, calcmax_psi+1-n] += Wigner_l[l+1+m,l+1-n]*(alm[1,p_alm])*(blm[1,p_blm])*(-1)^n
                    Tlmn_mtr[calcmax_phi+1-m, calcmax_psi+1-n] += Wigner_l[l+1-m,l+1-n]*conj(alm[1,p_alm])*(-1)^(m)*(blm[1,p_blm])*(-1)^n
                end
            end
            Wigner_l=0
            #GC.gc()
        end
        T_φψ[1:calcmax_phi,1:calcmax_psi] = ifft(Tlmn_mtr[1:calcmax_phi,1:calcmax_psi].*calcmax_phi.*calcmax_psi)
    end
    
    pixmin,pixmax=unique_theta_detect(idx, nside,npix)
    fid = h5open(dir, "w")
    resol = Resolution(nside)
    calcmax = pix_calcmax(idx,nside)
    for i in pixmin:pixmax
        ang=pix2angRing(resol,i)
        position =Int(round(ang[2]/(2pi)*(calcmax)))+1
        if position <= 0
            position = 1
        end
        write(fid,"i=$i",T_φψ[position, 1:calcmax_psi])
    end
    close(fid)
    #return QiU_φψ[:,:]
end
