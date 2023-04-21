function tod_convolution_idalhwp(cp, theta_tod, phi_tod, psi_tod, alpha)
    tod = zeros(Float64, length(theta_tod[:]))
    for l in cp.l_range[1]:cp.l_range[2]
        #@show l
        W = WignerD.wignerd(l,pi/2)
        alm_T = make_order_alm_2(cp.alm[1,:], cp.lmax, l, l)
        blm_T = make_order_alm_2(cp.blm[1,:], cp.lmax, l, l)
        alm_E= @views make_order_alm_2(cp.alm[2,:], cp.lmax, l, l)
        alm_B = @views make_order_alm_2(cp.alm[3,:], cp.lmax, l, l)
        Blm_E = @views make_order_alm_2(cp.blm[2,:], cp.lmax, l, l)
        Blm_B = @views make_order_alm_2(cp.blm[3,:], cp.lmax, l, l)
        _2alm = @views -(alm_E .+ 1im*alm_B)
        _2blm = @views -(Blm_E .+ 1im*Blm_B)
        ell_v = Vector(-l:1:l)
        φ_temp = exp.(-1im*ell_v*(pi./2))
        ψ_temp = exp.(-1im*ell_v*(pi./2))
        S0 = W*(alm_T.*φ_temp)
        B0 = W*(blm_T.*ψ_temp)
        SB0 = S0.*conj.(B0)
        S2 =  W*(_2alm.*φ_temp) #exp(pi*im)
        B2 = W*(_2blm.*ψ_temp)
        SB2 = S2.*conj.(B2)
        for i in 1:length(theta_tod[:])
            ψ2 = 2.0.*psi_tod[i]
            φ_temp = exp.(-1im*ell_v*pi./2).*exp.(1im*ell_v*phi_tod[i])
            ψ_temp = exp.(-1im*ell_v*psi_tod[i]) .* exp.(-1im*ell_v*pi./2)
            S0 = W*(alm_T.*φ_temp)
            B0 = W*(blm_T.*ψ_temp)
            SB0 = S0.*conj.(B0)
            S2 =  W*(_2alm.*φ_temp) #exp(pi*im)
            B2 = W*((exp(4im*alpha[i]))*_2blm.*ψ_temp)
            SB2 = S2.*conj.(B2)
            temp = real(sum(exp.(1im*ell_v*theta_tod[i]).*　SB0) .+ real(sum(exp.(1im*ell_v*theta_tod[i]).*　SB2)))
            tod[i] += @views temp
        end
    end
    return tod
end


function FFTConvolution_T(alm, blm, unique_θ, lmax, nside, idx, dir)
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
        Wigner_l = WignerD.wignerd(l, unique_θ[θ])
        Tlmn_mtr[1+m,1+n] += Wigner_l[1,1]*(alm[1,1])*conj(blm[1,1] )
    
        @time for l in lmax:-1:1
            @show l
            m = 0
            p_alm = for_healpy_order(l, m, lmax)
            n = 0
            p_blm = for_healpy_order(l, n, lmax)
            Wigner_l = WignerD.wignerd(l, unique_θ[θ])
            Tlmn_mtr[1+m,1+n] += Wigner_l[l+1,l+1]*(alm[1,p_alm])*conj(blm[1,p_blm])
            m=0
            for n in 1:l
                p_blm = for_healpy_order(l, n, lmax)
                Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1+n]*(alm[1,p_alm])*conj(blm[1,p_blm])
                Tlmn_mtr[1+m, calcmax_psi+1-n] += Wigner_l[l+1,l+1-n]*(alm[1,p_alm])*(blm[1,p_blm])*(-1)^n
            end
            n = 0
            p_blm = for_healpy_order(l, n, lmax)
            for m in 1:l
                p_alm = for_healpy_order(l, m, lmax)
                Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1]*(alm[1,p_alm])*conj(blm[1,p_blm])
                Tlmn_mtr[calcmax_phi+1-m,1+n] += Wigner_l[l+1-m,l+1]*conj(alm[1,p_alm])*(-1)^(m)*conj(blm[1,p_blm])
            end
            for m in 1:l
                for n in 1:l
                    p_alm = for_healpy_order(l, m, lmax)
                    p_blm = for_healpy_order(l, n, lmax)
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
end

function FFTConvolution_T(;alm, blm, nside, lmax, dir, idx)
    npix = nside2npix(nside)
    res = Resolution(nside)
    unique_θ = unique_theta(npix, res)
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
        Wigner_l = WignerD.wignerd(l, unique_θ[θ])
        Tlmn_mtr[1+m,1+n] += Wigner_l[1,1]*(alm[1,1])*conj(blm[1,1] )
    
        for l in lmax:-1:1
            #@show l
            m = 0
            p_alm = for_healpy_order(l, m, lmax)
            n = 0
            p_blm = for_healpy_order(l, n, lmax)
            Wigner_l = WignerD.wignerd(l, unique_θ[θ])
            Tlmn_mtr[1+m,1+n] += Wigner_l[l+1,l+1]*(alm[1,p_alm])*conj(blm[1,p_blm])
            m=0
            for n in 1:l
                p_blm = for_healpy_order(l, n, lmax)
                Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1+n]*(alm[1,p_alm])*conj(blm[1,p_blm])
                Tlmn_mtr[1+m, calcmax_psi+1-n] += Wigner_l[l+1,l+1-n]*(alm[1,p_alm])*(blm[1,p_blm])*(-1)^n
            end
            n = 0
            p_blm = for_healpy_order(l, n, lmax)
            for m in 1:l
                p_alm = for_healpy_order(l, m, lmax)
                Tlmn_mtr[1+m,1+n] += Wigner_l[l+1+m,l+1]*(alm[1,p_alm])*conj(blm[1,p_blm])
                Tlmn_mtr[calcmax_phi+1-m,1+n] += Wigner_l[l+1-m,l+1]*conj(alm[1,p_alm])*(-1)^(m)*conj(blm[1,p_blm])
            end
            for m in 1:l
                for n in 1:l
                    p_alm = for_healpy_order(l, m, lmax)
                    p_blm = for_healpy_order(l, n, lmax)
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
    fid = h5open(dir*"=$idx", "w")
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
end

function FFTConvolution_QU(alm, blm, unique_θ, lmax, nside, idx, dir)
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
            #@show l
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

function FFTConvolution_QU(;alm, blm, nside, lmax, dir, idx)
    npix = nside2npix(nside)
    res = Resolution(nside)
    unique_θ = unique_theta(npix, res)
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
    
        for l in lmax:-1:1
            #@show l
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
    fid = h5open(dir*"=$idx", "w")
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

