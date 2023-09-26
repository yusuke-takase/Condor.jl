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
            B2 = W*((exp(-4im*alpha[i]))*_2blm.*ψ_temp)
            SB2 = S2.*conj.(B2)
            temp = real(sum(exp.(1im*ell_v*theta_tod[i]).*　SB0) .+ real(sum(exp.(1im*ell_v*theta_tod[i]).*　SB2)))
            tod[i] += @views temp
        end
    end
    return tod
end

function tod_convolution_like_bc_new(cp, theta_tod, phi_tod, psi_tod, alpha, M_r)
    #M_c = complex_Muller(M_r)
    C = complex_Muller(M_r)
    #Rotate_HWP = R_complex.(alpha)
    tod = zeros(Float64, length(theta_tod[:]))
    alm_full = make_order_alm_3_2(cp.alm, cp.lmax)
    blm_full = make_order_alm_4(cp.blm, cp.lmax)
    sqrt2=sqrt(2)
    for l in cp.l_range[1]:cp.l_range[2]
        #@show l
        W = WignerD.wignerd(l,pi/2)
        ell_v = Vector(-l:1:l)
        φ_temp = exp.(-1im*ell_v*(pi./2))
        ψ_temp = exp.(-1im*ell_v*(pi./2))
        for i in 1:length(theta_tod[:])
            e2ia = exp(2 * 1im * alpha[i])
            e2iac = exp(-2 * 1im * alpha[i])
            e4ia = exp(4 * 1im * alpha[i])
            e4iac = exp(-4 * 1im * alpha[i])
            #M_rotate = Rotate_HWP[i]*M_c*conj.(Rotate_HWP[i])
            ψ2 = 2.0.*psi_tod[i]
            φ_temp = exp.(-1im*ell_v*pi./2).*exp.(1im*ell_v*phi_tod[i])
            ψ_temp = exp.(-1im*ell_v*psi_tod[i]) .* exp.(-1im*ell_v*pi./2)
            
            ### m<0
            S0mm = W[:,1:l]*(alm_full[1,l+1,-l+cp.lmax+1:cp.lmax].*φ_temp[1:l])
            B0mm = W[:,1:l]*((C[1,1].*blm_full[1,l+1,-l+cp.lmax+1+4:cp.lmax+4] .+ C[2,1].*blm_full[3,l+1,-l+cp.lmax+1+6:cp.lmax+6]/sqrt2*e2ia .+ C[3,1].*blm_full[2,l+1,-l+cp.lmax+1+2:cp.lmax+2]/sqrt2*e2iac).*ψ_temp[1:l])
            SB0mm = S0mm.*conj.(B0mm)
            #SB0 = 0
            S2mm =  W[:,1:l]*(alm_full[2,l+1,-l+cp.lmax+1:cp.lmax].*φ_temp[1:l]) #exp(pi*im)
            B2mm = W[:,1:l]*((C[1,2].*blm_full[1,l+1,-l+cp.lmax+1+6:cp.lmax+6].*e2iac .*sqrt2 .+ C[2,2].*blm_full[2,l+1,-l+cp.lmax+1+4:cp.lmax+4] .+ C[3,2].*blm_full[3,l+1,-l+cp.lmax+1+8:cp.lmax+8].*e4iac).*ψ_temp[1:l])
            SB2mm = S2mm.*conj.(B2mm)
            S3mm =  W[:,1:l]*(alm_full[3,l+1,-l+cp.lmax+1:cp.lmax].*φ_temp[1:l]) #exp(pi*im)
            B3mm = W[:,1:l]*((C[1,3].*blm_full[1,l+1,-l+cp.lmax+1+2:cp.lmax+2].*e2ia .*sqrt2 .+ C[2,3].*blm_full[2,l+1,-l+cp.lmax+1:cp.lmax].*e4ia .+ C[3,3].*blm_full[3,l+1,-l+cp.lmax+1+4:cp.lmax+4]).*ψ_temp[1:l])
            SB3mm = S3mm.*conj.(B3mm)
            #SB2=0
            tempmm = real(sum(exp.(1im*ell_v[1:l]*theta_tod[i]).*　SB0mm)) .+ real(sum(exp.(1im*ell_v[1:l]*theta_tod[i]).*　SB2mm))/2. + real(sum(exp.(1im*ell_v[1:l]*theta_tod[i]).*　SB3mm))/2.
            
            ### m >= 0
            S0 = W[:,l+1:2l+1]*(alm_full[1,l+1,cp.lmax+1:l+cp.lmax+1].*φ_temp[l+1:2l+1])
            B0 = W[:,l+1:2l+1]*((C[1,1].*blm_full[1,l+1,cp.lmax+1+4:l+cp.lmax+1+4] .+ C[2,1].*blm_full[3,l+1,cp.lmax+1+6:l+cp.lmax+1+6]/sqrt2*e2ia .+ C[3,1].*blm_full[2,l+1,cp.lmax+1+2:l+cp.lmax+1+2]/sqrt2*e2iac).*ψ_temp[l+1:2l+1])
            SB0 = S0.*conj.(B0)
            #SB0 = 0
            S2 =  W[:,l+1:2l+1]*(alm_full[2,l+1,cp.lmax+1:l+cp.lmax+1].*φ_temp[l+1:2l+1]) #exp(pi*im)
            B2 = W[:,l+1:2l+1]*((C[1,2].*blm_full[1,l+1,cp.lmax+1+6:l+cp.lmax+1+6].*e2iac .*sqrt2 .+ C[2,2].*blm_full[2,l+1,cp.lmax+1+4:l+cp.lmax+1+4] .+ C[3,2].*blm_full[3,l+1,cp.lmax+1+8:l+cp.lmax+1+8].*e4iac).*ψ_temp[l+1:2l+1])
            SB2 = S2.*conj.(B2)
            S3 =  W[:,l+1:2l+1]*(alm_full[3,l+1,cp.lmax+1:l+cp.lmax+1].*φ_temp[l+1:2l+1]) #exp(pi*im)
            B3 = W[:,l+1:2l+1]*((C[1,3].*blm_full[1,l+1,cp.lmax+1+2:l+cp.lmax+1+2].*e2ia .*sqrt2 .+ C[2,3].*blm_full[2,l+1,cp.lmax+1:l+cp.lmax+1].*e4ia .+ C[3,3].*blm_full[3,l+1,cp.lmax+1+4:l+cp.lmax+1+4]).*ψ_temp[l+1:2l+1])
            SB3 = S3.*conj.(B3)
            #SB2=0
            temp = real(sum(exp.(1im*ell_v[l+1:2l+1]*theta_tod[i]).*　SB0)) .+ real(sum(exp.(1im*ell_v[l+1:2l+1]*theta_tod[i]).*　SB2))/2. + real(sum(exp.(1im*ell_v[l+1:2l+1]*theta_tod[i]).*　SB3))/2.
            
            tod[i] += @views temp + tempmm
        end
    end
    return tod
end


function tod_convolution_like_mc(cp, theta_tod, phi_tod, psi_tod, alpha, M_r)
    #M_c = complex_Muller(M_r)
    M_rotate = complex_Muller(M_r)
    #Rotate_HWP = R_complex.(alpha)
    tod = zeros(Float64, length(theta_tod[:]))
    alm_full = make_order_alm_3(cp.alm, cp.lmax)
    blm_full = make_order_alm_4(cp.blm, cp.lmax)
    sqrt2=sqrt(2)
    for l in cp.l_range[1]:cp.l_range[2]
        #@show l
        W = WignerD.wignerd(l,pi/2)
        ell_v = Vector(-l:1:l)
        φ_temp = exp.(-1im*ell_v*(pi./2))
        ψ_temp = exp.(-1im*ell_v*(pi./2))
        for i in 1:length(theta_tod[:])
            e2ia = exp(2 * 1im * alpha[i])
            e2iac = exp(-2 * 1im * alpha[i])
            e4ia = exp(4 * 1im * alpha[i])
            e4iac = exp(-4 * 1im * alpha[i])
            #M_rotate = Rotate_HWP[i]*M_c*conj.(Rotate_HWP[i])
            ψ2 = 2.0.*psi_tod[i]
            φ_temp = exp.(-1im*ell_v*pi./2).*exp.(1im*ell_v*phi_tod[i])
            ψ_temp = exp.(-1im*ell_v*psi_tod[i]) .* exp.(-1im*ell_v*pi./2)
            S0 = W*(alm_full[1,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp)
            B0 = W*((M_rotate[1,1].*blm_full[1,l+1,-l+cp.lmax+1+4:l+cp.lmax+1+4] .+ M_rotate[1,2].*blm_full[2,l+1,-l+cp.lmax+1+6:l+cp.lmax+1+6]/sqrt2*e2ia .+ M_rotate[1,3].*blm_full[3,l+1,-l+cp.lmax+1+2:l+cp.lmax+1+2]/sqrt2*e2iac).*ψ_temp)
            SB0 = S0.*conj.(B0)
            #SB0 = 0
            S2 =  W*(alm_full[2,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp) #exp(pi*im)
            B2 = W*((M_rotate[2,1].*blm_full[1,l+1,-l+cp.lmax+1+6:l+cp.lmax+1+6].*e2iac .*sqrt2 .+ M_rotate[2,2].*blm_full[2,l+1,-l+cp.lmax+1+4:l+cp.lmax+1+4] .+ M_rotate[2,3].*blm_full[3,l+1,-l+cp.lmax+1+8:l+cp.lmax+1+8].*e4iac).*ψ_temp)
            SB2 = S2.*conj.(B2)
            #SB2=0
            temp = real(sum(exp.(1im*ell_v*theta_tod[i]).*　SB0) .+ real(sum(exp.(1im*ell_v*theta_tod[i]).*　SB2)))
            tod[i] += @views temp
        end
    end
    return tod
end

function tod_convolution_like_mc(cp, theta_tod, phi_tod, psi_tod, alpha, M_r, beam_mmax)
    #M_c = complex_Muller(M_r)
    M_rotate = complex_Muller(M_r)
    #Rotate_HWP = R_complex.(alpha)
    tod = zeros(Float64, length(theta_tod[:]))
    alm_full = make_order_alm_3(cp.alm, cp.lmax)
    blm_full = make_order_alm_4(cp.blm, cp.lmax)
    sqrt2=sqrt(2)
    e2ia = exp.(2 * 1im .* alpha)
    e2iac = exp.(-2 * 1im .* alpha)
    e4ia = exp.(4 * 1im .* alpha)
    e4iac = exp.(-4 * 1im .* alpha)
    ell_v = Vector(-cp.lmax:1:cp.lmax)
    for l in cp.l_range[1]:cp.l_range[2]
        #@show l
        if beam_mmax > l
            @views mmax = l
        else
            @views mmax = beam_mmax
        end
        W = WignerD.wignerd(l,pi/2)
        for i in 1:length(theta_tod[:])
            @views ψ2 = 2.0.*psi_tod[i]
            @views φ_temp = exp.(-1im*ell_v[cp.lmax+1-l:cp.lmax+1+l]*pi./2).*exp.(1im*ell_v[cp.lmax+1-l:cp.lmax+1+l]*phi_tod[i])
            @views ψ_temp = exp.(-1im*ell_v[cp.lmax+1-mmax:cp.lmax+1+mmax]*psi_tod[i]) .* exp.(-1im*ell_v[cp.lmax+1-mmax:cp.lmax+1+mmax]*pi./2)
            @views S0 = W*(alm_full[1,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp)
            @views B0 = W[:,-mmax+l+1:mmax+l+1]*((M_rotate[1,1].*blm_full[1,l+1,-mmax+cp.lmax+1+4:mmax+cp.lmax+1+4] .+ M_rotate[1,2].*blm_full[2,l+1,-mmax+cp.lmax+1+6:mmax+cp.lmax+1+6]/sqrt2*e2ia[i] .+ M_rotate[1,3].*blm_full[3,l+1,-mmax+cp.lmax+1+2:mmax+cp.lmax+1+2]/sqrt2*e2iac[i]).*ψ_temp)
            @views SB0 = S0.*conj.(B0)
            #SB0 = 0
            @views S2 =  W*(alm_full[2,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp) #exp(pi*im)
            @views B2 = W[:,-mmax+l+1:mmax+l+1]*((M_rotate[2,1].*blm_full[1,l+1,-mmax+cp.lmax+1+6:mmax+cp.lmax+1+6].*e2iac[i] .*sqrt2 .+ M_rotate[2,2].*blm_full[2,l+1,-mmax+cp.lmax+1+4:mmax+cp.lmax+1+4] .+ M_rotate[2,3].*blm_full[3,l+1,-mmax+cp.lmax+1+8:mmax+cp.lmax+1+8].*e4iac[i]).*ψ_temp)
            @views SB2 = S2.*conj.(B2)
            #SB2=0
            @views temp = real(sum(exp.(1im*ell_v[cp.lmax+1-l:cp.lmax+1+l]*theta_tod[i]).*　SB0) .+ real(sum(exp.(1im*ell_v[cp.lmax+1-l:cp.lmax+1+l]*theta_tod[i]).*　SB2)))
            @views tod[i] += temp
        end
    end
    return tod
end

function tod_convolution_like_mc_new(cp, theta_tod, phi_tod, psi_tod, alpha, M_r)
    #M_c = complex_Muller(M_r)
    C = complex_Muller(M_r)
    #Rotate_HWP = R_complex.(alpha)
    tod = zeros(Float64, length(theta_tod[:]))
    alm_full = make_order_alm_3_2(cp.alm, cp.lmax)
    blm_full = make_order_alm_4(cp.blm, cp.lmax)
    sqrt2=sqrt(2)
    for l in cp.l_range[1]:cp.l_range[2]
        #@show l
        W = WignerD.wignerd(l,pi/2)
        ell_v = Vector(-l:1:l)
        φ_temp = exp.(-1im*ell_v*(pi./2))
        ψ_temp = exp.(-1im*ell_v*(pi./2))
        for i in 1:length(theta_tod[:])
            e2ia = exp(2 * 1im * alpha[i])
            e2iac = exp(-2 * 1im * alpha[i])
            e4ia = exp(4 * 1im * alpha[i])
            e4iac = exp(-4 * 1im * alpha[i])
            #M_rotate = Rotate_HWP[i]*M_c*conj.(Rotate_HWP[i])
            ψ2 = 2.0.*psi_tod[i]
            φ_temp = exp.(-1im*ell_v*pi./2).*exp.(1im*ell_v*phi_tod[i])
            ψ_temp = exp.(-1im*ell_v*psi_tod[i]) .* exp.(-1im*ell_v*pi./2)
            S0 = W*(alm_full[1,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp)
            B0 = W*((C[1,1].*blm_full[1,l+1,-l+cp.lmax+1+4:l+cp.lmax+1+4] .+ C[2,1].*blm_full[3,l+1,-l+cp.lmax+1+6:l+cp.lmax+1+6]/sqrt2*e2ia .+ C[3,1].*blm_full[2,l+1,-l+cp.lmax+1+2:l+cp.lmax+1+2]/sqrt2*e2iac).*ψ_temp)
            SB0 = S0.*conj.(B0)
            #SB0 = 0
            S2 =  W*(alm_full[2,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp) #exp(pi*im)
            B2 = W*((C[1,2].*blm_full[1,l+1,-l+cp.lmax+1+6:l+cp.lmax+1+6].*e2iac .*sqrt2 .+ C[2,2].*blm_full[2,l+1,-l+cp.lmax+1+4:l+cp.lmax+1+4] .+ C[3,2].*blm_full[3,l+1,-l+cp.lmax+1+8:l+cp.lmax+1+8].*e4iac).*ψ_temp)
            SB2 = S2.*conj.(B2)
            S3 =  W*(alm_full[3,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp) #exp(pi*im)
            B3 = W*((C[1,3].*blm_full[1,l+1,-l+cp.lmax+1+2:l+cp.lmax+1+2].*e2ia .*sqrt2 .+ C[2,3].*blm_full[2,l+1,-l+cp.lmax+1:l+cp.lmax+1].*e4ia .+ C[3,3].*blm_full[3,l+1,-l+cp.lmax+1+4:l+cp.lmax+1+4]).*ψ_temp)
            SB3 = S3.*conj.(B3)
            #SB2=0
            temp = real(sum(exp.(1im*ell_v*theta_tod[i]).*　SB0)) .+ real(sum(exp.(1im*ell_v*theta_tod[i]).*　SB2))/2. + real(sum(exp.(1im*ell_v*theta_tod[i]).*　SB3))/2.
            tod[i] += @views temp
        end
    end
    return tod
end

function tod_convolution_like_mc_new(cp, theta_tod, phi_tod, psi_tod, alpha, M_r, beam_mmax)
    #M_c = complex_Muller(M_r)
    C = complex_Muller(M_r)
    #Rotate_HWP = R_complex.(alpha)
    tod = zeros(Float64, length(theta_tod[:]))
    alm_full = make_order_alm_3_2(cp.alm, cp.lmax)
    blm_full = make_order_alm_4(cp.blm, cp.lmax)
    sqrt2=sqrt(2)
    ell_v = Vector(-cp.lmax:1:cp.lmax)
    for l in cp.l_range[1]:cp.l_range[2]
        #@show l
        if beam_mmax > l
            @views mmax = l
        else
            @views mmax = beam_mmax
        end
        W = WignerD.wignerd(l,pi/2)
        #ell_v = Vector(-l:1:l)
        φ_temp = exp.(-1im*ell_v*(pi./2))
        ψ_temp = exp.(-1im*ell_v*(pi./2))
        for i in 1:length(theta_tod[:])
            e2ia = exp(2 * 1im * alpha[i])
            e2iac = exp(-2 * 1im * alpha[i])
            e4ia = exp(4 * 1im * alpha[i])
            e4iac = exp(-4 * 1im * alpha[i])
            #M_rotate = Rotate_HWP[i]*M_c*conj.(Rotate_HWP[i])
            #ψ2 = 2.0.*psi_tod[i]
            @views φ_temp = exp.(-1im*ell_v[cp.lmax+1-l:cp.lmax+1+l]*pi./2).*exp.(1im*ell_v[cp.lmax+1-l:cp.lmax+1+l]*phi_tod[i])
            ψ_temp = exp.(-1im*ell_v[cp.lmax+1-mmax:cp.lmax+1+mmax]*psi_tod[i]) .* exp.(-1im*ell_v[cp.lmax+1-mmax:cp.lmax+1+mmax]*pi./2)
            S0 = W*(alm_full[1,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp)
            B0 = W[:,-mmax+l+1:mmax+l+1]*((C[1,1].*blm_full[1,l+1,-mmax+cp.lmax+1+4:mmax+cp.lmax+1+4] .+ C[2,1].*blm_full[3,l+1,-mmax+cp.lmax+1+6:mmax+cp.lmax+1+6]/sqrt2*e2ia .+ C[3,1].*blm_full[2,l+1,-mmax+cp.lmax+1+2:mmax+cp.lmax+1+2]/sqrt2*e2iac).*ψ_temp)
            SB0 = S0.*conj.(B0)
            #SB0 = 0
            S2 =  W*(alm_full[2,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp) #exp(pi*im)
            B2 = W[:,-mmax+l+1:mmax+l+1]*((C[1,2].*blm_full[1,l+1,-mmax+cp.lmax+1+6:mmax+cp.lmax+1+6].*e2iac .*sqrt2 .+ C[2,2].*blm_full[2,l+1,-mmax+cp.lmax+1+4:mmax+cp.lmax+1+4] .+ C[3,2].*blm_full[3,l+1,-mmax+cp.lmax+1+8:mmax+cp.lmax+1+8].*e4iac).*ψ_temp)
            SB2 = S2.*conj.(B2)
            S3 =  W*(alm_full[3,l+1,-l+cp.lmax+1:l+cp.lmax+1].*φ_temp) #exp(pi*im)
            B3 = W[:,-mmax+l+1:mmax+l+1]*((C[1,3].*blm_full[1,mmax+1,-mmax+cp.lmax+1+2:mmax+cp.lmax+1+2].*e2ia .*sqrt2 .+ C[2,3].*blm_full[2,l+1,-mmax+cp.lmax+1:mmax+cp.lmax+1].*e4ia .+ C[3,3].*blm_full[3,l+1,-mmax+cp.lmax+1+4:mmax+cp.lmax+1+4]).*ψ_temp)
            SB3 = S3.*conj.(B3)
            #SB2=0
            temp = real(sum(exp.(1im*ell_v[cp.lmax+1-l:cp.lmax+1+l]*theta_tod[i]).*　SB0)) .+ real(sum(exp.(1im*ell_v[cp.lmax+1-l:cp.lmax+1+l]*theta_tod[i]).*　SB2))/2. + real(sum(exp.(1im*ell_v[cp.lmax+1-l:cp.lmax+1+l]*theta_tod[i]).*　SB3))/2.
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

