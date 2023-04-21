function get_psi_make_TOD_TQU_HWP(ss::ScanningStrategy,; division::Int, idx, map_div, dir_T, dir_QU)
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
    fid_T = h5open(dir_T*"=$idx", "r")
    fid_QU = h5open(dir_QU*"=$idx", "r")
    @views @inbounds for i = 1:division
        END = i .* month
        pix_tod, psi_tod, time_array = get_pointing_pixels(ss, BEGIN, END)
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds for k = eachindex(psi_tod[:,1])
                t = time_array[k]
                ipix = pix_tod_jth_det[k]
                if pixmin <= ipix <= pixmax
                    psi = - psi_tod_jth_det[k] .+pi
                    @views position_psi = Int(round(psi.*(calcmax_psi./(2 .*pi))))
                    
                    if position_psi == 0
                        position_psi = 1
                    end
                    #position_psi =1
                    @views re_psi = (position_psi.-1).*(2pi./calcmax_psi) .+ 2.0 .*mod2pi(ω_hwp).*t
                    #@show psi ,position_psi, re_psi 
                    #@show re_psi
                    position = 1 .+ipix .- pixmin
                    @views T_φψ=read(fid_T,"i=$ipix")
                    @views QiU_φψ=read(fid_QU,"i=$ipix")
                    @views d[1,position] += (real(T_φψ[position_psi]) +  real(QiU_φψ[position_psi].*exp(-4im .*mod2pi(ω_hwp).*t)))
                    @views d[2,position] += (real(T_φψ[position_psi]) +  real(QiU_φψ[position_psi].*exp(-4im.*mod2pi(ω_hwp).*t))).*exp(2im*re_psi)
                    @views d[3,position] += (real(T_φψ[position_psi]) +  real(QiU_φψ[position_psi].*exp(-4im.*mod2pi(ω_hwp).*t))).*exp(-2im*re_psi)
                    @views d[4,position] += 1
                    @views result_h[1,position] += exp(2 .*1im.*re_psi)
                    @views result_h[2,position] += exp(4 .*1im.*re_psi)
                end
            end
        end
        BEGIN = END
        next!(p)
    end
    close(fid_T)
    close(fid_QU)
    #return  psi_db, d
    return  [d[1,:]./Int.(d[4,:]) d[2,:]./Int.(d[4,:]) d[3,:]./Int.(d[4,:])],  [result_h[1,:]./Int.(d[4,:]) result_h[2,:]./Int.(d[4,:])]
end

function solver_matrix_TQU(dj, h)
    result_calc = zeros(ComplexF64, 3, length(h[:,1]))
    d_vector =  zeros(ComplexF64,3)
    h_matrix = zeros(ComplexF64, 3, 3)
    @inbounds for i in 1:length(h[:,1])
        h_matrix = @SMatrix [
        1.0 1.0/2*(h[i,1]) 1.0/2.0*conj(h[i,1]);
        1.0/2*(h[i,1]) 1.0/4.0*(h[i,2]) 1.0/4.0;
        1.0/2*conj(h[i,1]) 1.0/4.0 1.0/4.0*conj(h[i,2])
        ]
        d_vector = @SVector [dj[i,1], dj[i,2]/2.0, dj[i,3]/2.0]
        result_calc[:,i] .= h_matrix \ d_vector
    end
    return result_calc
end


function initial_state(ss::ScanningStrategy)
    ex = [1.0, 0.0, 0.0]
    ey = [0.0, 1.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    spin_axis = [cosd(ss.alpha), 0, sind(ss.alpha)]
    b₀ = [cosd(ss.alpha+ss.beta), 0, sind(ss.alpha+ss.beta)]
    anti_sun_axis = ex
    
    flip_angle = 0
    scan_direction = 1
    if ss.start_point == "equator"
        flip_angle = π
        scan_direction = -1
    end
    
    b₀ = vector_rotator(b₀, flip_angle, spin_axis)
    d₀ = vector_rotator(b₀, -π/2, ey)
    u₀ = scan_direction *  b₀ × d₀
    
    if ss.coord == "G"
        ex = ecliptic2galactic(ex)
        ey = ecliptic2galactic(ey)
        ez = ecliptic2galactic(ez)
        spin_axis = ecliptic2galactic(spin_axis)
        b₀ = ecliptic2galactic(b₀)
        anti_sun_axis = ex
        d₀ = ecliptic2galactic(d₀)
        u₀ = ecliptic2galactic(u₀)
    end
    
    qb₀ = Quaternion(0.0, b₀)
    qd₀ = Quaternion(0.0, d₀)
    qu₀ = Quaternion(0.0, u₀)
    qspin_axis = Quaternion(0.0, spin_axis)
    qanti_sun_axis = Quaternion(0.0, anti_sun_axis)

    q = quaternion_rotator(ss.start_angle, 1.0, ez)
    qb₀ = q * qb₀ / q
    qd₀ = q * qd₀ / q
    qu₀ = q * qu₀ / q
    qspin_axis = q * qspin_axis / q
    qanti_sun_axis = q * qanti_sun_axis / q

    spin_axis = vect(qspin_axis)
    anti_sun_axis = vect(qanti_sun_axis)
    
    return (qb₀, qd₀, qu₀, spin_axis, anti_sun_axis)
end

@inline function vec2ang_ver2(x, y, z)
    norm = sqrt(x^2 + y^2 + z^2)
    theta = acos(z / norm)
    phi = atan(y, x)
    phi = ifelse(phi > π, π-phi, phi)
    return (theta, phi)
end


@inline function quaternion_rotator(omega, t, rotate_axis)
    #= Generate a quaternion that rotates by the angle omega*t around the rotate_axis axis. =#
    rot_ang = omega * t
    Quaternion([cos(rot_ang/2.0), rotate_axis[1]*sin(rot_ang/2.0), rotate_axis[2]*sin(rot_ang/2.0), rotate_axis[3]*sin(rot_ang/2.0)])
end

@inline function vector_rotator(vec, rot_ang, rotate_axis)
    q = Quaternion([cos(rot_ang/2.0), rotate_axis[1]*sin(rot_ang/2.0), rotate_axis[2]*sin(rot_ang/2.0), rotate_axis[3]*sin(rot_ang/2.0)])
    vec_q = Quaternion(0.0, vec)
    rot_vec = q * vec_q / q
    return vect(rot_vec)
end

function ecliptic2galactic(v)
    β, λ = vec2ang(v[1], v[2], v[3])
    b, l = rot_E2G_ang(β, λ)
    x, y, z = ang2vec(b, l)
    return @SVector [x, y, z]
end

function rot_E2G_ang(β, λ)
    # reference:: https://aas.aanda.org/articles/aas/full/1998/01/ds1449/node3.html
    # equinox 1950
    # (theta,phi)
    # -π < λ,l < π  ,   0 < β < π
    β_NGP = deg2rad(29.81)
    λ_0 = deg2rad(269.32)
    l_1 = deg2rad(6.38)
    β = pi/2.0 - β # 0 < β < π →  -π/2 < β < π/2
    sin_b = (sin(β) * sin(β_NGP) - cos(β) * cos(β_NGP) * sin(λ - λ_0))
    b = asin(sin_b)

    cos_l_l_1 = (cos(λ - λ_0) * cos(β) / cos(b))
    sin_l_l_1 = ((sin(β) * cos(β_NGP) + cos(β) * sin(β_NGP) * sin(λ - λ_0)) / cos(b))
    l_l_1 = atan(sin_l_l_1, cos_l_l_1)
    l = l_l_1 + l_1
    b = pi/2.0 - b
    return b, l
end

function get_pointings_theta_phi_psi_alpha_pix_tod(ss::ScanningStrategy, start, stop)
    resol = Resolution(ss.nside)
    omega_spin = rpm2angfreq(ss.spin_rpm)
    omega_prec = rpm2angfreq(ss.prec_rpm)
    omega_revol = (2π) / (60.0 * 60.0 * 24.0 * 365.25)
    time_array = start:1/ss.sampling_rate:stop-1/ss.sampling_rate |> LinRange
    if start > stop-1/ss.sampling_rate
        error("ERROR: \n The `start` time of the calculation is greater than or equal to the `stop` time.")
    end
    loop_times = length(time_array)
    numof_det = length(ss.FP_theta)

    psi_tod = zeros(loop_times, numof_det)
    theta_tod = zeros(loop_times, numof_det)
    phi_tod = zeros(loop_times, numof_det)
    pix_tod = zeros(loop_times, numof_det)
    alpha_tod = zeros(loop_times, numof_det)
    ω_hwp = rpm2angfreq(ss.hwp_rpm)
    
    ey = @SVector [0.0, 1.0, 0.0]
    ez = @SVector [0.0, 0.0, 1.0]
    if ss.coord == "G"
        ey = ecliptic2galactic(ey)
        ez = ecliptic2galactic(ez)
    end

    qb₀, qd₀, qu₀, spin_axis, antisun_axis = initial_state(ss)
    
    @views @inbounds for j = eachindex(ss.FP_theta)
        qtheta_in_FP = quaternion_rotator(deg2rad(ss.FP_theta[j]), 1, ey)
        qphi_in_FP = quaternion_rotator(deg2rad(ss.FP_phi[j]), 1, vect(qb₀))
        qfp = qphi_in_FP * qtheta_in_FP
        qp₀ⱼ = qfp * qb₀ / qfp
        @views @inbounds @threads for i = eachindex(time_array)
            t = time_array[i]
            qᵣ = quaternion_rotator(omega_revol, t, ez)
            qₚ = quaternion_rotator(omega_prec, t, antisun_axis)
            qₛ = quaternion_rotator(omega_spin, t, spin_axis)
            Q = qᵣ * qₚ * qₛ
            qp = Q * qp₀ⱼ / Q
            qu = Q * qu₀ / Q
            
            p = vect(qp)
            u = vect(qu)
            
            ell = (p × ez) × p  
            
            θ, ϕ = vec2ang_ver2(p[1], p[2], p[3])
            θ, ϕ = vec2ang(p[1], p[2], p[3])

            theta_tod[i, j] = θ
            phi_tod[i, j] = ϕ
            pix_tod[i,j] = ang2pixRing(resol, θ,ϕ)

            k = ell × u
            cosk = dot(u, ell) / (norm(u) * norm(ell))
            cosk = ifelse(abs(cosk) > 1.0, sign(cosk), cosk)
            
            sign_kz = sign(k[3])
            sign_kz = ifelse(sign_kz==0, -1, sign_kz)
            sign_pz = sign(p[3])
            sign_pz = ifelse(sign_pz==0, 1, sign_pz)
            psi_tod[i, j] = acos(cosk) * sign_kz * sign_pz
            
            alpha_tod[i, j] = mod2pi(ω_hwp*t)
        end
    end
    return (theta_tod, phi_tod, psi_tod, alpha_tod, Int.(pix_tod))
end