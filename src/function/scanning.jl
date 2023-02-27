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
