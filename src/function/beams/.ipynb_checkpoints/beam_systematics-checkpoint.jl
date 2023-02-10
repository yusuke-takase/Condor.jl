function symmetrize(beam::bmpolar, atack_angles, average::Bool)
    sym_beam = deepcopy(beam)
    sym_beam.stokes = zeros(size(beam.stokes))
    nhits = length(atack_angles)
    phi_resol = 2π / beam.nphi
    if average == true
        for j in 1:4
            sym_beam.stokes[j,:,:] = mean(beam.stokes[j,:,:], dims=1) .* ones(beam.nphi)
        end
    else
        for j in 1:4
            for i in eachindex(atack_angles)
                rotidx = atack_angles[i]/phi_resol |> round |> Int
                sym_beam.stokes[j,:,:] .+= circshift(beam.stokes[j,:,:], (rotidx, 0))
            end
        end
    end
    sym_beam.stokes = sym_beam.stokes ./ nhits
    return sym_beam
end

angular_gaussbeam(θ, σ) = @. exp(-θ^2/(2*σ^2))/(σ*sqrt(2*pi))/(σ*sqrt(2*pi))
fwhm2std(fwhm) = fwhm/(2*sqrt(2*log(2))) 
absdbi(x) = @. 10log10(abs(x))
dbi(x) = @. 10log10(x)
