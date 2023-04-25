mutable struct bmcut{S<:String, I<:Int, F<:Float64, AAF<:AbstractArray{Float64}, AAC<:AbstractArray{Complex}}
    # eam type corresponding to Grasp "cut" file format.
    header::S
    ncomp::I # number of field components
    np::I
    ncut::I
    icomp::I
    icon::I
    sa::F
    da::F
    ca::AAF
    amp::AAC
end
function bm_cut_read(filename)
    sa, da, np, ca, icomp, icon, ncomp = 0.0, 0.0, 0, 0.0, 0, 0, 0
    ncut = 0
    header = ""
    open(filename, "r") do file
        header = readline(file)
        while !eof(file)
            line = readline(file)
            data = split(line)
            if length(data) == 7
                sa, da, np, ca, icomp, icon, ncomp = map(x -> parse(Float64, x), split(line))
                np, icomp, icon, ncomp = map(Int, (np, icomp, icon, ncomp))
                ncut += 1
            end
            if ncomp > 2
                error("Three field components present. Beam package can only handle two field components.")
            end

            if np % 2 == false
                error("The number of pixels in a cut must be odd.")
            end
        end
    end
    beam_amp = zeros(Complex, (ncomp, np, ncut))
    beam_ca = zeros(Float64, ncut)
    icut = 1
    open(filename, "r") do file
        while !eof(file)
            line = readline(file)
            data = split(line)
            if tryparse(Float64, data[1]) != nothing
                if length(data) == 7
                    sa, da, np, ca, icomp, icon, ncomp = map(x -> parse(Float64, x), split(line))
                    np, icomp, icon, ncomp = map(Int, (np, icomp, icon, ncomp))
                    beam_ca[icut] = ca
                end
                for i in 1:np
                    line = readline(file)
                    data = split(line)
                    tmp1, tmp2, tmp3, tmp4 = map(x -> parse(Float64, x), data)
                    beam_amp[1, i, icut] = Complex(tmp1, tmp2)
                    beam_amp[2, i, icut] = Complex(tmp3, tmp4)
                end
                icut += 1
            end
        end
    end
    bmcut(header,ncomp,np,ncut,icomp,icon,sa,da,beam_ca,beam_amp) 
end


function bm_cut_free(beam::bmcut)
    beam = nothing
end