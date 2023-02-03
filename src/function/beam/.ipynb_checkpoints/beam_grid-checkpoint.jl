mutable struct bmgrid{S<:String, I<:Int, F<:Float64,AA<:Array}
    """Create empty variables or lists of attributes for holding data for each dataset"""
    # Text Header
    header::S

    # File Type parameters
    ncomp::I  # number of field components
    nx::I
    ny::I
    kcomp::I # type of field components
    kgrid::I # type of coordinates
    ix::I
    iy::I
    xs::F
    ys::F
    xe::F
    ye::F
    
    # List of field objects
    amp::AA
end

function bm_grid_read(filename)
    header = ""
    ktype, nset, klimit = 0, 0, 0
    is, in = 0, 0
    tmp1, tmp2, tmp3, tmp4 = 0., 0., 0., 0.

    nset, kcomp, ncomp, kgrid = 0, 0, 0, 0
    ix, iy = 0, 0
    xs, ys, xe, ye = 0., 0., 0., 0.
    nx, ny = 0, 0
    
    open(filename, "r") do fi
        # Loop over initial lines before "++++" getting text
    
        while ! eof(fi)
            line = readline(fi)
            if line[1:4] == "++++"
                break
            else
                header = string(header, "\n", line )
            end
        end
        ktype = parse(Int, readline(fi))
        @assert ktype==1 "Unknown Grasp grid format, ktype != 1"
        line = split(readline(fi))
        nset = parse(Int, line[1])
        kcomp = parse(Int, line[2])
        ncomp = parse(Int, line[3])
        kgrid = parse(Int, line[4])
        if nset > 1
            println("Warning: nset > 1, only reading first beam in file")
        end
        line = split(readline(fi))
        ix = parse(Int, line[1])
        iy = parse(Int, line[2])
        i=2
        while i<=nset
            fi.readline()
        end
        line = split(readline(fi))
        xs = parse(Float64, line[1])
        ys = parse(Float64, line[2])
        xe = parse(Float64, line[3])
        ye = parse(Float64, line[4])  
        
        line = split(readline(fi))
        nx = parse(Int, line[1])
        ny = parse(Int, line[2])
        klimit = parse(Int, line[3])

        amp = zeros(Complex, (ncomp, nx, ny))
        for j in 1:ny
            if (klimit ==1)
                line = split(readline(fi))
                #TODO not sure what klimit is
                is = parse(Int, line[1])
                in = parse(Int, line[2])
            else
                is = 1
                in = nx 
            end
            for i in range(is, in)
                line = split(readline(fi))
                amp[1, i, j] = parse(Float64, line[1]) + parse(Float64, line[2]) * 1im
                amp[2, i, j] = parse(Float64, line[3]) + parse(Float64, line[4]) * 1im
            end
        end
        # return header, ncomp, nx, ny, kcomp, kgrid, ix, iy, xs, ys, xe, ye, amp
        beam = bmgrid(header, ncomp, nx, ny, kcomp, kgrid, ix, iy, xs, ys, xe, ye, amp)
        return beam
    end
end

function bm_grid_free(beam::bmgrid)
    beam = nothing
end