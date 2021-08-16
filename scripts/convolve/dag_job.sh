out=/group/cmb/litebird/usr/ytakase/beam_convolution/devided_maps/outputs/std

bsub -q l -J 'Convolve[1-1000]' -o "$out/stdout.%J.%I" -e "$out/stderr.%J.%I" "julia effective_beam_convolution01.jl \$LSB_JOBINDEX"

bsub -q h -J 'Convolve[1001-1300]' -o "$out/stdout.%J.%I" -e "$out/stderr.%J.%I" "julia effective_beam_convolution01.jl \$LSB_JOBINDEX"

bsub -q cmb_p -J 'Convolve[1301-1536]' -o "$out/stdout.%J.%I" -e "$out/stderr.%J.%I" "julia effective_beam_convolution01.jl \$LSB_JOBINDEX"

bsub -q cmb_px -w 'done("Convolve[*]")' -J 'Mearging' julia merge.jl
