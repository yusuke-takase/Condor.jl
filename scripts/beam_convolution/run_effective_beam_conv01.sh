for i in `seq 1 1500`
do
    bsub -q l julia effective_beam_convolution01.jl $i
done

for i in `seq 1501 1536`
do
    bsub -q cmb_p julia effective_beam_convolution01.jl $i
done