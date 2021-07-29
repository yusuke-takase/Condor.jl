for i in `seq 1 1536`
do
    bsub -q s julia get_psiDataBase.jl $i
done
