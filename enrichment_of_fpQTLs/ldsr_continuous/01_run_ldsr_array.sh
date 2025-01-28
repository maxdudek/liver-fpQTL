#!/bin/bash

# array=(annotations/*/)

# for i in {1..10}; do 
#     dir=${array[i]}
#     echo "$dir"; 
# done

for n in annotations/* ; do
    arrIN=(${n//\// })
    NAME=${arrIN[-1]}

    echo "Submitting job for $NAME..."
    
    sbatch \
        --export=NAME="$NAME" \
        --job-name="ldsr_$NAME" \
        -c 1 \
        --mem=8G \
        -t 2-00:00:00 \
        -o job_out/slurm.%x.%j.out \
        -e job_out/slurm.%x.%j.out \
        ldsr_array.sh

done
