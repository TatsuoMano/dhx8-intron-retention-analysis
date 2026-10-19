#!/bin/bash
function func(){
    i=${1}
    samtools sort -@ 32 output/mapping/${i}_Aligned.sortedByCoord.out.bam -o output/mapping/${i}_sorted.bam
    samtools index -@ 24 output/mapping/${i}_sorted.bam
}

export -f func
cat 000_sample_list.txt | xargs -P4 -I % bash -c "func %"
