#!/bin/bash
mkdir -p output/mapping
ulimit -n 4096
function func(){
    . ./001_config.txt
    i=${1}
    STAR --runThreadN 24 \
    	 --genomeDir ${STAR_index} \
	 --readFilesIn output/trimmed/${i}_1_trimmed.fq.gz output/trimmed/${i}_2_trimmed.fq.gz \
	 --outSAMtype BAM \
	 SortedByCoordinate \
	 --readFilesCommand zcat \
	 --outFileNamePrefix output/mapping/${i}_
}

export -f func

cat 000_sample_list.txt | xargs -P4 -I % bash -c "func %"
