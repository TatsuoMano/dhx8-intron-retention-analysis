#!/bin/bash
input_dir=/mnt/NAS/rawdata/RH23097820/rawdata
output_dir=output
mkdir -p ${output_dir}
mkdir -p output/fastqc_pre
mkdir -p output/trimmed
mkdir -p output/fastqc_post

function func(){
    . ./001_config.txt
    sample=${1}
    input_dir=/mnt/NAS/rawdata/RH23097820/rawdata
    output_dir=output

    fastqc --nogroup -t 24 \
           -o ${output_dir}/fastqc_pre \
           ${input_dir}/${sample}_1.fq.gz \
           ${input_dir}/${sample}_2.fq.gz \
        &>>output/020_fastqc_pre.log
    
    java -jar ${trimmomatic} PE -threads 24 -phred33 \
         ${input_dir}/${sample}_1.fq.gz \
         ${input_dir}/${sample}_2.fq.gz \
         ${output_dir}/trimmed/${sample}_1_trimmed.fq.gz ${output_dir}/trimmed/${sample}_1_unpaired.fq.gz \
	 ${output_dir}/trimmed/${sample}_2_trimmed.fq.gz ${output_dir}/trimmed/${sample}_2_unpaired.fq.gz \
	 ILLUMINACLIP:/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:10:10 \
         LEADING:20 \
         TRAILING:20 \
         SLIDINGWINDOW:4:15 \
         MINLEN:30 \
    &>>output/030_trim.log

    # Run fastqc post trimming
    fastqc --nogroup -t 24 \
           -o ${output_dir}/fastqc_post \
           ${output_dir}/trimmed/${sample}_1_trimmed.fq.gz \
           ${output_dir}/trimmed/${sample}_2_trimmed.fq.gz \
        &>>output/040_fastqc_post.log

}

export -f func

cat 000_sample_list.txt | xargs -P4 -I % bash -c "func %"
