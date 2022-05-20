#!/bin/bash
# This script takes two fastq files (paired-end) of RNA-seq data, runs fastp, STAR, Qualimap and Salmon.
# USAGE: 
# sh RNA-seq.sh <name of fastq file 1> <name of fastq file 2>
# qsub -N <samplename> -l nodes=1:ppn=3 -F "<name of fastq file 1> <name of fastq file 2>" RNA-seq.sh
# take consideration of names of fastq files to modify ${samplename}
# take consideration of analyzed species to modify genome,transcriptome index files and gtf file
# change ${output_dir} when analyzing different projects

##PBS configure
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=40G
#PBS -N for_RNA
echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
set -x  ## Log everything
START=$(date +%s.%N)


# Enter the fastq path
file_path=/public/home/changjianye/project/magang/gz/RNAseq/gz/
cd $file_path

# initialize a variable with an intuitive name to store the name of the input fastq file
fq_1=$1
fq_2=$2

# grab base of filename for naming outputs
samplename=${fq_1%_1.fq.gz*}
echo "Sample name is $samplename"       

# specify the number of cores to use
cores=2

# directory with the genome and transcriptome index files + name of the gene annotation file
 /public/software/env01/bin/STAR  \
 --runMode genomeGenerate \
 --genomeDir /public/home/changjianye/project/duck/STAR_index \
 --runThreadN 1 \
 --genomeFastaFiles     /public/home/changjianye/project/duck/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.fna \
 --sjdbGTFfile /public/home/changjianye/project/duck/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gtf

#salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i homo38_index

genome=/public/home/changjianye/project/duck/STAR_index
transcriptome=/public/home/changjianye/project/duck/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic_index
gtf= /public/home/changjianye/project/duck/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gtf

# make all of the output directories
# The -p option means no error if existing, make parent directories as needed
output_dir=/public/home/changjianye/project/magang/gz/RNAseq/gz/$samplename/
mkdir -p ${output_dir}fastp
mkdir -p ${output_dir}STAR
mkdir -p ${output_dir}salmon
mkdir -p ${output_dir}qualimap

# set up output directories
fastp_out=${output_dir}fastp/
STAR_out=${output_dir}STAR/
salmon_out=${output_dir}salmon/
qualimap_out=${output_dir}qualimap/


## set not to open a GUI when running Qualimap
unset DISPLAY

echo "Starting QC and trimming for $samplename"

## Quality control and read trimming by fastp
/public/home/changjianye/anaconda3/envs/cuttag/bin/fastp -i ${fq_1} \
	-I ${fq_2} \
	-o ${fastp_out}trimmed_${samplename}_R1.fastq.gz \
	-O ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
	--trim_poly_x \
	-h ${fastp_out}${samplename}_fastp.html \
	-j ${fastp_out}${samplename}_fastp.json

## Mapping reads to the genome by STAR
/public/software/env01/bin/STAR --runThreadN ${cores} \
	--genomeDir ${genome} \
	--readFilesIn ${fastp_out}trimmed_${samplename}_R1.fastq.gz ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix ${STAR_out}${samplename} \
	--outSAMtype BAM Unsorted

## Sort BAM files
/public/software/env01/bin/samtools sort -@ 4 -O bam -o ${STAR_out}${samplename}_sorted.bam ${STAR_out}${samplename}Aligned.out.bam
rm ${STAR_out}${samplename}Aligned.out.bam


## Quality control by Qualimap
if [! -d ${qualimap_out}${samplename} ];then
	mkdir ${qualimap_out}${samplename}
fi

/public/home/changjianye/anaconda3/bin/qualimap rnaseq \
	-bam ${STAR_out}${samplename}_sorted.bam \
	-gtf ${gtf} \
	-outdir ${qualimap_out}${samplename} \
	-p non-strand-specific \
	-pe --java-mem-size=64G

echo "Starting Salmon run for $samplename"


## Quantifying transcript abundance by salmon
if [! -d ${salmon_out}${samplename} ];then
	mkdir ${salmon_out}${samplename}
fi

/public/home/changjianye/anaconda3/bin/salmon quant -i ${transcriptome} -l A \
	-1 ${fastp_out}trimmed_${samplename}_R1.fastq.gz \
	-2 ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
	--gcBias --validateMappings --seqBias -p ${cores} -o ${salmon_out}${samplename}


END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration
