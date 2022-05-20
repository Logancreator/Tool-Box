#!/bin/bash
# This script takes two fastq files (paired-end) of ATAC-seq data,runs fastp (QC and trimming),bowtie2 (mapping trimmed reads),deeptools (shifting reads and creating bigwig) and MACS2 (calling peaks).
# USAGE: 
# sh ATAC-seq.sh <name of fastq file 1> <name of fastq file 2>
# qsub -N <samplename> -l nodes=1:ppn=3 -F "<name of fastq file 1> <name of fastq file 2>" ATAC-seq.sh
# take consideration of names of fastq files to modify ${samplename}
# take consideration of analyzed species to modify genome and blacklist file
# change ${output_dir} when analyzing different projects

##PBS configure
#PBS -N for_ATAC
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=40G

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
set -x  ## Log everything
START=$(date +%s.%N)

# initialize a variable with an intuitive name to store the name of the input fastq file
fq_1=$1
fq_2=$2
# chrM name define
chrM="NC_023832.1"
# grab base of filename for naming outputs
samplename=${fq_1%.R*}
echo "Sample name is $samplename"     

# specify the number of cores to use
cores=15

# directory with the bowtie2 genome index
genome=/md01/changjy/data/Goose/Goose_bulk/ref/ncbi_ref/bowtie2_index/index
gtf=/md01/changjy/data/Goose/Goose_bulk/ref/ncbi_ref/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gtf
#blacklist=/public/home/shipy3/DB/dm6/blacklist/dm6-blacklist.v2.bed

# make all of the output directories
# The -p option means no error if existing, make parent directories as needed
output_dir=/md01/changjy/data/Goose/Goose_bulk/map/$samplename/
mkdir -p ${output_dir}fastp
mkdir -p ${output_dir}bam
mkdir -p ${output_dir}MACS2
mkdir -p ${output_dir}bigwig
mkdir -p ${output_dir}ATACFragQC
# set up output directories
fastp_out=${output_dir}fastp/
bam_out=${output_dir}bam/
peak_out=${output_dir}MACS2/
bigwig_out=${output_dir}bigwig/
ATACFragQC_out=${output_dir}ATACFragQC/
## Activate conda python2.7 environment for running deeptools
source /md01/changjy/software/miniconda2/bin/activate atac

echo "Starting QC and trimming for $samplename"
## Quality control and read trimming by fastp
~/software/miniconda2/envs/atac/bin/fastp -i ${fq_1} \
	-I ${fq_2} \
	-o ${fastp_out}trimmed_${samplename}_R1.fastq.gz \
	-O ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
	-h ${fastp_out}${samplename}_fastp.html \
	-j ${fastp_out}${samplename}_fastp.json

## Map reads to reference genome by bowtie2
echo "Starting mapping for $samplename"
~/software/miniconda2/envs/atac/bin/bowtie2 --very-sensitive -X 2000 -p ${cores} \
	-x ${genome} \
	-1 ${fastp_out}trimmed_${samplename}_R1.fastq.gz \
	-2 ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
	-S ${bam_out}${samplename}.sam &> ${bam_out}${samplename}_bowtie2_summary.txt

## Compress sam files to bam files 
~/software/miniconda2/envs/atac/bin/samtools view -S -b ${bam_out}${samplename}.sam -o ${bam_out}${samplename}.bam 
rm ${bam_out}${samplename}.sam

## Remove the mitochondrial reads after alignment
echo "Remove the mitochondrial reads of $samplename"
chrM="NC_023832.1"
echo "chrM : $chrM"
~/software/miniconda2/envs/atac/bin/samtools view -h ${bam_out}${samplename}.bam \
	| grep -v $chrM \
	| ~/software/miniconda2/envs/atac/bin/samtools sort -@ 10 -O bam -o ${bam_out}${samplename}_sorted_rmChrM.bam

#source /md01/changjy/software/miniconda2/bin/activate
#conda activate atac
## Remove duplicates by using picard
echo "Remove duplicates of $samplename"
~/software/miniconda2/envs/atac/bin/picard MarkDuplicates\
	I=${bam_out}${samplename}_sorted_rmChrM.bam \
	O=${bam_out}${samplename}_sorted_rmChrM_rmDup.bam \
	M=${bam_out}${samplename}_rmDup_metrics.txt \
	REMOVE_DUPLICATES=true

## Sorting BAM files by genomic coordinates
mv ${bam_out}${samplename}_sorted_rmChrM_rmDup.bam ${bam_out}${samplename}_rmChrM_rmDup.bam 
~/software/miniconda2/envs/atac/bin/samtools sort -@ 10 -O bam -o ${bam_out}${samplename}_sorted_rmChrM_rmDup.bam ${bam_out}${samplename}_rmChrM_rmDup.bam
rm ${bam_out}${samplename}_rmChrM_rmDup.bam

## Filter and keep the uniquely mapped reads
# filtered out multimappers by specifying `[XS] == null`
# filtered out unmapped reads by specifying in the filter `not unmapped`
echo "Keep the uniquely mapped reads of $samplename"
~/software/miniconda2/envs/atac/bin/sambamba view -h -t ${cores} -f bam -F "[XS] == null and not unmapped" \
	${bam_out}${samplename}_sorted_rmChrM_rmDup.bam > ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped.bam

## index BAM files
cd ${bam_out}
~/software/miniconda2/envs/atac/bin/samtools index ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped.bam

## Filter out reads in Blacklist Regions
echo "Filter out reads in Blacklist Regions for $samplename and Rename"
mv ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped.bam \
${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl.bam

## index BAM files
~/software/miniconda2/envs/atac/bin/samtools index ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl.bam


##ATACFragQC
echo "QC for ${bam_out}${samplename}_sorted_rmChrM.bam start"
~/software/miniconda2/bin/ATACFragQC -i ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl.bam -r $gtf -n 0 -p -o
echo "QC sucessfully"


## shift BAM files by deeptools alignmentSieve
~/software/miniconda2/envs/atac/bin/alignmentSieve -b ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl.bam \
	-p ${cores} --ATACshift \
	-o ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam


## generate bigwig files
mv ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam ${bam_out}${samplename}_rmChrM_rmDup_mapped_rmbl_shift.bam
~/software/miniconda2/envs/atac/bin/samtools sort -@ 4 -O bam -o ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam \
	${bam_out}${samplename}_rmChrM_rmDup_mapped_rmbl_shift.bam
rm ${bam_out}${samplename}_rmChrM_rmDup_mapped_rmbl_shift.bam

~/software/miniconda2/envs/atac/bin/samtools index ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam

~/software/miniconda2/envs/atac/bin/bamCoverage --bam ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam \
	--numberOfProcessors ${cores} --binSize 10 \
	--normalizeUsing RPKM \
	-o ${bigwig_out}${samplename}_RPKM_normalized.bw

## Peak calling by 

##python effective_genome_size.py  \
##/md01/changjy/data/Goose/Goose_bulk/ref/ncbi_ref/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.fna
genome_size=1083772221
~/software/miniconda2/envs/atac/bin/macs2 callpeak \
	-t ${bam_out}${samplename}_sorted_rmChrM_rmDup_mapped_rmbl_shift.bam \
	-f BAMPE -g $genome_size --keep-dup all \
	-n ${samplename} \
	-q 0.01 \
	--outdir ${peak_out} \
	&>${peak_out}${samplename}_MACS2Peaks_summary.txt


END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed!"
echo $Duration

