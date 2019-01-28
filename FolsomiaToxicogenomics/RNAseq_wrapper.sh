#!/bin/bash
#### RUN RNAseq ALIGNMENTS ON FOLDER OF FASTQ FILES USING STAR, THEN HANDOFF DATA TO QoRTS, EdgeR, and DESeq2 ####

# Path to genome
genomeDir="${HOME}/dbs/folsomia/GCF_002217175.1_ASM221717v1/ensembl"
# GTF="${HOME}/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.corrected.no_gene_name.added_transcript_ids_final.gtf"
GTF="${HOME}/dbs/folsomia/GCF_002217175.1_ASM221717v1/ensembl/fcand_genes.full.gtf"
# GENOME="${HOME}/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.QoRTs_format.fna"
GENOME="${HOME}/dbs/folsomia/GCF_002217175.1_ASM221717v1/ensembl/fcand_genome.QoRTs_format.multiline.fa"
echo $genomeDir

#### RUN STAR FOR EACH FASTQ SEPARATELY ####

mkdir ./STAR_output
mkdir ./STAR_output/QC

for fastq in *fastq
do echo ${fastq}
name=${fastq%fastq}
fname=${name%.}
echo ${name}
echo ${fname}
~/programs/STAR/source/STAR --genomeDir ${genomeDir} --runThreadN 15 --readFilesIn ${fastq} --outFileNamePrefix ./STAR_output/${name} --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --limitBAMsortRAM 80000000000
java -Xmx24G -jar ~/programs/QoRTs-1.3.0/QoRTs.jar QC --maxReadLength 800 --stranded --stranded_fr_secondstrand --singleEnded --maxPhredScore 76 --generatePlots --genomeFA ${GENOME} --rawfastq ${fastq} ./STAR_output/${name}Aligned.sortedByCoord.out.bam ${GTF} ./STAR_output/QC/${fname}
done

#### RUN STAR FOR ALL FASTQs TOGETHER ####

# Initialize array for filenames
files=()
# Populate with FASTQ files
files=(*fastq)
# Convert to a string
fileString="${files[@]}"
# Add commas to string
star_input="${fileString// /,}"
echo $star_input

# Initialize array and string for read groups
numFiles=$(echo "${#files[@]}")
RGstring=$(for (( j=0; j<${numFiles}; j++ )); do echo "ID:""${files[$j]/%.fastq/}"" SM:""${files[$j]/%.fastq/}"" , " ; done)
RGstring=$(echo ${RGstring})
RGstring=$(echo ${RGstring% ,})
echo $RGstring

# Run STAR
mkdir STAR_output_all
~/programs/STAR/source/STAR --outSAMattrRGline ${RGstring} --genomeDir ${genomeDir} --runThreadN 14 --readFilesIn ${star_input} --outFileNamePrefix ./STAR_output_all/ --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts
## Option??? --sjdbGTFtagExonParentGene gene_name

QCfolder="${PWD}/STAR_output/QC/"
echo "QC directory is ${QCfolder}"

# Finish DEG analysis in R
Rscript ${HOME}/scripts/ECCC/FolsomiaToxicogenomics/QoRTs_DESeq_Report.R ${QCfolder} ${GTF}  2>&1 | tee RNASeq_Rscript_output.log

# Done
