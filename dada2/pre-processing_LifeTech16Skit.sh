#!/bin/bash
# A script that will take each raw FASTQ file from a LifeTech 16S Metagenomics Kit experiment and outputs separate FASTQ files for each V region amplicon.

# Specify database for mothur alignments
db=~/dbs/mothur/silva/silva.nr_v128.align

# Convert FASTQ files to FASTA
# Quality score encoding set for Ion Torrent output with Q33 flag
for fq in *.fastq
do
name=${fq%.fastq}
~/programs/fastx_toolkit/bin/fastq_to_fasta -i ${fq} -o ${name}.fasta -Q33
done

# Align FASTA to database (SEE ALTERNATE METHOD BELOW). This is done to identify reads that map to each V region.
#for fasta in *fasta
#do
#align_seqs.py -i ${fasta} -t ${db} -o $PWD/pynast_aligned_defaults/ &
#done

# This is done better with mothur, rather than pynast.
# First make sure FASTA files have no dashes in file names
# The mothur syntax interprets dashes as file delimiters so they must be removed from file names

for file in *fasta
	do echo $file
	mv "${file}" "${file//-/\.}"
done

files=(*fasta)
files="${files[@]}"
mothur_input=${files// /-}
echo $mothur_input

mothur "#align.seqs(fasta=${mothur_input}, reference=${db}, flip=t, processors=10)"

# List sequence names corresponding to each variable region
# This can be done using the custom R script
# Input: each .align.report file from mothur
# Output: an align.regions file with appended column for each V region identified
for input in *align.report
	do echo ${input}
	Rscript ~/scripts/ECCC/dada2/pre-processing_LifeTech16Skit.R ${input}
done

# Output a new FASTQ file for each V region specified
for regionFile in *.regions
	do echo ${regionFile}
	regionFileName=${regionFile/regions/}
	echo ${regionFileName}
	for variableRegion in $(echo V2 V3 V4 V6-7 V8 V9)
		do echo ${variableRegion}
		cat ${regionFile} | grep ${variableRegion} | awk '{print $2}' > ${regionFileName}${variableRegion}.readlist
	done
done


# Mothur could also be used to get sequences for each variable region
# mothur > pcr.seqs(fasta=${fasta}, start=V2_start, end=V2_end)



# Output new FASTQ file for each sample for each amplicon

# Primer sequences

#>V3_forward
#ACTGAGACACGGTCCARACT
#>V3_reverse
#GTATTACCGCGGCTGCTG
#>V67_forward
#ACAAGCGGHGGARCATGT
#>V67_reverse
#GACGTCATCCCCACCTTCC
#>V9_forward
#GTTACGACTTCACCCCAGTCA
#>V9_reverse
#GCGTCGTAGTCCGGATTGG
#>V2_forward
#GGCGSACGGGTGAGTAA
#>V2_reverse
#GCTGCCTCCCGTAGGAGT
#>V4_forward
#CCAGCAGCCGCGGTAATA
#>V4_reverse
#GGACTACCAGGGTATCTAATCCTGT
#>V8_forward
#GYTGTCGTCAGCTCGTGT
#>V8_reverse
#CGATTACTAGCGAYTCCGACTTCA


