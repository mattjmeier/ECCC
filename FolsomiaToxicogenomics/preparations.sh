# Make GTF file compatible with STAR
# 1) Use gffread to convert from GFF to GTF
~/programs/cufflinks-2.2.1.Linux_x86_64/gffread -T GCF_002217175.1_ASM221717v1_genomic.gff -o GCF_002217175.1_ASM221717v1_genomic.gtf
# 3) Use sed to remove incorrect gene_id and* transcript_id fields, then change Name field to gene_name field
cat GCF_002217175.1_ASM221717v1_genomic.gtf | sed 's/transcript_id.*gene_name/gene_id/' | sed 's/Name/gene_name/' | sed 's/transcript_id.*gene_name[^;]*;//'> GCF_002217175.1_ASM221717v1_genomic.corrected.no_gene_name.gtf
# The following biotypes do not have a transcript ID so must be added with an overly complex substitution:
# 'gene_biotype "tRNA"'
# 'gene_biotype "pseudogene"' 
cat GCF_002217175.1_ASM221717v1_genomic.corrected.no_gene_name.gtf  | sed   '/gene_biotype "pseudogene";/ s/gene_id \([^;]*\)\(;.* pseudo "true";\)/gene_id \1 ; transcript_id \1\2/' > GCF_002217175.1_ASM221717v1_genomic.corrected.no_gene_name.added_transcript_ids.gtf
cat GCF_002217175.1_ASM221717v1_genomic.corrected.no_gene_name.added_transcript_ids.gtf  | sed   '/gene_biotype "tRNA";/ s/gene_id \([^;]*\)\(;.* gene_biotype "tRNA";\)/gene_id \1 ; transcript_id \1\2/'  > GCF_002217175.1_ASM221717v1_genomic.corrected.no_gene_name.added_transcript_ids_final.gtf

# Generate STAR index for Folsomia candida
~/programs/STAR/source/STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/ --genomeFastaFiles ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.fna --sjdbGTFfile GCF_002217175.1_ASM221717v1_genomic.corrected.no_gene_name.added_transcript_ids_final.gtf --genomeSAindexNbases 13

# Make reference genome compatible with QoRTs by removing extra info in FASTA headers
cat ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.fna | sed 's/\s.*$//' >  ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.QoRTs_format.fna


# This turned out to be too clunky... Attempting to rearrange fields within column 9 of the GTF
# cat  GCF_002217175.1_ASM221717v1_genomic.corrected.gtf | awk   'BEGIN{OFS="\t"; FS="\t"} {n=split($9,A,";"); print $1,$2,$3,$4,$5,$6,$7,$8,A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9],A[10,A[11],A[12],A[13],A[14],A[15]} ' | less
# Remove problematic "Name" field ### THIS IS WRONG, DON'T DO THIS
# cat GCF_002217175.1_ASM221717v1_genomic.gff | sed 's/Name=.[^;]*//' > GCF_002217175.1_ASM221717v1_genomic.no_name.gff
###  MANDATORY GTF FIELDS
# gene_id "ENSBTAG00000020601"; transcript_id "ENSBTAT00000027448"; gene_name "ZNF366";


