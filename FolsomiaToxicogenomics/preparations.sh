# Make GTF file compatible with STAR
#~/programs/cufflinks-2.2.1.Linux_x86_64/gffread -F -T GCF_002217175.1_ASM221717v1_genomic.gff -o GCF_002217175.1_ASM221717v1_genomic.full.gtf
~/programs/cufflinks-2.2.1.Linux_x86_64/gffread -T GCF_002217175.1_ASM221717v1_genomic.gff -o GCF_002217175.1_ASM221717v1_genomic.gtf

# Generate STAR index for Folsomia candida
# /data/results/home/mmeier/programs/STAR/source/STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/ --genomeFastaFiles ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.fna --sjdbGTFfile ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.full.gtf --genomeSAindexNbases 13
/data/results/home/mmeier/programs/STAR/source/STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/ --genomeFastaFiles ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.fna --sjdbGTFfile ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.gtf --genomeSAindexNbases 13

# Make reference genome compatible with QoRTs by removing extra info in FASTA headers
cat ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.fna | sed 's/\s.*$//' >  ~/dbs/folsomia/GCF_002217175.1_ASM221717v1/GCF_002217175.1_ASM221717v1_genomic.QoRTs_format.fna

