### Barebones commands to run DADA2 on Ion Torrent data ###
 
library(dada2)
library("DECIPHER")
# library(msa)
path <- "/input/files/path"
filt_path <- file.path("/output/files/path", "filtered")
plots_path <- file.path("/output/files/path", "plots")
dir.create(plots_path)

fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)

pdf(file=paste0(plots_path,"/Read_quality.pdf"), width=8.5, height=11)
plotQualityProfile(fnFs[1:length(fnFs)])
dev.off()

# filtFs <- file.path(filt_path, paste0(sample.names, "_filt.fastq.gz"))
filtFs <- list.files(filt_path, full.names=TRUE)

out <-filterAndTrim(fnFs, filtFs,  truncLen=275, maxN=0, maxEE=c(2,2), truncQ=2, trimLeft=15, rm.phix=FALSE, compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
pdf(file=paste0(plots_path,"/Error_profiles.pdf"), width=8.5, height=11)
plotErrors(errF, nominalQ=TRUE)
dev.off()

derepFs <- derepFastq(filtFs, verbose=TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "\\."), `[`, 1)
names(derepFs) <- sample.names

dadaFs <- dada(derepFs, err=errF, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, multithread=TRUE)
seqtab <- makeSequenceTable(dadaFs)
# dim(seqtab)
# table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seqtab.nochim)
# sum(seqtab.nochim)/sum(seqtab)

taxa_silva <- assignTaxonomy(seqtab.nochim, "~/dbs/mothur/silva/silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
taxa_silva <- addSpecies(taxa_silva, "~/dbs/mothur/silva/silva_species_assignment_v132.fa.gz")

taxa_rdp <- assignTaxonomy(seqtab.nochim, "~/dbs/mothur/rdp/rdp_train_set_16.fa.gz", multithread=TRUE, tryRC=TRUE)
taxa_rdp <- addSpecies(taxa_rdp, "~/dbs/mothur/rdp/rdp_species_assignment_16.fa.gz")

# taxa_gg <- assignTaxonomy(seqtab.nochim, "~/dbs/mothur/greengenes/gg_13_8_train_set_97.fa.gz, multithread=TRUE, tryRC=TRUE)
# taxa_gg <- addSpecies(taxa, "~/dbs/")

# taxa.print <- taxa # Removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
# mult <- msa(seqs, method="ClustalW", type="dna", order="input")
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

save.image(file=paste0(filt_path,"/R_output.IonTorrent.Rdata"))
savehistory(file = paste0(filt_path,"/DADA2.IonTorrent.Rhistory"))

##### COMPLETE IN WINDOWS #####
