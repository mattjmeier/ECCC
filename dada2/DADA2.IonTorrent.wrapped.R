### Barebones commands to run DADA2 on Ion Torrent data ###
### Argument 1 should be where your input files are
### Argument 2 should be where you want the output
### Argument 3 should glob for a list of files for DADA2 input

args<-commandArgs(TRUE)
library(dada2)
library("DECIPHER")
# library(msa)

path <- args[1]
filt_path <- file.path(args[2], "filtered")
plots_path <- file.path(args[2], "plots")
dir.create(plots_path)

paste("Output paths:",path,filt_path,plots_path,sep="\n")

fnFs <- sort(list.files(path, pattern=args[3], full.names = TRUE))
print(fnFs)
sample.names <- sapply(strsplit(basename(fnFs), "\\.fastq"), `[`, 1)
print(sample.names)

if(!file.exists(paste0(plots_path,"/Read_quality.pdf"))) {
	pdf(file=paste0(plots_path,"/Read_quality.pdf"), width=8.5, height=11)
	plotQualityProfile(fnFs[1:length(fnFs)])
	dev.off()
}

# Make list of files for filtering
filtFs <- file.path(filt_path, paste0(sample.names, "_filtered.fastq.gz"))
print(filtFs)

# Filter reads
out <-filterAndTrim(fnFs, filtFs,  truncLen=240, maxN=0, maxEE=2, truncQ=2, trimLeft=15, rm.phix=FALSE, compress=TRUE, multithread=TRUE)

# Collect actual list of filtered FASTQ for which files were created (in case some had no reads passing filter)
filtFs <- list.files(filt_path, full.names=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
pdf(file=paste0(plots_path,"/Error_profiles.pdf"), width=8.5, height=11)
plotErrors(errF, nominalQ=TRUE)
dev.off()

derepFs <- derepFastq(filtFs, verbose=TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "\\_filtered.fastq.gz"), `[`, 1)
print(sample.names)
names(derepFs) <- sample.names

dadaFs <- dada(derepFs, err=errF, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, multithread=TRUE)
seqtab <- makeSequenceTable(dadaFs)
print(dim(seqtab))
print(table(nchar(getSequences(seqtab))))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
print(dim(seqtab.nochim))
print(sum(seqtab.nochim)/sum(seqtab))

taxa_silva <- assignTaxonomy(seqtab.nochim, "~/dbs/mothur/silva/silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
taxa_silva <- addSpecies(taxa_silva, "~/dbs/mothur/silva/silva_species_assignment_v132.fa.gz")

taxa_rdp <- assignTaxonomy(seqtab.nochim, "~/dbs/mothur/rdp/rdp_train_set_16.fa.gz", multithread=TRUE, tryRC=TRUE)
taxa_rdp <- addSpecies(taxa_rdp, "~/dbs/mothur/rdp/rdp_species_assignment_16.fa.gz")

taxa_gg <- assignTaxonomy(seqtab.nochim, "~/dbs/mothur/greengenes/gg_13_8_train_set_97.fa.gz", multithread=TRUE, tryRC=TRUE)
#taxa_gg <- addSpecies(taxa_gg, "~/dbs/DOESNOTEXIST?")

taxa.print <- taxa_silva # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

## Create multisequence alignment of all unique sequences
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
# Other method:
# library(msa)
# mult <- msa(seqs, method="ClustalW", type="dna", order="input")

# Use phangorn library to create tree from multisequence alignment
library(phangorn)
phang.align <- phyDat(as(alignment,"matrix"),type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) ## Tip order is not the same as sequence order!
fit = pml(treeNJ, data=phang.align)
# negative edges lengths changed to 0
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

# Create a basic data frame for sample data. This will need to be reconfigured for further analysis.
sampledata.df <- data.frame(row.names = row.names(seqtab.nochim), number=c(1:nrow(seqtab.nochim)))
sampleOrder = row.names(sampledata.df)

# Do some preliminary phyloseq work
library(phyloseq)
# Import tree from database
x = read_tree_greengenes("~/dbs/trees/97_otus.tree")
# Create phyloseq object
# ps_silva <- phyloseq(tax_table(taxa_silva), sample_data(sampledata.df), otu_table(seqtab.nochim, taxa_are_rows=FALSE), phy_tree(fitGTR$tree))
# ps_gg <- phyloseq(tax_table(taxa_gg), sample_data(sampledata.df), otu_table(seqtab.nochim, taxa_are_rows=FALSE), phy_tree(x))

# Save current state
save.image(file=paste0(args[2],"/R_output.IonTorrent.FULL.RData"))
# savehistory(file = paste0(path,"/DADA2.IonTorrent.Rhistory"))

# Save lightweight version
dadaFs <- NULL
derepFs <- NULL
x <- NULL

save.image(file=paste0(args[2],"/R_output.IonTorrent.light.RData"))

##### COMPLETE IN WINDOWS #####
