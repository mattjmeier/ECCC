### Barebones commands to run DADA2 on Illumina data ###

library(dada2)
# library(msa)
library("DECIPHER")

# Define folders
path <- "path/to/fastqs"
filt_path <- file.path("/output/of/dada2", "filtered")
plots_path <- file.path("/output/of/dada2", "plots")
dir.create(plots_path)

# Get forward and reverse FQ files, ensure globbing pattern is correct for your case!
fnFs <- sort( list.files(path, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort( list.files(path, pattern="R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
sample.names <- sapply(strsplit(sample.names, "_L001_R1_001"), `[`, 1)

pdf(file=paste0(plots_path,"/Read_quality.pdf"), width=8.5, height=11)
plotQualityProfile(fnFs[1:length(fnFs)])
plotQualityProfile(fnRs[1:length(fnRs)])
dev.off()

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <-filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(241,231), trimLeft=1, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf(file=paste0(plots_path,"/Error_profiles.pdf"), width=8.5, height=11)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
write.table(dim(seqtab), paste0(plots_path,"/dim_seqtab.txt"), quote=FALSE, sep="\t")
write.table(table(nchar(getSequences(seqtab))), paste0(plots_path,"/sequence_lengths_after_merging.txt"), quote=FALSE, sep="\t")
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(432,447)]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
write.table(dim(seqtab.nochim), paste0(plots_path,"/dim_seqtab.nochim.txt"), quote=FALSE, sep="\t")
write.table(sum(seqtab.nochim)/sum(seqtab), paste0(plots_path,"/proportion_of_bimera.txt"), quote=FALSE, sep="\t")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
#head(track)
write.table(track, paste0(plots_path,"/read_tracking.txt"), quote=FALSE, sep="\t")

taxa_silva <- assignTaxonomy(seqtab.nochim, "~/dbs/mothur/silva/silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
taxa_silva <- addSpecies(taxa_silva, "~/dbs/mothur/silva/silva_species_assignment_v132.fa.gz")

taxa_rdp <- assignTaxonomy(seqtab.nochim, "~/dbs/mothur/rdp/rdp_train_set_16.fa.gz", multithread=TRUE)
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

save.image(file=paste0(filt_path,"/DADA2_Illumina_output.Mouse_AgNP.Date.Rdata"))
savehistory(file = paste0(filt_path,"/DADA2_Illumina.Mouse_AgNP.Date.Rhistory"))

##### COMPLETE IN WINDOWS #####

library(phangorn)
phang.align <- phyDat(as(alignment,"matrix"),type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) ## Tip order is not the same as sequence order!
fit = pml(treeNJ, data=phang.align)
# negative edges lengths changed to 0
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

detatch("package:phangorn", unload=TRUE)