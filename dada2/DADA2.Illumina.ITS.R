### Barebones commands to run DADA2 on Illumina data ###
### Looking for ITS in FASTQ name, and using ITS taxonomic assignment
 
library(dada2)
# library(msa)
library("DECIPHER")

# Define folders
path <- "/input/files/path"
filt_path <- file.path("/output/files/path", "filtered")
plots_path <- file.path("/output/files/path", "plots")
dir.create(plots_path)

# Get forward and reverse FQ files, ensure globbing pattern is correct
fnFs <- sort( list.files(path, pattern=".*ITS.*_R1.fq.gz", full.names = TRUE))
fnRs <- sort( list.files(path, pattern=".*ITS.*_R2.fq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
sample.names <- sapply(strsplit(sample.names, "_R1"), `[`, 1)

pdf(file=paste0(plots_path,"/Read_quality.pdf"), width=8.5, height=11)
plotQualityProfile(fnFs[1:length(fnFs)])
plotQualityProfile(fnRs[1:length(fnRs)])
dev.off()

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <-filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

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
# dim(seqtab)
# table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(100,500)]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seqtab.nochim)
# sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
#head(track)
write.table(track, paste0(plots_path,"/read_tracking.txt"), quote=FALSE, sep="\t")

taxa_its <- assignTaxonomy(seqtab.nochim, "~/dbs/unite/sh_general_release_dynamic_s_01.12.2017.fasta", multithread=TRUE, tryRC=TRUE)
# taxa_its <- addSpecies(taxa, "~/dbs/")

# taxa_gg <- assignTaxonomy(seqtab.nochim, "~/dbs/mothur/greengenes/gg_13_8_train_set_97.fa.gz, multithread=TRUE, tryRC=TRUE)
# taxa_gg <- addSpecies(taxa, "~/dbs/")

# taxa.print <- taxa # Removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
#library(msa)
#mult <- msa(seqs, method="ClustalW", type="dna", order="input")
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

save.image(file=paste0(filt_path,"/DADA2_Illumina_output.Rdata"))
savehistory(file = paste0(filt_path,"/DADA2_Illumina.Rhistory")

# for (thing in ls()) { message(thing); print(object.size(get(thing)), units='auto') }

rm(derepRs)
rm(derepFs)
rm(dadaFs)
rm(dadaRs)
save.image(file=paste0(filt_path,"/DADA2_Illumina_output.compressed.Rdata"))

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