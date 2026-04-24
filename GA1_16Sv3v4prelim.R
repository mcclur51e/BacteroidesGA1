##########################################################################################################
##### DESCRIPTION                                                                                    #####
##### Code for preliminary processing of .fastq files from Illumina sequencing. Processing includes: #####
##### 1. standard DADA2 pipeline                                                                     #####
##### 2. assigning taxonomy using a Bacteroides-specific database                                    #####
##### 3. passing data into phyloseq object                                                           #####
##### 4. correcting a few mis-identified important taxa                                              #####
##### 5. adding columns to mapping file based on calculations from processed sequencing data         #####
##### 6. identify and remove contaminants using decontam                                             #####
##### Use this file to prepare initial raw phyloseq object for all future downstream analysis.       #####
###########################################################################################################################################################
##### This code requires as input:                                                                                                                    #####
##### 1. A working repository (default "~/Desktop/Wannigan/")                                                                                   #####
##### 2. Folder within working repository containing all .fastq files (default "~/Desktop/Wannigan/raw/16S/")                                   #####
##### 3. .csv file containing sample metadata listed with unique sampleIDs in first column (default "~/Desktop/Wannigan/raw/table_map.csv")     #####
##### 4. .fa file containing reference taxonomy database (default "~/Desktop/Wannigan/raw/dbBacteroides.fa")                                     #####                                   
###########################################################################################################################################################

#####################################################################
########## Processing Set −up #######################################
#####################################################################
# R version 4.5.2 (2025-10-31) -- "[Not] Part in a Rumble"
##### All code written by Emily Ann McClure github.com/mcclur51e   #####
##### No AI was used at any stage in writing or editing this code. #####

########## Call libraries for use ##########
# Packages from CRAN
library("ggplot2") # version 4.0.0
library("data.table") # version 1.17.8
library("dplyr") # version 1.1.4

# Packages from Bioconductor
library("phyloseq") # version 1.52.0
library("decontam") # version 1.28.0 # identify contaminant ASVs 
library("dada2") # version 1.36.0

################################################################################
########## If running by line (i.e. RStudio) use the following block ###########
################################################################################
path <- "~/Desktop/Wannigans/Wannigan_Casey/" # change this to directory where you will be working
setwd(path) # set working directory
dir.create("outputPrelim") # create directory for output files to go
dir.create("plots") # create directory for plots to go

################################################################################
##### ASV assignment with DADA2 ################################################
################################################################################
fnFs <- sort(list.files(paste0(path,"raw/16S/"), pattern="R1_001", full.names = TRUE))
fnRs <- sort(list.files(paste0(path,"raw/16S/"), pattern="R2_001", full.names = TRUE))

if(length(fnFs) != length(fnRs)) stop("At least one sample is unpaired. Please check forward and reverse reads are present for all samples")

sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
#sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1) # to check if reverse reads have same sample names as forward
#write.csv(sample.names,"outputPrelim/sampleNames.csv") # print list of sample names to check if they match map file
table_map <- data.frame(read.csv("raw/table_map.csv", header = TRUE, row.names = 1, check.names=FALSE)) # reads csv file into data.frame with row names in column 1

save(fnFs,file=("outputPrelim/output_fnFs.RData")) # Save in a .RData file 
save(fnRs,file=("outputPrelim/output_fnRs.RData")) # Save in a .RData file 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
save(filtFs,file=("outputPrelim/output_filtFs.RData")) # Save in a .RData file
save(filtRs,file=("outputPrelim/output_filtRs.RData")) # Save in a .RData file

filtered <- filterAndTrim(fnFs, filtFs, 
                          fnRs, filtRs,
                          trimLeft=13, trimRight=17,
                          maxLen=530, minLen = 150,
                          maxN=0, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, verbose=TRUE, multithread=TRUE) # filtering values set here are specific to this study. Modify as appropriate when sequencing other regions.
save(filtered,file=("outputPrelim/output_filtered.RData")) # Save .RData file

### Subset filtFs and filtRs to include files with > 0 reads (some files may have been emptied during error analysis) 
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]
save(filtFs,file=("outputPrelim/output_filtFs.RData")) # Save .RData file
save(filtRs,file=("outputPrelim/output_filtRs.RData")) # Save .RData file 

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
save(errF,file=("outputPrelim/output_errF.RData")) # Save .RData file
save(errR,file=("outputPrelim/output_errR.RData")) # Save .RData file
#plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
save(dadaFs,file=("outputPrelim/output_dadaFs.RData")) # Save in a .RData file
save(dadaRs,file=("outputPrelim/output_dadaRs.RData")) # Save in a .RData file
save(mergers,file=("outputPrelim/output_mergers.RData")) # Save in a .RData file

seqtab <- makeSequenceTable(mergers)
# table(nchar(getSequences(seqtab))) # use to check distribution of sequence lengths
seqtab.all <- seqtab[,nchar(colnames(seqtab)) %in% 410:450] # these cut-off values have been chosen specific to this study. Modify as appropriate when sequencing other regions.
#seqtab.all <- seqtab[,nchar(colnames(seqtab)) %in% val.minLength:val.maxLength]
seqtab.noBim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE)

getNreads <- function(x) sum(getUniques(x))
track <- cbind(filtered, sapply(dadaFs, getNreads), sapply(dadaRs, getNreads), sapply(mergers, getNreads), rowSums(seqtab), rowSums(seqtab.noBim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "tabled", "nonBim")
track

#####################################################################
########## Transfer data into Phyloseq ##############################
##### This preliminary Phyloseq object will be used to trim #########
##### singletons from the dataset before assigning taxonomy #########
#####################################################################

OTUraw = otu_table(seqtab.noBim, taxa_are_rows=FALSE) # assigns ASV table from DADA output
MAP = sample_data(table_map) # reads csv file into data.frame with row names in column 1
phyW = phyloseq(OTUraw, MAP) # prepare phyloseq object
save(phyW, file=("outputPrelim/physeq_raw.RData"))
save(OTUraw,file=("outputPrelim/output_otuRaw.RData")) # Save the phyloseq data object in a .RData file 

taxRaw <- assignTaxonomy(OTUraw, "raw/dbBacteroides.fa", multithread=TRUE) # assign taxonomy to trimmed dataset
save(taxRaw,file=("outputPrelim/output_taxRaw.RData")) # Save the phyloseq data object in a .RData file 

#####################################################################
########## Transfer data into Phyloseq ##############################
#####################################################################
#load("outputPrelim/output_otuRaw.RData")
#load("outputPrelim/output_taxRaw.RData")
#table_map <- data.frame(read.csv("raw/table_map.csv", header = TRUE, row.names = 1, check.names=FALSE)) # reads csv file into data.frame with row names in column 1

OTU = otu_table(OTUraw, taxa_are_rows=FALSE) # assigns ASV table from trimmed DADA output
TAX = tax_table(taxRaw) # assigns taxonomy table
MAP = sample_data(table_map) # assigns metadata table
physeq = phyloseq(OTU, TAX, MAP) # prepare phyloseq object

### add data to taxonomy table ###
#table_spp <- data.frame(read.csv("~/Desktop/Wannigan_METRC/raw/table_SppDescription.csv", header = TRUE, check.names=FALSE)) # reads csv file into data.frame with row names in column 1
bind.spp <- data.frame(cbind(tax_table(physeq), paste(tax_table(physeq)[,c("Genus")],tax_table(physeq)[,c("Species")]))) # add column listing ASVs numerically
bind.asv <- cbind(bind.spp, paste0("ASV", seq(ntaxa(physeq)))) # add column listing ASVs numerically
bind.seq <- cbind(bind.asv, row.names(tax_table(physeq))) # add column listing ASV sequences
bind.o <- cbind(bind.seq, seq(ntaxa(physeq))) # adding a column for sorting table later
bind.length <- cbind(bind.o, nchar(bind.o[,10])) # add a column listing length of sequences
bind.count <- cbind(bind.length, taxa_sums(physeq)) # add a column listing total reads for ASV
bind.prev <- cbind(bind.count, rowSums(t(otu_table(physeq)) != 0)) # add a column listing total number of samples in which ASV appears
bind.prevSample <- cbind(bind.prev, rowSums(t(otu_table(subset_samples(physeq, Control%in%c("sample","positive")))) != 0)) # add a column listing total number of samples (excluding controls) in which ASV appears
bind.CP <- cbind(bind.prevSample, as.numeric(bind.prevSample[,13]) / as.numeric(bind.prevSample[,14])) # calculate ~average reads/sample
colnames(bind.CP) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain","Spp","ASV",
                       "Sequence","Sort","Length","Count","Prevalence","PrevSamples","CP") # rename columns
rownames(bind.CP) <- bind.CP$ASV # reassign ASV names
TAX2 = tax_table(as.matrix(bind.CP)) # define new taxonomy table
taxa_names(physeq) <- tax_table(TAX2)[,"ASV"]
physeqR1 <- phyloseq(otu_table(physeq), TAX2, MAP) # prepare phyloseq object modified with ASV numbers

#######################################
##### add columns to mapping data #####
#######################################
df.map <- as.data.frame(sample_data(physeqR1))
df.map$sampleID <- rownames(df.map) # sample names
df.map$LibrarySize <- sample_sums(physeqR1) # total reads count per sample
df.map <- cbind(df.map, estimate_richness(physeqR1, split = TRUE, measures = NULL))
df.map$mutant <- with(df.map, ifelse(genotype%in%c("WT1","WT2"), "WT", "dtssBC")) 
physeqR = phyloseq(tax_table(physeqR1),otu_table(physeqR1),sample_data(df.map)) # return map to phyloseq object

#############################################################
########## Identify and remove likely contaminants ##########
#############################################################
### Calculate ASVs to remove based on counts in the dataset
phyR.neg2 <- subset_samples(physeqR, Control%in%c("negExtract","negSeq","negative")) # subset data to negative controls
phyR.neg <- subset_taxa(physeqR, taxa_sums(phyR.neg2) > 1 & Kingdom=="Bacteria") # subset to taxa present in negative controls
df.neg <- data.frame(otu_table(phyR.neg)) # otu table for negative controls
phyR.taxTrim <- subset_taxa(physeqR, taxa_sums(physeqR) > quantile(as.matrix(df.neg[df.neg>0]), na.rm = T, probs = c(0.7))) # keep taxa present at counts above the 70% quantile from negative controls (i.e. have to have more reads than most ASVs from negative control samples)

### identify contaminants using decontam
ps <- physeqR
sample_data(ps)$is.neg <- sample_data(ps)$Control == "negative" # define negative samples
sample_data(ps)$QubitConc <- as.numeric(sample_data(ps)$QubitConc) + 0.001 # convert values to numeric and >0
contamdf.either <- isContaminant(ps, method="either", conc="QubitConc", neg="is.neg", threshold=0.5) # determine contaminants
#table(contamdf.either$contaminant)
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control%in%c("negative","negEtract","negSeq"), ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control=="sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.either$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

phyR.noContam <- subset_taxa(phyR.taxTrim, !ASV%in%rownames(subset(df.pa, contaminant=="TRUE")))

########################################################################################
##### save phyloseq object and its components for easy access in future processing #####
########################################################################################
save(phyR.noContam,file=("outputPrelim/phyRnoContam.RData")) # Save the phyloseq data object in a .RData file 
write.csv(tax_table(phyR.noContam),"outputPrelim/table_tax.csv") # Save taxonomy table as .csv
write.csv(otu_table(phyR.noContam),"outputPrelim/table_otu.csv") # Save ASV table as .csv
write.csv(data.frame(sample_data(phyR.noContam)),"outputPrelim/table_mapModified.csv") # Save modified sample data as .csv
ps.noContam <- psmelt(phyR.noContam)
save(ps.noContam, file="outputPrelim/psNoContam.csv") # Save total data as .csv

##########################################################################################################
##### You have now successfully imported all the raw data into R and performed preliminary clean-up. #####
##### You are ready to proceed with further analysis                                                 #####
##########################################################################################################

##### All code written by Emily Ann McClure github.com/mcclur51e   #####
##### No AI was used at any stage in writing or editing this code. #####