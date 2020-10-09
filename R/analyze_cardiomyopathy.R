# Analyze SJLIFE data
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])

library(dplyr)
library(magrittr)
library(GWASTools)
library(SNPRelate)
library(data.table)
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/ICSKAT.R")
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/clean_geno.R")
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/ICSKAT_fit_null.R")
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/make_IC_dmat.R")

output_dir <- "/rsrch3/home/biostatistics/rsun3/sjlife/output"
max_gene_size <- 5000

# load the gene information
setwd("/rsrch3/home/biostatistics/rsun3/github/LungCancerAssoc/Data")
load(file="ensembl_refgene_hg19_20180109.rda")
all_genes <- ensembl_refgene_hg19_20180109

# load outcome data
setwd("/rsrch3/home/biostatistics/rsun3/sjlife/data")
cmTimes <- fread("cardiomyopathyTimes.txt", data.table=FALSE) %>%
	filter(!is.na(U) & !is.na(V) & everChecked == TRUE)

# only use full covariates
# no black people in WGS?
cmCovar <- fread("cardiomyopathyCovar.txt", data.table=FALSE) %>%
	filter(as.character(sjlid) %in% as.character(cmTimes$sjlid)) %>%
	filter(as.character(race) %in% c("White", "Black")) %>%
	filter(!(cyclophos == 1 & is.na(cyclophosDose))) %>%
	filter(!(ifos == 1 & is.na(ifosDose))) %>%
	filter(!(chest == 1 & is.na(chestDose))) %>%
	mutate(cyclophosDose = ifelse(is.na(cyclophosDose), 0, cyclophosDose)) %>%
	mutate(ifosDose = ifelse(is.na(ifosDose), 0, ifosDose)) %>%
	mutate(chestDose = ifelse(is.na(chestDose), 0, chestDose)) %>%
	mutate(anthDose = ifelse(is.na(anthDose), 0, anthDose)) %>%	
	mutate(race = ifelse(as.character(race) == "Black", 1, 0)) %>%
	mutate(gender = ifelse(as.character(gender) == "Male", 1, 0)) %>%	
	mutate(anth = ifelse(is.na(anth), 0, 1)) %>%	
	select(sjlid, agediag, gender, cyclophos, ifos, race, chest, anth)

# merge covariates, times
cmData <- merge(cmTimes, cmCovar, by="sjlid")

# Pick a section of genes to analyze based on aID argument
genes_per_chunk <- 200
nchunks <- ceiling(table(all_genes$Chr) / genes_per_chunk)
chunk_pts <- cumsum(nchunks)
if (aID <= chunk_pts[1]) {
	chr <- 1
	start_row <- genes_per_chunk*(aID - 1) + 1
	end_row <- genes_per_chunk*aID
} else {
	chr <- min(which(chunk_pts >= aID))
	start_row <- (aID - chunk_pts[chr-1] -1) * genes_per_chunk + 1
	end_row <- (aID - chunk_pts[chr-1]) * genes_per_chunk
}
buffer <- 5000

#------------------------------------------------------------------#
# End possibly varying parameters
#------------------------------------------------------------------#

# Cut to that chromosome
gene_info <- all_genes %>% filter(Chr == chr) %>% 
	arrange(txStart) %>%
	slice(start_row:min(nrow(.), end_row))

# if the .gds does not exist, will create a new one
# warning, this will be the same size as the WGS data, could be very big
setwd("/rsrch3/home/biostatistics/rsun3/sjlife/data/")
fname <- paste0("SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr", chr, ".PASS.decomposed.sjlid_eur")
gdsfile <- paste0(fname, ".gds")
if (!file.exists(gdsfile)) {
    bed_file <- paste0(fname, ".bed")
    fam_file <- paste0(fname, ".fam")
    bim_file <- paste0(fname, ".bim")
    snpgdsBED2GDS(bed_file, fam_file, bim_file, gdsfile, family = TRUE,
                  cvt.chr = "int", cvt.snpid = "int")
}
# read in the gds file
gds <- GdsGenotypeReader(gdsfile, YchromCode=24L, XYchromCode=25L)
# matches indices with positions
pos <- getPosition(gds)
# the scan IDs for those in the outcome dataset
all_IDs <- getScanID(gds)
scan_idx <- which(all_IDs %in% as.character(cmTimes$sjlid))
scan_ids <- all_IDs[scan_idx]
# other parameters
n <- length(scan_ids)
impute <- FALSE
rareMAF <- 0.05

# results data frame
resultsDF <- data.frame(gene = gene_info$HGNC_name, start=gene_info$txStart,
                        end=gene_info$txEnd, q=NA, cleanq=NA, rareq=NA, skatp=NA,
                        burdenp=NA, complex=NA, Wskatp=NA, Wburdenp=NA, Wcomplex=NA,
                        Rskatp=NA, Rburdenp=NA, Rcomplex=NA, RWskatp=NA, RWburdenp=NA,
                        RWcomplex=NA)

# loop through genes
for (gene_it in 1:nrow(gene_info)) {

    # pick out the indices
    start_idx <- min(which(pos >= gene_info$txStart[gene_it] - buffer))
    end_idx <- max(which(pos <= gene_info$txEnd[gene_it] + buffer))

    # stop if only one SNP
    q <- end_idx - start_idx + 1
    resultsDF$q[gene_it] <- q
    if (q <= 1) (next)

    # extract genotypes
    geno_data <- getGenotypeSelection(gds, snp=(start_idx:end_idx), scan=scan_idx)
    alleleA <- getAlleleA(gds, index=(start_idx:(start_idx+q-1)))
    alleleB <- getAlleleB(gds, index=(start_idx:(start_idx+q-1)))

    # make a map file
    map_file <- data.frame(BP=pos[start_idx:end_idx], A=alleleA, B=alleleB)

    # remove non-SNPs
    simpleSNPs <- which(!(nchar(alleleA) > 1 | alleleA == "-" | nchar(alleleB) > 1 | alleleB == "-"))
    geno_data <- geno_data[simpleSNPs, ]
    map_file <- map_file[simpleSNPs, ]

    # clean the genotypes
    cleaned_output <- clean_geno(n=n, geno_data=geno_data, map_file=map_file)
    cleanG = cleaned_output$cleanG
    cleanMap <- cleaned_output$cleanMap
    if (nrow(cleanMap) < 2){
        resultsDF$cleanq[gene_it] <- nrow(cleanMap)
        next
    }

    # some will have MAF 0
    MAFs <- apply(cleanG, 2, mean) / 2
    MAFpos <- which(MAFs > 0)
    cleanG <- cleanG[, MAFpos]
    cleanMap <- cleanMap[MAFpos, ]
    MAFs <- MAFs[MAFpos]
    resultsDF$cleanq[gene_it] <- length(MAFs)
    if (length(MAFs) < 2 | length(MAFs) > max_gene_size) {next}
    # weights by beta
    weights <- dbeta(MAFs, 1, 25)
    GW <- cleanG %*% diag(weights)

    # order the genotype and outcome matrices
		cleanG <- cleanG[which(rownames(cleanG) %in% cmData$sjlid), ]		
		GW <- GW[which(rownames(GW) %in% cmData$sjlid), ]
		cmData <- cmData[match(rownames(GW), as.character(cmData$sjlid)), ]

    # make design matrices
    dmats <- make_IC_dmat(X=as.matrix(cmData %>% select(agediag, gender, cyclophos, ifos, chest, anth)), 
				lt=cmData$U, rt=cmData$V)

    # fit null model
    obs_ind <- as.numeric(cmData$V < 999)
    tpos_ind <- as.numeric(cmData$U > 0)
    null_fit <- ICSKAT_fit_null(init_beta=c(rep(0, 6), 1, 1, 0.5), lt=cmData$U, rt=cmData$V,
                                left_dmat=dmats$left_dmat,
                                right_dmat=dmats$right_dmat,
                                obs_ind=obs_ind, tpos_ind=tpos_ind)
    # get pvalue
    skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                          right_dmat=dmats$right_dmat, G=cleanG, lt=cmData$U, rt=cmData$V,
                          null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
    # get weighted p-value
    weighted_skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                                   right_dmat=dmats$right_dmat, G=GW, lt=cmData$U, rt=cmData$V,
                                   null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
    # regular skat results
    resultsDF$skatp[gene_it] <- skat_output$p_SKAT
    resultsDF$burdenp[gene_it] <- skat_output$p_burden
    resultsDF$complex[gene_it] <- skat_output$complex
    # weighted results
    resultsDF$Wskatp[gene_it] <- weighted_skat_output$p_SKAT
    resultsDF$Wburdenp[gene_it] <- weighted_skat_output$p_burden
    resultsDF$Wcomplex[gene_it] <- weighted_skat_output$complex

    # only rare SNPs
    rareSNPs <- which(MAFs <= 0.05)
    resultsDF$rareq[gene_it] <- length(rareSNPs)
    if (length(rareSNPs) < 2) {next}
    MAFs_rare <- MAFs[rareSNPs]
    Gmat_rare <- cleanG[, rareSNPs]
    weights_rare <- dbeta(MAFs_rare, 1, 25)
    GW_rare <- Gmat_rare %*% diag(weights_rare)

    # get pvalue
    rare_skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                          right_dmat=dmats$right_dmat, G=Gmat_rare, lt=cmData$U, rt=cmData$V,
                          null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
    # get weighted p-value
    RW_skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                                   right_dmat=dmats$right_dmat, G=GW_rare, lt=cmData$U, rt=cmData$V,
                                   null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)

    # rare skat results
    resultsDF$Rskatp[gene_it] <- rare_skat_output$p_SKAT
    resultsDF$Rburdenp[gene_it] <- rare_skat_output$p_burden
    resultsDF$Rcomplex[gene_it] <- rare_skat_output$complex
    # weighted results
    resultsDF$RWskatp[gene_it] <- RW_skat_output$p_SKAT
    resultsDF$RWburdenp[gene_it] <- RW_skat_output$p_burden
    resultsDF$RWcomplex[gene_it] <- RW_skat_output$complex
    cat(gene_it)
    #if (gene_it%%10 == 0) {cat(gene_it)}
}

# close the gds connection
close(gds)

# write results
setwd(output_dir)
write.table(resultsDF, paste0("CM2results", aID, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep="\t")


