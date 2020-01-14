# Analyze SJLIFE data
library(dplyr)
library(magrittr)
library(GWASTools)
library(SNPRelate)

# load outcome data
setwd('/users/ryansun/desktop/skat/ghd')
ghd_dat <- read.csv("time_678_999.csv")

# Pick a chromosome
chr <- 14
buffer <- 5000
setwd('/users/ryansun/desktop/skat/chr14')

# load the gene data on the selected chromosome
library(LungCancerAssoc)
data("ensembl_refgene_hg19_20180109")
gene_info <- ensembl_refgene_hg19_20180109 %>% filter(Chr == chr) %>% arrange(txStart)

# if the .gds does not exist, will create a new one
# warning, this will be the same size as the WGS data, could be very big
fname <- "SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr14.PASS.decomposed.sjlid_eur"
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
scan_idx <- which(all_IDs %in% as.character(ghd_dat$SJLIFEID))
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
    if (length(MAFs) < 2) {next}
    # weights by beta
    weights <- dbeta(MAFs, 1, 25)
    GW <- cleanG %*% diag(weights)

    # order the genotype and outcome matrices
    temp_ghd <- ghd_dat %>% filter(SJLIFEID %in% scan_ids)
    temp_ghd <- temp_ghd[match(rownames(cleanG), as.character(temp_ghd$SJLIFEID)), ]

    # make design matrices
    dmats <- make_IC_dmat(X=NULL, lt=temp_ghd$finalU, rt=temp_ghd$finalV)

    # fit null model
    obs_ind <- as.numeric(temp_ghd$finalV < 999)
    tpos_ind <- as.numeric(temp_ghd$finalU > 0)
    null_fit <- ICSKAT_fit_null(init_beta=c(1, 1, 0.5), lt=temp_ghd$finalU, rt=temp_ghd$finalV,
                                left_dmat=dmats$left_dmat,
                                right_dmat=dmats$right_dmat,
                                obs_ind=obs_ind, tpos_ind=tpos_ind)
    # get pvalue
    skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                          right_dmat=dmats$right_dmat, G=cleanG, lt=temp_ghd$finalU, rt=temp_ghd$finalV,
                          null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
    # get weighted p-value
    weighted_skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                                   right_dmat=dmats$right_dmat, G=GW, lt=temp_ghd$finalU, rt=temp_ghd$finalV,
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
                          right_dmat=dmats$right_dmat, G=Gmat_rare, lt=temp_ghd$finalU, rt=temp_ghd$finalV,
                          null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
    # get weighted p-value
    RW_skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                                   right_dmat=dmats$right_dmat, G=GW_rare, lt=temp_ghd$finalU, rt=temp_ghd$finalV,
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



#' Clean gds genotype function, called by geno_sjlife_gds.R
#'
#' @param n only keep SNPs with at least n non-missing entries
#' @param geno_data genotype matrix from the getGenotype() function
#'
#' @return An n*m matrix where m depends on the number of SNPs that passed QC.

clean_geno <- function(n, geno_data, map_file) {

    # loop through SNPs, remove missing and add two halves
    Gmat <- c()
    removed <- c()
    flipped <- c()
    for (snp_it in 1:nrow(geno_data)) {
        # remove missing
        temp_snp <- geno_data[snp_it, ]
        temp_snp <- temp_snp[which(!is.na(temp_snp))]
        if (length(temp_snp) < n) {
            removed <- c(removed, snp_it)
            next
        }

        # flip to minor allele
        if (mean(temp_snp) > 1) {
            temp_snp <- 2 - temp_snp
            flipped <- c(flipped, snp_it)
        }

        # add SNP to data matrix
        Gmat <- cbind(Gmat, temp_snp[1:n])
    }

    # flip before remove in map file
    temp_new_A <- map_file$A
    temp_new_A[flipped] <- map_file$B[flipped]
    map_file$B[flipped] <- map_file$A[flipped]
    map_file$A <- temp_new_A
    # remove
    map_file <- map_file[-removed, ]

    return(list(cleanG=Gmat, cleanMap = map_file))
}


#' Clean gds genotype function, called by geno_sjlife_gds.R.
#' Imputes missing data.
#'
#' @param n only keep SNPs with at least n non-missing entries
#' @param geno_data genotype matrix from the getGenotype() function
#'
#' @return An n*m matrix where m depends on the number of SNPs that passed QC.

clean_geno_impute <- function(n, geno_data) {

    # impute each column of geno_data
    Gmat <- t(geno_data)
    for (snp_it in 1:ncol(Gmat)) {

        # identify missing first
        miss_idx <- which(is.na(Gmat[, snp_it]))

        # flip to minor allele
        temp_MAF <- mean(Gmat[, snp_it], na.rm=TRUE)
        if (temp_MAF > 1) {
            Gmat[, snp_it] <- 2 - Gmat[, snp_it]
            temp_MAF <- mean(Gmat[, snp_it], na.rm=TRUE)
        }

        # impute
        if (length(miss_idx) == 0) {next}
        Gmat[miss_idx, snp_it] <- rbinom(n=length(miss_idx), size=2, prob=temp_MAF)

        # flip to minor allele again after impute
        if (mean(Gmat[, snp_it]) > 1) {Gmat[, snp_it] <- 2 - Gmat[, snp_it]}
    }

    return(Gmat)
}
